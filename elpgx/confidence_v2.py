"""
Enhanced Confidence Score Calculator for PGx Pipeline v2.

Improvements over v1:
1. Coverage-aware: Checks depth at ALL defining positions, not just variant calls
2. Completeness-aware: Properly handles reference calls
3. CNV-aware: Incorporates CNV detection confidence for applicable genes
4. Gene-specific: Different weights for different gene complexity levels

Formula:
  confidence = w1 * coverage_score + w2 * quality_score + w3 * completeness_score + w4 * cnv_score

Where weights depend on gene type:
  - Simple genes (VKORC1, IFNL3): coverage=0.50, quality=0.35, completeness=0.15, cnv=0.00
  - Complex genes (CYP2C19, CYP2C9): coverage=0.40, quality=0.30, completeness=0.30, cnv=0.00
  - CNV genes (CYP2D6, CYP2B6): coverage=0.30, quality=0.20, completeness=0.25, cnv=0.25
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional

from elpgx.models import ConfidenceFlag, ConfidenceMetrics, Variant, VariantDepthStats


class GeneComplexity(str, Enum):
    """Gene complexity levels affecting confidence calculation."""
    SIMPLE = "simple"          # Single variant genes (VKORC1, IFNL3)
    MODERATE = "moderate"      # Multiple defining variants (CYP2C19, CYP2C9)
    COMPLEX = "complex"        # CNV-prone genes (CYP2D6, CYP2B6)
    HLA = "hla"               # HLA genes (special handling)


# Gene to complexity mapping
GENE_COMPLEXITY = {
    # Simple - single key variant
    "VKORC1": GeneComplexity.SIMPLE,
    "IFNL3": GeneComplexity.SIMPLE,
    "SLCO1B1": GeneComplexity.SIMPLE,

    # Moderate - multiple defining variants
    "CYP2C19": GeneComplexity.MODERATE,
    "CYP2C9": GeneComplexity.MODERATE,
    "CYP3A5": GeneComplexity.MODERATE,
    "TPMT": GeneComplexity.MODERATE,
    "NUDT15": GeneComplexity.MODERATE,
    "DPYD": GeneComplexity.MODERATE,
    "UGT1A1": GeneComplexity.MODERATE,
    "NAT2": GeneComplexity.MODERATE,

    # Complex - CNV-prone
    "CYP2D6": GeneComplexity.COMPLEX,
    "CYP2B6": GeneComplexity.COMPLEX,
    "G6PD": GeneComplexity.COMPLEX,

    # HLA
    "HLA-A": GeneComplexity.HLA,
    "HLA-B": GeneComplexity.HLA,
    "CACNA1S": GeneComplexity.MODERATE,
    "RYR1": GeneComplexity.MODERATE,
}


@dataclass
class ConfidenceWeights:
    """Weights for confidence components based on gene complexity."""
    coverage: float = 0.40
    quality: float = 0.30
    completeness: float = 0.30
    cnv: float = 0.00

    @classmethod
    def for_gene(cls, gene: str) -> "ConfidenceWeights":
        """Get appropriate weights for a gene."""
        complexity = GENE_COMPLEXITY.get(gene, GeneComplexity.MODERATE)

        if complexity == GeneComplexity.SIMPLE:
            return cls(coverage=0.50, quality=0.35, completeness=0.15, cnv=0.00)
        elif complexity == GeneComplexity.MODERATE:
            return cls(coverage=0.40, quality=0.30, completeness=0.30, cnv=0.00)
        elif complexity == GeneComplexity.COMPLEX:
            return cls(coverage=0.30, quality=0.20, completeness=0.25, cnv=0.25)
        elif complexity == GeneComplexity.HLA:
            # HLA typing has different confidence model
            return cls(coverage=0.30, quality=0.30, completeness=0.40, cnv=0.00)
        else:
            return cls()


@dataclass
class ConfidenceConfigV2:
    """Enhanced configuration for confidence calculation."""

    # Thresholds
    target_depth: int = 30
    min_depth_threshold: int = 10
    max_quality: float = 100.0
    min_quality_threshold: float = 20.0

    # CNV thresholds
    cnv_min_ratio: float = 0.7   # Min CN ratio for confident deletion
    cnv_max_ratio: float = 1.3   # Max CN ratio for confident single copy

    # Position coverage requirements
    min_position_coverage: int = 5  # Minimum depth to consider position "covered"


@dataclass
class CoverageAtPosition:
    """Coverage information at a specific genomic position."""
    chrom: str
    pos: int
    depth: int
    is_variant: bool = False
    variant_quality: float = 0.0


@dataclass
class CNVEvidence:
    """Copy number variation evidence for a gene."""
    gene: str
    estimated_copy_number: float
    confidence: float  # 0-1 confidence in CNV call
    method: str  # "read_depth", "paired_end", "split_read", "none"
    is_deletion: bool = False
    is_duplication: bool = False

    @classmethod
    def no_cnv_data(cls, gene: str) -> "CNVEvidence":
        """Create CNV evidence when no CNV analysis was performed."""
        return cls(
            gene=gene,
            estimated_copy_number=2.0,  # Assume diploid
            confidence=0.5,  # Uncertain
            method="none",
        )


class ConfidenceFlag_V2(str, Enum):
    """Extended quality flags for v2 confidence."""
    # Original flags
    LOW_DEPTH = "LOW_DEPTH"
    LOW_QUALITY = "LOW_QUALITY"
    INCOMPLETE_ALLELE = "INCOMPLETE_ALLELE"
    AMBIGUOUS_CALL = "AMBIGUOUS_CALL"
    NO_CALL = "NO_CALL"

    # New flags
    POSITION_NOT_COVERED = "POSITION_NOT_COVERED"  # Defining position has no coverage
    CNV_DETECTED = "CNV_DETECTED"                   # CNV may affect interpretation
    CNV_UNCERTAIN = "CNV_UNCERTAIN"                 # CNV status unknown
    REFERENCE_ASSUMED = "REFERENCE_ASSUMED"         # Reference assumed due to no variant call
    PHASE_UNKNOWN = "PHASE_UNKNOWN"                 # Haplotype phasing uncertain
    HYBRID_ALLELE = "HYBRID_ALLELE"                # Possible hybrid/recombinant allele


def calculate_coverage_score(
    positions_coverage: list[CoverageAtPosition],
    config: ConfidenceConfigV2,
) -> tuple[float, float, int, int]:
    """
    Calculate coverage score across ALL defining positions.

    Unlike v1 which only checked variant positions, this checks
    coverage at every defining position (variant or reference).

    Returns:
        Tuple of (coverage_score, mean_depth, min_depth, uncovered_positions)
    """
    if not positions_coverage:
        return 0.0, 0.0, 0, 0

    depths = [p.depth for p in positions_coverage]
    mean_depth = sum(depths) / len(depths)
    min_depth = min(depths)

    # Count uncovered positions
    uncovered = sum(1 for p in positions_coverage if p.depth < config.min_position_coverage)

    # Score based on:
    # 1. Mean depth relative to target
    # 2. Penalty for uncovered positions
    depth_ratio = min(1.0, mean_depth / config.target_depth)
    coverage_penalty = uncovered / len(positions_coverage)

    coverage_score = depth_ratio * (1.0 - coverage_penalty)

    return coverage_score, mean_depth, min_depth, uncovered


def calculate_quality_score_v2(
    positions_coverage: list[CoverageAtPosition],
    config: ConfidenceConfigV2,
) -> float:
    """
    Calculate quality score from variant positions.

    For reference positions (no variant), quality is inferred from depth.
    """
    if not positions_coverage:
        return 0.0

    scores = []
    for p in positions_coverage:
        if p.is_variant:
            # Use actual variant quality
            score = min(1.0, p.variant_quality / config.max_quality)
        else:
            # For reference positions, infer quality from depth
            # High depth = high confidence in reference call
            score = min(1.0, p.depth / config.target_depth)
        scores.append(score)

    return sum(scores) / len(scores)


def calculate_completeness_score_v2(
    positions_coverage: list[CoverageAtPosition],
    total_defining_positions: int,
    config: ConfidenceConfigV2,
) -> float:
    """
    Calculate completeness score with proper handling of reference calls.

    A position is "complete" if:
    - It has a variant call, OR
    - It has sufficient coverage to confidently call reference
    """
    if total_defining_positions == 0:
        return 1.0

    # Count positions that are confidently genotyped
    confident_positions = 0
    for p in positions_coverage:
        if p.is_variant:
            # Variant was called
            confident_positions += 1
        elif p.depth >= config.min_position_coverage:
            # Reference call with sufficient coverage
            confident_positions += 1
        # else: position not confidently genotyped

    return confident_positions / total_defining_positions


def calculate_cnv_score(
    cnv_evidence: CNVEvidence,
    config: ConfidenceConfigV2,
) -> float:
    """
    Calculate CNV confidence score.

    High score = confident in copy number state
    Low score = uncertain about CNVs
    """
    if cnv_evidence.method == "none":
        # No CNV analysis performed
        return 0.5  # Neutral/uncertain

    # Check if copy number is in expected diploid range
    cn = cnv_evidence.estimated_copy_number

    if config.cnv_min_ratio * 2 <= cn <= config.cnv_max_ratio * 2:
        # Looks like normal diploid (CN ~2)
        return cnv_evidence.confidence
    else:
        # Abnormal copy number detected
        # Score reflects confidence in that detection
        return cnv_evidence.confidence


def determine_flags_v2(
    positions_coverage: list[CoverageAtPosition],
    completeness_score: float,
    call_quality: str,
    cnv_evidence: Optional[CNVEvidence],
    config: ConfidenceConfigV2,
) -> list[ConfidenceFlag]:
    """
    Determine quality flags with enhanced checks.
    """
    flags: list[ConfidenceFlag] = []

    # Check for uncovered positions
    uncovered = [p for p in positions_coverage if p.depth < config.min_position_coverage]
    if uncovered:
        flags.append(ConfidenceFlag.LOW_DEPTH)

    # Check variant quality
    low_qual_variants = [
        p for p in positions_coverage
        if p.is_variant and p.variant_quality < config.min_quality_threshold
    ]
    if low_qual_variants:
        flags.append(ConfidenceFlag.LOW_QUALITY)

    # Check completeness
    if completeness_score < 1.0:
        flags.append(ConfidenceFlag.INCOMPLETE_ALLELE)

    # Check call quality
    if call_quality == "ambiguous":
        flags.append(ConfidenceFlag.AMBIGUOUS_CALL)
    elif call_quality == "no_call":
        flags.append(ConfidenceFlag.NO_CALL)

    # Note: Extended flags (POSITION_NOT_COVERED, CNV_*, etc.)
    # would use ConfidenceFlag_V2 enum

    return flags


def calculate_confidence_v2(
    gene: str,
    positions_coverage: list[CoverageAtPosition],
    total_defining_positions: int,
    call_quality: str = "definitive",
    cnv_evidence: Optional[CNVEvidence] = None,
    config: Optional[ConfidenceConfigV2] = None,
) -> ConfidenceMetrics:
    """
    Enhanced confidence calculation with gene-specific weights.

    Args:
        gene: Gene name (e.g., "CYP2D6")
        positions_coverage: Coverage at all defining positions
        total_defining_positions: Total number of defining positions
        call_quality: Quality from diplotype caller
        cnv_evidence: Optional CNV analysis results
        config: Configuration options

    Returns:
        ConfidenceMetrics with score and flags
    """
    if config is None:
        config = ConfidenceConfigV2()

    # Get gene-specific weights
    weights = ConfidenceWeights.for_gene(gene)

    # Calculate component scores
    coverage_score, mean_depth, min_depth, uncovered = calculate_coverage_score(
        positions_coverage, config
    )

    quality_score = calculate_quality_score_v2(positions_coverage, config)

    completeness_score = calculate_completeness_score_v2(
        positions_coverage, total_defining_positions, config
    )

    # CNV score (only for complex genes)
    if cnv_evidence is None:
        cnv_evidence = CNVEvidence.no_cnv_data(gene)
    cnv_score = calculate_cnv_score(cnv_evidence, config)

    # Calculate weighted composite score
    composite_score = (
        weights.coverage * coverage_score +
        weights.quality * quality_score +
        weights.completeness * completeness_score +
        weights.cnv * cnv_score
    )
    composite_score = round(composite_score, 2)

    # Determine flags
    flags = determine_flags_v2(
        positions_coverage, completeness_score, call_quality, cnv_evidence, config
    )

    # Build depth stats
    depth_stats = None
    if positions_coverage:
        depth_stats = VariantDepthStats(
            mean_depth=round(mean_depth, 1),
            min_depth=min_depth,
        )

    return ConfidenceMetrics(
        score=composite_score,
        flags=flags,
        variant_depth_stats=depth_stats,
    )


# Utility function to get coverage at positions from BAM
def get_coverage_at_positions(
    bam_path: str,
    positions: list[tuple[str, int]],  # List of (chrom, pos)
) -> list[CoverageAtPosition]:
    """
    Get read depth at specified genomic positions from BAM file.

    This would use pysam to query the BAM file.
    """
    try:
        import pysam
    except ImportError:
        # Return empty if pysam not available
        return []

    coverage = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for chrom, pos in positions:
            # Get pileup at position
            depth = 0
            for pileup_column in bam.pileup(chrom, pos - 1, pos, truncate=True):
                if pileup_column.pos == pos - 1:  # pysam is 0-based
                    depth = pileup_column.n
                    break

            coverage.append(CoverageAtPosition(
                chrom=chrom,
                pos=pos,
                depth=depth,
                is_variant=False,  # Will be updated when matching with VCF
            ))

    return coverage
