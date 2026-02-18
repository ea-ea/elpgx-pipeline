"""
Allele caller module for PGx Pipeline.

Matches observed variants to known star alleles and determines diplotypes.
This is the core engine for pharmacogenetic allele assignment.
"""

import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from elpgx.models import (
    AlleleDefinition,
    AlleleVariant,
    BiomarkerType,
    ConfidenceFlag,
    ConfidenceMetrics,
    GeneAlleleDefinitions,
    GeneTier,
    Variant,
    VariantDepthStats,
)
from elpgx.vcf_parser import ParsedVariant, has_variant


@dataclass
class AlleleMatch:
    """Result of matching variants to a single allele."""

    allele: str
    function: str
    activity_value: float | None
    matched_variants: list[ParsedVariant]
    total_defining_variants: int
    is_complete_match: bool  # All defining variants found


@dataclass
class DiplotypeCall:
    """Result of diplotype calling for a gene."""

    gene: str
    allele1: str
    allele2: str
    diplotype: str  # "{allele1}/{allele2}"
    variants: list[Variant]
    allele1_function: str
    allele2_function: str
    allele1_activity: float | None
    allele2_activity: float | None
    confidence: ConfidenceMetrics = field(default_factory=lambda: ConfidenceMetrics(score=0.95, flags=[]))
    is_phased: bool = False
    call_quality: str = "definitive"  # "definitive", "ambiguous", "no_call"


def load_gene_definitions(allele_defs_dir: Path) -> dict[str, GeneAlleleDefinitions]:
    """
    Load all gene allele definitions from JSON files.

    Args:
        allele_defs_dir: Path to directory containing allele definition JSONs.

    Returns:
        Dict mapping gene name to GeneAlleleDefinitions.
    """
    definitions: dict[str, GeneAlleleDefinitions] = {}

    for json_file in allele_defs_dir.glob("*.json"):
        with open(json_file) as f:
            data = json.load(f)

        gene = data.get("gene")
        if not gene:
            continue

        # Parse alleles
        alleles = []
        for allele_data in data.get("alleles", []):
            defining_variants = []
            for var_data in allele_data.get("defining_variants", []):
                defining_variants.append(
                    AlleleVariant(
                        rsid=var_data.get("rsid"),
                        chrom=var_data["chrom"],
                        pos=var_data["pos"],
                        ref=var_data["ref"],
                        alt=var_data["alt"],
                        functional_impact=var_data.get("functional_impact"),
                    )
                )

            # Parse activity value, handling string cases like "≥3.0"
            activity_raw = allele_data.get("activity_value")
            activity_value = None
            if activity_raw is not None:
                if isinstance(activity_raw, (int, float)):
                    activity_value = float(activity_raw)
                elif isinstance(activity_raw, str):
                    # Extract numeric portion from strings like "≥3.0"
                    match = re.search(r'[\d.]+', activity_raw)
                    if match:
                        try:
                            activity_value = float(match.group())
                        except ValueError:
                            pass

            alleles.append(
                AlleleDefinition(
                    allele=allele_data["allele"],
                    function=allele_data.get("function") or "Unknown",
                    activity_value=activity_value,
                    defining_variants=defining_variants,
                    evidence_strength=allele_data.get("evidence_strength"),
                    citations=allele_data.get("citations"),
                )
            )

        definitions[gene] = GeneAlleleDefinitions(
            gene=gene,
            chromosome=data.get("chromosome", ""),
            gene_tier=GeneTier(data.get("gene_tier", "group1")),
            biomarker_type=BiomarkerType(data.get("biomarker_type", "metabolic_enzyme")),
            reference_allele=data.get("reference_allele", "*1"),
            alleles=alleles,
        )

    return definitions


def match_allele(
    allele_def: AlleleDefinition,
    variants: list[ParsedVariant],
) -> AlleleMatch:
    """
    Check if observed variants match a star allele definition.

    Args:
        allele_def: Allele definition with defining variants.
        variants: List of observed variants from VCF.

    Returns:
        AlleleMatch with matched variants and completeness status.
    """
    matched_variants: list[ParsedVariant] = []
    total_defining = len(allele_def.defining_variants)

    for def_var in allele_def.defining_variants:
        found = has_variant(
            variants,
            chrom=def_var.chrom,
            pos=def_var.pos,
            ref=def_var.ref,
            alt=def_var.alt,
        )
        if found:
            matched_variants.append(found)

    is_complete = len(matched_variants) == total_defining and total_defining > 0

    return AlleleMatch(
        allele=allele_def.allele,
        function=allele_def.function,
        activity_value=allele_def.activity_value,
        matched_variants=matched_variants,
        total_defining_variants=total_defining,
        is_complete_match=is_complete,
    )


def find_matching_alleles(
    gene_def: GeneAlleleDefinitions,
    variants: list[ParsedVariant],
) -> list[AlleleMatch]:
    """
    Find all alleles that match the observed variants for a gene.

    Args:
        gene_def: Gene allele definitions.
        variants: List of observed variants.

    Returns:
        List of AlleleMatch for alleles with at least one matched variant.
    """
    matches = []

    for allele_def in gene_def.alleles:
        # Skip reference allele (has no defining variants)
        if not allele_def.defining_variants:
            continue

        # Skip tandem/hybrid alleles (e.g., *13+*2, *36+*10) — these require
        # structural variant detection (long-read or specialized caller) and
        # cannot be reliably called from short-read SNV data alone
        if "+" in allele_def.allele:
            continue

        match = match_allele(allele_def, variants)
        if match.matched_variants:
            matches.append(match)

    return matches


def count_allele_copies(
    match: AlleleMatch,
    variants: list[ParsedVariant],
) -> int:
    """
    Determine how many copies of an allele are present (0, 1, or 2).

    For single-variant alleles, uses genotype (0/1 = 1 copy, 1/1 = 2 copies).
    For multi-variant alleles, requires phasing (simplified to 1 copy if any match).

    Args:
        match: AlleleMatch result.
        variants: List of observed variants.

    Returns:
        Number of copies (0, 1, or 2).
    """
    if not match.matched_variants:
        return 0

    # For single-variant alleles, use genotype directly
    if match.total_defining_variants == 1 and len(match.matched_variants) == 1:
        gt = match.matched_variants[0].genotype
        if gt in ("1/1", "1|1"):
            return 2
        elif "1" in gt:
            return 1
        return 0

    # For multi-variant alleles, simplified: if all defining variants present, 1 copy
    # (Full phasing would require haplotype reconstruction)
    if match.is_complete_match:
        return 1

    return 0


def determine_diplotype(
    gene_def: GeneAlleleDefinitions,
    variants: list[ParsedVariant],
) -> DiplotypeCall:
    """
    Determine diplotype for a gene based on observed variants.

    Algorithm:
    1. Find all matching non-reference alleles
    2. Count copies of each matched allele
    3. Assign diplotype based on copy numbers
    4. Fill remaining with reference allele

    Args:
        gene_def: Gene allele definitions.
        variants: List of observed variants.

    Returns:
        DiplotypeCall with diplotype assignment.
    """
    ref_allele = gene_def.reference_allele
    ref_function = "Normal Function"
    ref_activity: float | None = None

    # Find reference allele definition to get its function
    for allele_def in gene_def.alleles:
        if allele_def.allele == ref_allele:
            ref_function = allele_def.function
            ref_activity = allele_def.activity_value
            break

    # Find matching alleles
    matches = find_matching_alleles(gene_def, variants)

    if not matches:
        # No variants found -> homozygous reference
        # Flag that this is inferred from absence of defining variants
        return DiplotypeCall(
            gene=gene_def.gene,
            allele1=ref_allele,
            allele2=ref_allele,
            diplotype=f"{ref_allele}/{ref_allele}",
            variants=[],
            allele1_function=ref_function,
            allele2_function=ref_function,
            allele1_activity=ref_activity,
            allele2_activity=ref_activity,
            confidence=ConfidenceMetrics(
                score=0.90,
                flags=[ConfidenceFlag.NO_DEFINING_VARIANTS],
            ),
            call_quality="definitive",
        )

    # Count copies of each matched allele
    allele_copies: list[tuple[AlleleMatch, int]] = []
    for match in matches:
        copies = count_allele_copies(match, variants)
        if copies > 0:
            allele_copies.append((match, copies))

    # Sort by specificity (more defining variants = more specific)
    allele_copies.sort(key=lambda x: (-x[0].total_defining_variants, x[0].allele))

    # Assign diplotype
    assigned_alleles: list[tuple[str, str, float | None]] = []  # (allele, function, activity)
    total_assigned = 0

    for match, copies in allele_copies:
        for _ in range(copies):
            if total_assigned < 2:
                assigned_alleles.append(
                    (match.allele, match.function, match.activity_value)
                )
                total_assigned += 1

    # Fill remaining with reference allele
    while total_assigned < 2:
        assigned_alleles.append((ref_allele, ref_function, ref_activity))
        total_assigned += 1

    # Sort alleles for consistent representation
    assigned_alleles.sort(key=lambda x: x[0])

    allele1, func1, act1 = assigned_alleles[0]
    allele2, func2, act2 = assigned_alleles[1]

    # Collect all variants that contributed to the call, deduplicated by position
    all_variants: list[Variant] = []
    seen_variants: set[tuple[str, int, str, str, str]] = set()
    for match, _ in allele_copies:
        for pv in match.matched_variants:
            var_key = (pv.chrom, pv.pos, pv.ref, pv.alt, pv.genotype)
            if var_key in seen_variants:
                continue
            seen_variants.add(var_key)

            # Find the allele definition for functional impact
            functional_impact = None
            for allele_def in gene_def.alleles:
                if allele_def.allele == match.allele:
                    for def_var in allele_def.defining_variants:
                        if def_var.pos == pv.pos:
                            functional_impact = def_var.functional_impact
                            break

            all_variants.append(
                Variant(
                    rsid=pv.rsid,
                    chrom=pv.chrom,
                    pos=pv.pos,
                    ref=pv.ref,
                    alt=pv.alt,
                    genotype=pv.genotype,
                    quality=pv.quality,
                    depth=pv.depth,
                    allele_defining=True,
                    star_allele=match.allele,
                    functional_impact=functional_impact,
                )
            )

    # Determine call quality and confidence with informative flags
    call_quality = "definitive"
    confidence_score = 0.95
    confidence_flags: list[ConfidenceFlag] = []
    depth_stats: VariantDepthStats | None = None

    if len(allele_copies) > 2:
        call_quality = "ambiguous"
        confidence_score = 0.70
        confidence_flags.append(ConfidenceFlag.AMBIGUOUS_CALL)

    if len(all_variants) == 0:
        confidence_score = min(confidence_score, 0.85)
        confidence_flags.append(ConfidenceFlag.NO_DEFINING_VARIANTS)

    # Calculate depth stats and depth-based flags
    if all_variants:
        depths = [v.depth for v in all_variants if v.depth is not None]
        qualities = [v.quality for v in all_variants if v.quality is not None]

        if depths:
            mean_depth = sum(depths) / len(depths)
            min_depth = min(depths)
            depth_stats = VariantDepthStats(
                mean_depth=round(mean_depth, 1),
                min_depth=min_depth,
            )

            if mean_depth < 10:
                confidence_score = min(confidence_score, 0.60)
                confidence_flags.append(ConfidenceFlag.LOW_DEPTH)
            elif mean_depth < 20:
                confidence_score = min(confidence_score, 0.80)
                confidence_flags.append(ConfidenceFlag.LOW_DEPTH)

        if qualities:
            min_quality = min(qualities)
            if min_quality < 20:
                confidence_score = min(confidence_score, 0.75)
                confidence_flags.append(ConfidenceFlag.LOW_QUALITY)

    return DiplotypeCall(
        gene=gene_def.gene,
        allele1=allele1,
        allele2=allele2,
        diplotype=f"{allele1}/{allele2}",
        variants=all_variants,
        allele1_function=func1,
        allele2_function=func2,
        allele1_activity=act1,
        allele2_activity=act2,
        confidence=ConfidenceMetrics(
            score=confidence_score,
            flags=confidence_flags,
            variant_depth_stats=depth_stats,
        ),
        call_quality=call_quality,
    )


def call_alleles_for_gene(
    gene: str,
    variants: list[ParsedVariant],
    allele_defs_dir: Path,
) -> DiplotypeCall | None:
    """
    Call diplotype for a single gene.

    Args:
        gene: Gene name (e.g., "CYP2C19").
        variants: List of observed variants.
        allele_defs_dir: Path to allele definitions directory.

    Returns:
        DiplotypeCall or None if gene definition not found.
    """
    json_file = allele_defs_dir / f"{gene}.json"
    if not json_file.exists():
        return None

    with open(json_file) as f:
        data = json.load(f)

    # Parse to GeneAlleleDefinitions
    alleles = []
    for allele_data in data.get("alleles", []):
        defining_variants = []
        for var_data in allele_data.get("defining_variants", []):
            defining_variants.append(
                AlleleVariant(
                    rsid=var_data.get("rsid"),
                    chrom=var_data["chrom"],
                    pos=var_data["pos"],
                    ref=var_data["ref"],
                    alt=var_data["alt"],
                    functional_impact=var_data.get("functional_impact"),
                )
            )

        # Parse activity value, handling string cases like "≥3.0"
        activity_raw = allele_data.get("activity_value")
        activity_value = None
        if activity_raw is not None:
            if isinstance(activity_raw, (int, float)):
                activity_value = float(activity_raw)
            elif isinstance(activity_raw, str):
                match = re.search(r'[\d.]+', activity_raw)
                if match:
                    try:
                        activity_value = float(match.group())
                    except ValueError:
                        pass

        alleles.append(
            AlleleDefinition(
                allele=allele_data["allele"],
                function=allele_data.get("function") or "Unknown",
                activity_value=activity_value,
                defining_variants=defining_variants,
                evidence_strength=allele_data.get("evidence_strength"),
                citations=allele_data.get("citations"),
            )
        )

    gene_def = GeneAlleleDefinitions(
        gene=data["gene"],
        chromosome=data.get("chromosome", ""),
        gene_tier=GeneTier(data.get("gene_tier", "group1")),
        biomarker_type=BiomarkerType(data.get("biomarker_type", "metabolic_enzyme")),
        reference_allele=data.get("reference_allele", "*1"),
        alleles=alleles,
    )

    return determine_diplotype(gene_def, variants)


def call_all_genes(
    variants: list[ParsedVariant],
    allele_defs_dir: Path,
    genes: list[str] | None = None,
) -> list[DiplotypeCall]:
    """
    Call diplotypes for all genes (or specified subset).

    Args:
        variants: List of observed variants from VCF.
        allele_defs_dir: Path to allele definitions directory.
        genes: Optional list of genes to call. If None, calls all.

    Returns:
        List of DiplotypeCall for each gene.
    """
    definitions = load_gene_definitions(allele_defs_dir)

    if genes is None:
        genes = list(definitions.keys())

    calls = []
    for gene in genes:
        if gene not in definitions:
            continue
        call = determine_diplotype(definitions[gene], variants)
        calls.append(call)

    return calls
