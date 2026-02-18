"""
Confidence score calculator for PGx Pipeline.

Computes confidence metrics based on:
- Read depth at variant positions
- Variant call quality
- Completeness of allele-defining variant coverage

Formula: confidence = w1 * depth_score + w2 * quality_score + w3 * completeness_score
"""

from dataclasses import dataclass

from elpgx.models import ConfidenceFlag, ConfidenceMetrics, Variant, VariantDepthStats


@dataclass
class ConfidenceConfig:
    """Configuration for confidence calculation."""

    # Weights for composite score
    depth_weight: float = 0.40
    quality_weight: float = 0.35
    completeness_weight: float = 0.25

    # Thresholds
    target_depth: int = 30
    max_quality: float = 100.0
    min_depth_threshold: int = 10
    min_quality_threshold: float = 20.0


def calculate_depth_score(
    variants: list[Variant],
    target_depth: int = 30,
) -> tuple[float, float, int]:
    """
    Calculate depth score component.

    Args:
        variants: List of variants with depth information.
        target_depth: Target depth for full score.

    Returns:
        Tuple of (depth_score, mean_depth, min_depth).
    """
    if not variants:
        return 0.0, 0.0, 0

    depths = [v.depth for v in variants]
    mean_depth = sum(depths) / len(depths)
    min_depth = min(depths)

    # Score: ratio of mean depth to target, capped at 1.0
    depth_score = min(1.0, mean_depth / target_depth)

    return depth_score, mean_depth, min_depth


def calculate_quality_score(
    variants: list[Variant],
    max_quality: float = 100.0,
) -> float:
    """
    Calculate quality score component.

    Average of individual variant quality scores, normalized.

    Args:
        variants: List of variants with quality information.
        max_quality: Maximum quality for full score.

    Returns:
        Quality score (0.0-1.0).
    """
    if not variants:
        return 0.0

    # Normalize each variant's quality and average
    normalized_scores = []
    for v in variants:
        score = min(1.0, v.quality / max_quality)
        normalized_scores.append(score)

    return sum(normalized_scores) / len(normalized_scores)


def calculate_completeness_score(
    genotyped_positions: int,
    total_defining_positions: int,
) -> float:
    """
    Calculate completeness score component.

    Ratio of genotyped defining positions to total defining positions.

    Args:
        genotyped_positions: Number of positions successfully genotyped.
        total_defining_positions: Total number of defining positions for allele.

    Returns:
        Completeness score (0.0-1.0).
    """
    if total_defining_positions == 0:
        return 1.0  # No defining variants = reference allele, fully complete

    return genotyped_positions / total_defining_positions


def determine_flags(
    variants: list[Variant],
    completeness_score: float,
    call_quality: str,
    config: ConfidenceConfig,
) -> list[ConfidenceFlag]:
    """
    Determine quality flags based on metrics.

    Args:
        variants: List of variants.
        completeness_score: Calculated completeness score.
        call_quality: Call quality from diplotype caller ("definitive", "ambiguous", "no_call").
        config: Confidence configuration.

    Returns:
        List of applicable ConfidenceFlag values.
    """
    flags: list[ConfidenceFlag] = []

    # Check depth
    for v in variants:
        if v.depth < config.min_depth_threshold:
            flags.append(ConfidenceFlag.LOW_DEPTH)
            break

    # Check quality
    for v in variants:
        if v.quality < config.min_quality_threshold:
            flags.append(ConfidenceFlag.LOW_QUALITY)
            break

    # Check completeness
    if completeness_score < 1.0:
        flags.append(ConfidenceFlag.INCOMPLETE_ALLELE)

    # Check call quality
    if call_quality == "ambiguous":
        flags.append(ConfidenceFlag.AMBIGUOUS_CALL)
    elif call_quality == "no_call":
        flags.append(ConfidenceFlag.NO_CALL)

    return flags


def calculate_confidence(
    variants: list[Variant],
    total_defining_positions: int = 0,
    call_quality: str = "definitive",
    config: ConfidenceConfig | None = None,
) -> ConfidenceMetrics:
    """
    Calculate confidence metrics for a gene call.

    Composite score formula:
    confidence = w1 * depth_score + w2 * quality_score + w3 * completeness_score

    Args:
        variants: List of variants that contributed to the call.
        total_defining_positions: Total defining positions for the called alleles.
        call_quality: Quality indicator from diplotype caller.
        config: Optional confidence configuration.

    Returns:
        ConfidenceMetrics with score, flags, and depth stats.
    """
    if config is None:
        config = ConfidenceConfig()

    # Calculate component scores
    depth_score, mean_depth, min_depth = calculate_depth_score(
        variants, config.target_depth
    )
    quality_score = calculate_quality_score(variants, config.max_quality)

    # Completeness: genotyped = number of variants found
    genotyped = len(variants)
    if total_defining_positions == 0:
        # Reference allele case: use number of variants as proxy
        total_defining_positions = max(genotyped, 1)
    completeness_score = calculate_completeness_score(genotyped, total_defining_positions)

    # Calculate composite score
    composite_score = (
        config.depth_weight * depth_score
        + config.quality_weight * quality_score
        + config.completeness_weight * completeness_score
    )

    # Round to 2 decimal places
    composite_score = round(composite_score, 2)

    # Determine flags
    flags = determine_flags(variants, completeness_score, call_quality, config)

    # Build depth stats
    depth_stats = None
    if variants:
        depth_stats = VariantDepthStats(
            mean_depth=round(mean_depth, 1),
            min_depth=min_depth,
        )

    return ConfidenceMetrics(
        score=composite_score,
        flags=flags,
        variant_depth_stats=depth_stats,
    )


def calculate_confidence_for_reference_call(
    coverage_at_positions: dict[int, int],
    config: ConfidenceConfig | None = None,
) -> ConfidenceMetrics:
    """
    Calculate confidence for a reference allele call (no variants found).

    Uses coverage data at defining positions to assess confidence.

    Args:
        coverage_at_positions: Dict mapping position to read depth.
        config: Optional confidence configuration.

    Returns:
        ConfidenceMetrics.
    """
    if config is None:
        config = ConfidenceConfig()

    if not coverage_at_positions:
        # No coverage data available
        return ConfidenceMetrics(
            score=0.5,  # Uncertain
            flags=[ConfidenceFlag.INCOMPLETE_ALLELE],
            variant_depth_stats=None,
        )

    depths = list(coverage_at_positions.values())
    mean_depth = sum(depths) / len(depths)
    min_depth = min(depths)

    # Depth score
    depth_score = min(1.0, mean_depth / config.target_depth)

    # Quality score: assume good quality for reference positions
    quality_score = 1.0

    # Completeness: all positions checked
    completeness_score = 1.0

    # Composite score
    composite_score = (
        config.depth_weight * depth_score
        + config.quality_weight * quality_score
        + config.completeness_weight * completeness_score
    )
    composite_score = round(composite_score, 2)

    # Flags
    flags: list[ConfidenceFlag] = []
    if min_depth < config.min_depth_threshold:
        flags.append(ConfidenceFlag.LOW_DEPTH)

    return ConfidenceMetrics(
        score=composite_score,
        flags=flags,
        variant_depth_stats=VariantDepthStats(
            mean_depth=round(mean_depth, 1),
            min_depth=min_depth,
        ),
    )
