"""
Report generation module for PGx Pipeline.

Builds the final JSON report by aggregating all results:
- Gene results with diplotypes and phenotypes
- Confidence metrics
- Pipeline metadata with full provenance
- Warnings
"""

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from elpgx.allele_caller import DiplotypeCall
from elpgx.confidence import ConfidenceConfig, calculate_confidence
from elpgx.models import (
    BiomarkerType,
    ConfidenceFlag,
    ConfidenceMetrics,
    GeneResult,
    GeneTier,
    PgxReport,
    PipelineMeta,
    ProvenanceInfo,
)
from elpgx.phenotype import PhenotypeResult
from elpgx.provenance import collect_provenance, validate_provenance


def load_gene_tier_info(gene_tiers_file: Path) -> dict[str, dict]:
    """
    Load gene tier and biomarker type information.

    Args:
        gene_tiers_file: Path to gene_tiers.json.

    Returns:
        Dict mapping gene name to tier info.
    """
    if not gene_tiers_file.exists():
        return {}

    with open(gene_tiers_file) as f:
        data = json.load(f)

    return data.get("genes", {})


def build_gene_result(
    diplotype_call: DiplotypeCall,
    phenotype_result: PhenotypeResult,
    gene_tier_info: dict,
    confidence_config: ConfidenceConfig | None = None,
) -> GeneResult:
    """
    Build a GeneResult from diplotype call and phenotype result.

    Args:
        diplotype_call: DiplotypeCall from allele caller.
        phenotype_result: PhenotypeResult from phenotype translator.
        gene_tier_info: Tier info for this gene.
        confidence_config: Optional confidence calculation config.

    Returns:
        GeneResult model.
    """
    gene = diplotype_call.gene

    # Get tier and biomarker type from config, with defaults
    tier_str = gene_tier_info.get("tier", "group1")
    biomarker_str = gene_tier_info.get("biomarker_type", "metabolic_enzyme")

    try:
        gene_tier = GeneTier(tier_str)
    except ValueError:
        gene_tier = GeneTier.GROUP1

    try:
        biomarker_type = BiomarkerType(biomarker_str)
    except ValueError:
        biomarker_type = BiomarkerType.METABOLIC_ENZYME

    # Use confidence from diplotype call (already calculated during allele calling)
    # This handles reference calls properly where no variants are found
    confidence = diplotype_call.confidence

    return GeneResult(
        gene=gene,
        gene_tier=gene_tier,
        biomarker_type=biomarker_type,
        diplotype=diplotype_call.diplotype,
        functional_status=phenotype_result.functional_status,
        phenotype=phenotype_result.phenotype,
        activity_score=phenotype_result.activity_score,
        prescribing_implication=phenotype_result.prescribing_implication,
        confidence=confidence,
        variants=diplotype_call.variants,
    )


def load_gene_strategy(data_dir: Path, gene: str) -> Optional[dict]:
    """
    Load gene-specific strategy configuration.

    Args:
        data_dir: Path to data directory.
        gene: Gene name.

    Returns:
        Strategy dict or None if not found.
    """
    strategy_file = data_dir / "gene_strategies" / f"{gene}_STRATEGY.json"
    if strategy_file.exists():
        with open(strategy_file) as f:
            return json.load(f)
    return None


def apply_phenotype_blocking(
    gene_result: GeneResult,
    strategy: Optional[dict],
) -> tuple[GeneResult, list[str]]:
    """
    Apply phenotype blocking rules based on gene strategy.

    For low-confidence calls, phenotype may be blocked and set to "Indeterminate".

    Args:
        gene_result: The gene result to potentially modify.
        strategy: Gene-specific strategy config.

    Returns:
        Tuple of (possibly modified GeneResult, list of warnings).
    """
    warnings = []

    if strategy is None:
        return gene_result, warnings

    confidence_thresholds = strategy.get("confidence_thresholds", {})
    blocking_rules = strategy.get("phenotype_blocking_rules", [])

    # Check blocking rules
    for rule in blocking_rules:
        condition = rule.get("condition", "")
        action = rule.get("action", "")

        should_block = False

        # Parse simple conditions
        if "confidence_score <" in condition:
            try:
                threshold = float(condition.split("<")[1].strip())
                if gene_result.confidence.score < threshold:
                    should_block = True
            except (ValueError, IndexError):
                pass

        if "flag ==" in condition:
            flag_name = condition.split("==")[1].strip().strip("'\"")
            flag_values = [f.value for f in gene_result.confidence.flags]
            if flag_name in flag_values:
                should_block = True

        if should_block and action == "block_phenotype":
            # Create modified result with blocked phenotype
            blocked_phenotype = rule.get("phenotype_value", "Indeterminate")
            message = rule.get("message", "Phenotype blocked due to low confidence")

            warnings.append(f"{gene_result.gene}: {message}")

            # Create new GeneResult with blocked phenotype
            gene_result = GeneResult(
                gene=gene_result.gene,
                gene_tier=gene_result.gene_tier,
                biomarker_type=gene_result.biomarker_type,
                diplotype=gene_result.diplotype,
                functional_status="Indeterminate",
                phenotype=blocked_phenotype,
                activity_score=None,
                prescribing_implication="Consult clinical pharmacist - phenotype could not be determined with sufficient confidence",
                confidence=gene_result.confidence,
                variants=gene_result.variants,
            )
            break  # Stop after first blocking rule matches

    return gene_result, warnings


def apply_complex_locus_rules(
    gene_result: GeneResult,
    strategy: Optional[dict],
) -> tuple[GeneResult, list[str]]:
    """
    Apply complex locus rules for genes with unsupported SV features.

    For genes like CYP2D6 where CNV/hybrid detection is not supported,
    cap confidence and add warning flags so the result is clearly
    marked as limited by short-read sequencing.

    Args:
        gene_result: The gene result to potentially modify.
        strategy: Gene-specific strategy config.

    Returns:
        Tuple of (possibly modified GeneResult, list of warnings).
    """
    warnings: list[str] = []

    if strategy is None:
        return gene_result, warnings

    unsupported = strategy.get("unsupported_features", {})
    cnv = unsupported.get("cnv_detection", {})
    hybrid = unsupported.get("hybrid_alleles", {})

    has_cnv_limitation = cnv.get("supported") is False and cnv.get("clinical_impact") in ("critical", "high")
    has_hybrid_limitation = hybrid.get("supported") is False and hybrid.get("clinical_impact") in ("critical", "high")

    if not (has_cnv_limitation or has_hybrid_limitation):
        return gene_result, warnings

    # Cap confidence score — CNV status unknown means result is incomplete
    current_flags = list(gene_result.confidence.flags)
    if ConfidenceFlag.COMPLEX_LOCUS_SHORT_READ not in current_flags:
        current_flags.append(ConfidenceFlag.COMPLEX_LOCUS_SHORT_READ)

    capped_score = min(gene_result.confidence.score, 0.70)

    gene_result = GeneResult(
        gene=gene_result.gene,
        gene_tier=gene_result.gene_tier,
        biomarker_type=gene_result.biomarker_type,
        diplotype=gene_result.diplotype,
        functional_status=gene_result.functional_status,
        phenotype=gene_result.phenotype,
        activity_score=gene_result.activity_score,
        prescribing_implication=gene_result.prescribing_implication,
        confidence=ConfidenceMetrics(
            score=capped_score,
            flags=current_flags,
            variant_depth_stats=gene_result.confidence.variant_depth_stats,
        ),
        variants=gene_result.variants,
    )

    # Add pipeline-level warnings from the strategy
    gene = gene_result.gene
    if has_cnv_limitation:
        cnv_warning = cnv.get("warning_message", f"{gene}: CNV detection not supported — results assume diploid copy number")
        warnings.append(f"{gene}: {cnv_warning}")
    if has_hybrid_limitation:
        warnings.append(f"{gene}: Hybrid/chimeric allele detection not supported — requires long-read sequencing")

    known_limitations = strategy.get("known_limitations", [])
    if known_limitations:
        warnings.append(f"{gene}: Informational only — {'; '.join(known_limitations[:2])}")

    return gene_result, warnings


def build_report(
    sample_id: str,
    diplotype_calls: list[DiplotypeCall],
    phenotype_results: list[PhenotypeResult],
    gene_tiers_file: Path,
    gene_tiers_included: list[str],
    hla_enabled: bool = False,
    pgx_db_version: str = "CPIC v1.24",
    pipeline_version: str = "1.0.0",
    reference_build: str = "GRCh38",
    confidence_config: ConfidenceConfig | None = None,
    provenance: ProvenanceInfo | None = None,
    data_dir: Path | None = None,
    reference_path: Path | None = None,
    input_bam_path: Path | None = None,
) -> PgxReport:
    """
    Build the complete PGx report.

    Args:
        sample_id: Sample identifier.
        diplotype_calls: List of DiplotypeCall from allele caller.
        phenotype_results: List of PhenotypeResult from phenotype translator.
        gene_tiers_file: Path to gene_tiers.json.
        gene_tiers_included: List of tier names included in analysis.
        hla_enabled: Whether HLA module was enabled.
        pgx_db_version: Version of PGx database used.
        pipeline_version: Pipeline version string.
        confidence_config: Optional confidence calculation config.
        provenance: Optional pre-collected provenance info.
        data_dir: Path to data directory (for strategy loading and provenance).
        reference_path: Path to reference FASTA (for provenance).
        input_bam_path: Path to input BAM (for provenance).

    Returns:
        PgxReport model.
    """
    # Load gene tier info
    gene_tier_info = load_gene_tier_info(gene_tiers_file)

    # Map phenotype results by gene for easy lookup
    phenotype_by_gene = {r.gene: r for r in phenotype_results}

    # Build gene results
    results: list[GeneResult] = []
    warnings: list[str] = []

    for call in diplotype_calls:
        gene = call.gene
        phenotype = phenotype_by_gene.get(gene)

        if phenotype is None:
            warnings.append(f"No phenotype translation available for {gene}")
            continue

        tier_info = gene_tier_info.get(gene, {})
        gene_result = build_gene_result(
            diplotype_call=call,
            phenotype_result=phenotype,
            gene_tier_info=tier_info,
            confidence_config=confidence_config,
        )

        # Apply gene strategy rules (phenotype blocking + complex locus warnings)
        if data_dir:
            strategy = load_gene_strategy(data_dir, gene)
            gene_result, blocking_warnings = apply_phenotype_blocking(gene_result, strategy)
            warnings.extend(blocking_warnings)
            gene_result, locus_warnings = apply_complex_locus_rules(gene_result, strategy)
            warnings.extend(locus_warnings)

        results.append(gene_result)

    # Build metadata
    tier_enums = []
    for tier_name in gene_tiers_included:
        try:
            tier_enums.append(GeneTier(tier_name))
        except ValueError:
            pass

    # Collect provenance if not provided
    if provenance is None:
        docker_images = [
            f"elifffarslan/elpgx-python:{pipeline_version}",
            f"elifffarslan/elpgx-gatk:{pipeline_version}",
        ]
        if hla_enabled:
            docker_images.append(f"elifffarslan/elpgx-hla:{pipeline_version}")

        provenance = collect_provenance(
            reference_path=reference_path,
            input_bam_path=input_bam_path,
            data_dir=data_dir,
            docker_images=docker_images,
        )

    # Validate provenance and add warnings for missing info
    provenance_issues = validate_provenance(provenance)
    for issue in provenance_issues:
        warnings.append(f"Provenance: {issue}")

    meta = PipelineMeta(
        sample_id=sample_id,
        pipeline_version=pipeline_version,
        run_date=datetime.now(timezone.utc),
        reference_build=reference_build,
        pgx_db_version=pgx_db_version,
        gene_tiers_included=tier_enums,
        hla_module_enabled=hla_enabled,
        provenance=provenance,
    )

    return PgxReport(
        meta=meta,
        results=results,
        warnings=warnings,
    )


def write_report(report: PgxReport, output_path: Path) -> None:
    """
    Write report to JSON file.

    Args:
        report: PgxReport to write.
        output_path: Path to output file.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as f:
        f.write(report.to_json(indent=2))


def generate_report(
    sample_id: str,
    diplotype_calls: list[DiplotypeCall],
    phenotype_results: list[PhenotypeResult],
    output_path: Path,
    data_dir: Path,
    gene_tiers_included: list[str] | None = None,
    hla_enabled: bool = False,
    pgx_db_version: str = "CPIC v1.24",
    pipeline_version: str = "1.0.0",
    reference_build: str = "GRCh38",
    reference_path: Path | None = None,
    input_bam_path: Path | None = None,
) -> PgxReport:
    """
    High-level function to generate and write the PGx report.

    Args:
        sample_id: Sample identifier.
        diplotype_calls: List of DiplotypeCall from allele caller.
        phenotype_results: List of PhenotypeResult from phenotype translator.
        output_path: Path to write JSON report.
        data_dir: Path to data directory containing gene_tiers.json.
        gene_tiers_included: List of tier names (default: ["group1"]).
        hla_enabled: Whether HLA module was enabled.
        pgx_db_version: Version of PGx database.
        pipeline_version: Pipeline version.
        reference_build: Reference genome build.
        reference_path: Path to reference FASTA (for provenance).
        input_bam_path: Path to input BAM (for provenance).

    Returns:
        Generated PgxReport.
    """
    if gene_tiers_included is None:
        gene_tiers_included = ["group1"]

    gene_tiers_file = data_dir / "gene_tiers.json"

    report = build_report(
        sample_id=sample_id,
        diplotype_calls=diplotype_calls,
        phenotype_results=phenotype_results,
        gene_tiers_file=gene_tiers_file,
        gene_tiers_included=gene_tiers_included,
        hla_enabled=hla_enabled,
        pgx_db_version=pgx_db_version,
        pipeline_version=pipeline_version,
        reference_build=reference_build,
        data_dir=data_dir,
        reference_path=reference_path,
        input_bam_path=input_bam_path,
    )

    write_report(report, output_path)

    return report
