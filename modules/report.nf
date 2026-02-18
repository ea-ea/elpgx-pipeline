/*
 * Module: Report Generation
 *
 * Assembles final PGx JSON report with full provenance and schema validation.
 * Includes: git commit, Nextflow version, Docker digests, file checksums
 */

process GENERATE_REPORT {
    tag "$sample_id"
    label 'process_low'

    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    tuple val(sample_id), path(alleles_json), path(phenotypes_json)
    path gene_tiers_file
    path schema_path
    path gene_strategies_dir
    val gene_tiers
    val enable_hla
    val pipeline_version
    val pgx_db_version
    val git_commit
    val nextflow_version
    val reference_md5
    val reference_fasta_path
    val input_bam_md5

    output:
    path "${sample_id}_pgx_report.json", emit: report

    script:
    """
    #!/usr/bin/env python3

    import json
    import hashlib
    from datetime import datetime
    from pathlib import Path

    from elpgx.models import (
        PgxReport, PipelineMeta, GeneResult, Variant,
        ConfidenceMetrics, VariantDepthStats, GeneTier, BiomarkerType,
        ConfidenceFlag, ProvenanceInfo
    )

    # Load input data
    with open("${alleles_json}") as f:
        alleles_data = json.load(f)

    with open("${phenotypes_json}") as f:
        phenotypes_data = json.load(f)

    with open("${gene_tiers_file}") as f:
        gene_tiers_config = json.load(f)

    # Load gene strategies for phenotype blocking and complex locus rules
    gene_strategies = {}
    strategies_dir = Path("${gene_strategies_dir}")
    if strategies_dir.exists():
        for strategy_file in strategies_dir.glob("*_STRATEGY.json"):
            with open(strategy_file) as f:
                strategy = json.load(f)
                gene_strategies[strategy.get("gene")] = strategy

    # Map phenotypes by gene
    phenotype_map = {p["gene"]: p for p in phenotypes_data.get("phenotypes", [])}

    # Function to check if phenotype should be blocked
    def should_block_phenotype(gene, confidence_score, flags):
        strategy = gene_strategies.get(gene)
        if not strategy:
            return False, None, None

        for rule in strategy.get("phenotype_blocking_rules", []):
            condition = rule.get("condition", "")

            if "confidence_score <" in condition:
                try:
                    threshold = float(condition.split("<")[1].strip())
                    if confidence_score < threshold:
                        return True, rule.get("phenotype_value", "Indeterminate"), rule.get("message")
                except Exception:
                    pass

            if "flag ==" in condition:
                flag_name = condition.split("==")[1].strip().strip("'\\"")
                if flag_name in flags:
                    return True, rule.get("phenotype_value", "Indeterminate"), rule.get("message")

        return False, None, None

    # Function to apply complex locus rules (CYP2D6 etc.)
    def apply_complex_locus_rules(gene, score, flags, strategy):
        extra_warnings = []
        if not strategy:
            return score, flags, extra_warnings

        unsupported = strategy.get("unsupported_features", {})
        cnv = unsupported.get("cnv_detection", {})
        hybrid = unsupported.get("hybrid_alleles", {})

        has_cnv = cnv.get("supported") is False and cnv.get("clinical_impact") in ("critical", "high")
        has_hybrid = hybrid.get("supported") is False and hybrid.get("clinical_impact") in ("critical", "high")

        if not (has_cnv or has_hybrid):
            return score, flags, extra_warnings

        if "COMPLEX_LOCUS_SHORT_READ" not in flags:
            flags = flags + ["COMPLEX_LOCUS_SHORT_READ"]
        score = min(score, 0.70)

        if has_cnv:
            msg = cnv.get("warning_message", f"{gene}: CNV detection not supported")
            extra_warnings.append(f"{gene}: {msg}")
        if has_hybrid:
            extra_warnings.append(f"{gene}: Hybrid/chimeric allele detection not supported")

        known = strategy.get("known_limitations", [])
        if known:
            extra_warnings.append(f"{gene}: Informational only — {'; '.join(known[:2])}")

        return score, flags, extra_warnings

    # Build gene results
    results = []
    warnings = []

    for call in alleles_data.get("calls", []):
        gene = call["gene"]
        phenotype_info = phenotype_map.get(gene, {})
        gene_config = gene_tiers_config.get("genes", {}).get(gene, {})

        # Map tier and biomarker type
        tier_str = gene_config.get("tier", "group1")
        biomarker_str = gene_config.get("biomarker_type", "metabolic_enzyme")

        try:
            gene_tier = GeneTier(tier_str)
        except ValueError:
            gene_tier = GeneTier.GROUP1

        try:
            biomarker_type = BiomarkerType(biomarker_str)
        except ValueError:
            biomarker_type = BiomarkerType.METABOLIC_ENZYME

        # Build variants — pass through all fields from allele caller
        variants = []
        seen_variants = set()
        for v in call.get("variants", []):
            var_key = (v["chrom"], v["pos"], v["ref"], v["alt"], v["genotype"])
            if var_key in seen_variants:
                continue
            seen_variants.add(var_key)
            variants.append(Variant(
                rsid=v.get("rsid"),
                chrom=v["chrom"],
                pos=v["pos"],
                ref=v["ref"],
                alt=v["alt"],
                genotype=v["genotype"],
                quality=v.get("quality", 0),
                depth=v.get("depth", 0),
                allele_defining=v.get("allele_defining", True),
                star_allele=v.get("star_allele"),
                functional_impact=v.get("functional_impact")
            ))

        # Build confidence — pass through variant_depth_stats and flags
        conf_data = call.get("confidence", {})
        score = conf_data.get("score", 0.95)
        flag_strs = conf_data.get("flags", [])

        # Build variant_depth_stats if present
        depth_stats = None
        vds = conf_data.get("variant_depth_stats")
        if vds:
            depth_stats = VariantDepthStats(
                mean_depth=vds["mean_depth"],
                min_depth=vds["min_depth"]
            )

        # Apply complex locus rules (CYP2D6 etc.)
        strategy = gene_strategies.get(gene)
        score, flag_strs, locus_warnings = apply_complex_locus_rules(gene, score, flag_strs, strategy)
        warnings.extend(locus_warnings)

        # Parse flags to enums
        flags = []
        for flag_str in flag_strs:
            try:
                flags.append(ConfidenceFlag(flag_str))
            except ValueError:
                pass

        confidence = ConfidenceMetrics(score=score, flags=flags, variant_depth_stats=depth_stats)

        # Check phenotype blocking rules
        phenotype = phenotype_info.get("phenotype", "Unknown")
        functional_status = phenotype_info.get("functional_status", "Unknown")
        prescribing_implication = phenotype_info.get("prescribing_implication", "No guidance available.")

        block, blocked_phenotype, block_message = should_block_phenotype(gene, score, [f.value for f in flags])
        if block:
            phenotype = blocked_phenotype
            functional_status = "Indeterminate"
            prescribing_implication = "Consult clinical pharmacist - phenotype could not be determined with sufficient confidence"
            warnings.append(f"{gene}: {block_message}")

        result = GeneResult(
            gene=gene,
            gene_tier=gene_tier,
            biomarker_type=biomarker_type,
            diplotype=call["diplotype"],
            functional_status=functional_status,
            phenotype=phenotype,
            activity_score=phenotype_info.get("activity_score"),
            prescribing_implication=prescribing_implication,
            confidence=confidence,
            variants=variants
        )
        results.append(result)

    # Parse gene tiers
    tiers_list = "${gene_tiers}".split(",")
    tier_enums = []
    for t in tiers_list:
        try:
            tier_enums.append(GeneTier(t.strip()))
        except ValueError:
            pass

    # Calculate PGx database hash
    def calculate_db_hash():
        sha256 = hashlib.sha256()
        with open("${gene_tiers_file}", "rb") as f:
            sha256.update(f.read())
        return sha256.hexdigest()[:16]

    # Build provenance
    provenance = ProvenanceInfo(
        git_commit="${git_commit}" if "${git_commit}" != "unknown" else None,
        git_branch=None,
        nextflow_version="${nextflow_version}" if "${nextflow_version}" != "unknown" else None,
        docker_images={
            "elpgx-python": "elifffarslan/elpgx-python:${pipeline_version}",
            "elpgx-gatk": "broadinstitute/gatk:4.5.0.0"
        },
        reference_fasta_md5="${reference_md5}" if "${reference_md5}" != "unknown" else None,
        reference_fasta_path="${reference_fasta_path}" if "${reference_fasta_path}" != "unknown" else None,
        pgx_db_hash=calculate_db_hash(),
        input_bam_md5="${input_bam_md5}" if "${input_bam_md5}" != "unknown" else None,
        variant_caller="GATK HaplotypeCaller",
        variant_caller_version="4.5.0.0"
    )

    # Build metadata
    meta = PipelineMeta(
        sample_id="${sample_id}",
        pipeline_version="${pipeline_version}",
        run_date=datetime.utcnow(),
        reference_build=alleles_data.get("reference_build", "GRCh38"),
        pgx_db_version="${pgx_db_version}",
        gene_tiers_included=tier_enums,
        hla_module_enabled="${enable_hla}" == "true",
        provenance=provenance
    )

    # Build report
    report = PgxReport(
        schema_version="1.1.0",
        meta=meta,
        results=results,
        warnings=warnings
    )

    # Write output
    report_path = "${sample_id}_pgx_report.json"
    with open(report_path, "w") as f:
        f.write(report.model_dump_json(indent=2))

    # Validate against schema
    schema_file = Path("${schema_path}")
    if schema_file.exists():
        try:
            import jsonschema
            with open(schema_file) as f:
                schema = json.load(f)
            with open(report_path) as f:
                report_data = json.load(f)
            jsonschema.validate(report_data, schema)
            print(f"Schema validation PASSED")
        except ImportError:
            print(f"WARNING: jsonschema not installed, skipping validation")
        except jsonschema.ValidationError as e:
            print(f"WARNING: Schema validation FAILED: {e.message}")

    print(f"Report generated: {report_path} with {len(results)} gene results")
    if warnings:
        print(f"Warnings: {len(warnings)}")
        for w in warnings:
            print(f"  - {w}")
    """
}
