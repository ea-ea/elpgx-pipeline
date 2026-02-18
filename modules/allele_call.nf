/*
 * Module: Allele Calling
 *
 * Determines star-allele diplotypes from VCF variants.
 * Supports GRCh37 input with automatic coordinate translation.
 */

process ALLELE_CALL {
    tag "$sample_id"
    label 'process_medium'

    publishDir "${params.outdir}/alleles", mode: params.publish_dir_mode

    input:
    tuple val(sample_id), path(vcf), path(vcf_index)
    path allele_defs_dir
    val gene_tiers
    path validation_json
    path grch37_coord_map

    output:
    tuple val(sample_id), path("alleles.json"), emit: alleles

    script:
    """
    #!/usr/bin/env python3

    import json
    import sys
    from pathlib import Path

    from elpgx.vcf_parser import parse_vcf, load_target_positions, filter_by_positions
    from elpgx.allele_caller import load_gene_definitions, determine_diplotype

    # Detect build from validation result
    with open("${validation_json}") as f:
        validation = json.load(f)
    detected_build = validation.get("reference_build", "GRCh38")

    # Parse gene tiers
    tiers = "${gene_tiers}".split(",")

    # Parse VCF
    vcf_path = Path("${vcf}")
    all_variants = parse_vcf(vcf_path)

    # If GRCh37, translate coordinates to GRCh38 before allele matching
    if detected_build == "GRCh37":
        from elpgx.vcf_parser import translate_variants_grch37_to_grch38
        coord_map_path = Path("${grch37_coord_map}")
        with open(coord_map_path) as f:
            raw_map = json.load(f)
        # Build grch37->grch38 lookup
        coord_map = {}
        for gene_name, positions in raw_map.items():
            for grch38_pos_str, mapping in positions.items():
                grch37_chrom = mapping["grch37_chrom"]
                grch37_pos = mapping["grch37_pos"]
                grch38_pos = int(grch38_pos_str)
                grch38_chrom = "chr" + grch37_chrom if not grch37_chrom.startswith("chr") else grch37_chrom
                if grch37_chrom not in coord_map:
                    coord_map[grch37_chrom] = {}
                coord_map[grch37_chrom][grch37_pos] = (grch38_chrom, grch38_pos)
        all_variants = translate_variants_grch37_to_grch38(all_variants, coord_map)

    # Load gene definitions from directory (returns dict of gene_name -> GeneAlleleDefinitions)
    allele_defs_path = Path("${allele_defs_dir}")
    gene_defs = load_gene_definitions(allele_defs_path)

    # Load target positions and filter variants
    target_positions = load_target_positions(allele_defs_path)
    pgx_variants = filter_by_positions(all_variants, target_positions)

    # Call alleles for each gene in selected tiers
    calls = []
    for gene_name, gene_def in gene_defs.items():
        if gene_def.gene_tier.value not in tiers:
            continue

        # Get variants on this gene's chromosome
        gene_variants = [v for v in pgx_variants if v.chrom == gene_def.chromosome]
        diplotype_call = determine_diplotype(gene_def, gene_variants)

        call_entry = {
            "gene": diplotype_call.gene,
            "diplotype": diplotype_call.diplotype,
            "haplotypes": [diplotype_call.allele1, diplotype_call.allele2],
            "allele1_function": diplotype_call.allele1_function,
            "allele2_function": diplotype_call.allele2_function,
            "allele1_activity": diplotype_call.allele1_activity,
            "allele2_activity": diplotype_call.allele2_activity,
            "call_quality": diplotype_call.call_quality,
            "confidence": {
                "score": diplotype_call.confidence.score,
                "flags": [f.value for f in diplotype_call.confidence.flags],
                "variant_depth_stats": {
                    "mean_depth": diplotype_call.confidence.variant_depth_stats.mean_depth,
                    "min_depth": diplotype_call.confidence.variant_depth_stats.min_depth,
                } if diplotype_call.confidence.variant_depth_stats else None
            },
            "variants": [
                {
                    "rsid": v.rsid,
                    "chrom": v.chrom,
                    "pos": v.pos,
                    "ref": v.ref,
                    "alt": v.alt,
                    "genotype": v.genotype,
                    "depth": v.depth,
                    "quality": v.quality,
                    "allele_defining": v.allele_defining,
                    "star_allele": v.star_allele,
                    "functional_impact": v.functional_impact
                }
                for v in diplotype_call.variants
            ]
        }
        calls.append(call_entry)

    # Build output
    output = {
        "sample_id": "${sample_id}",
        "reference_build": detected_build,
        "gene_tiers": tiers,
        "calls": calls
    }

    with open("alleles.json", "w") as f:
        json.dump(output, f, indent=2)

    print(f"Called {len(calls)} genes for sample ${sample_id}")
    """
}
