/*
 * Module: Variant Calling
 *
 * Calls variants in PGx target regions using GATK HaplotypeCaller.
 * Automatically selects GRCh37 or GRCh38 BED based on validation result.
 */

process VARIANT_CALL {
    tag "$sample_id"
    label 'process_high'

    publishDir "${params.outdir}/variants", mode: params.publish_dir_mode

    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference
    path reference_fai
    path reference_dict
    path pgx_regions
    path pgx_regions_grch37
    path validation_json

    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz"), path("${sample_id}.vcf.gz.tbi"), emit: vcf
    path "variants.meta.json", emit: meta

    script:
    """
    #!/bin/bash
    set -euo pipefail

    # Detect build from validation result
    BUILD=\$(python3 -c "import json; print(json.load(open('${validation_json}'))['reference_build'])" 2>/dev/null || echo "GRCh38")

    # Select correct BED file
    if [ "\$BUILD" = "GRCh37" ]; then
        REGIONS="${pgx_regions_grch37}"
    else
        REGIONS="${pgx_regions}"
    fi

    # Run GATK HaplotypeCaller on PGx regions
    gatk HaplotypeCaller \\
        -R ${reference} \\
        -I ${bam} \\
        -L \$REGIONS \\
        -O ${sample_id}.vcf.gz \\
        --min-base-quality-score ${params.min_base_quality} \\
        --minimum-mapping-quality ${params.min_mapping_quality} \\
        --disable-sequence-dictionary-validation true

    # Index VCF
    gatk IndexFeatureFile -I ${sample_id}.vcf.gz

    # Create metadata file
    cat > variants.meta.json << EOF
    {
        "caller": "GATK HaplotypeCaller",
        "caller_version": "\$(gatk --version | head -1 | cut -d' ' -f2)",
        "reference_build": "\$BUILD",
        "sample_id": "${sample_id}",
        "regions_bed": "\$REGIONS",
        "min_base_quality": ${params.min_base_quality},
        "min_mapping_quality": ${params.min_mapping_quality}
    }
EOF
    """
}
