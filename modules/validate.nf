/*
 * Module: Input Validation
 *
 * Validates BAM/CRAM input files and extracts sample metadata.
 */

process VALIDATE_INPUT {
    tag "$sample_id"
    label 'process_low'

    publishDir "${params.outdir}/validation", mode: params.publish_dir_mode

    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference
    path reference_fai
    path pgx_regions

    output:
    tuple val(sample_id), path(bam), path(bai), emit: validated_bam
    path "validation_result.json", emit: validation_json

    script:
    """
    #!/usr/bin/env python3

    import json
    import sys
    from pathlib import Path

    from elpgx.validate import validate_bam_file

    bam_path = Path("${bam}")
    ref_path = Path("${reference}")

    result = validate_bam_file(bam_path, ref_path)

    output = {
        "sample_id": result.sample_id or "${sample_id}",
        "is_valid": result.is_valid,
        "reference_build": result.reference_build,
        "errors": result.errors,
        "warnings": result.warnings
    }

    with open("validation_result.json", "w") as f:
        json.dump(output, f, indent=2)

    if not result.is_valid:
        print("Validation failed:", result.errors, file=sys.stderr)
        sys.exit(1)

    print(f"Validation passed: {result.sample_id}, {result.reference_build}")
    """
}
