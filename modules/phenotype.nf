/*
 * Module: Phenotype Translation
 *
 * Translates diplotype calls to phenotypes and prescribing implications.
 */

process PHENOTYPE_TRANSLATE {
    tag "$sample_id"
    label 'process_low'

    publishDir "${params.outdir}/phenotypes", mode: params.publish_dir_mode

    input:
    tuple val(sample_id), path(alleles_json)
    path phenotype_maps_dir
    path implications_dir

    output:
    tuple val(sample_id), path("phenotypes.json"), emit: phenotypes

    script:
    """
    #!/usr/bin/env python3

    import json
    from pathlib import Path

    from elpgx.allele_caller import DiplotypeCall
    from elpgx.phenotype import translate_diplotype_to_phenotype
    from elpgx.models import Variant, ConfidenceMetrics

    # Load alleles.json
    with open("${alleles_json}") as f:
        alleles_data = json.load(f)

    phenotype_maps_path = Path("${phenotype_maps_dir}")
    implications_path = Path("${implications_dir}")

    phenotypes = []

    for call_data in alleles_data.get("calls", []):
        # Reconstruct Variant objects
        variants = []
        for v in call_data.get("variants", []):
            variants.append(Variant(
                rsid=v.get("rsid"),
                chrom=v["chrom"],
                pos=v["pos"],
                ref=v["ref"],
                alt=v["alt"],
                genotype=v["genotype"],
                quality=v.get("quality", 0),
                depth=v.get("depth", 0),
                allele_defining=True
            ))

        # Reconstruct confidence
        conf_data = call_data.get("confidence", {})
        confidence = ConfidenceMetrics(
            score=conf_data.get("score", 0.95),
            flags=[]
        )

        # Create DiplotypeCall
        diplotype_call = DiplotypeCall(
            gene=call_data["gene"],
            allele1=call_data["haplotypes"][0],
            allele2=call_data["haplotypes"][1],
            diplotype=call_data["diplotype"],
            variants=variants,
            allele1_function=call_data.get("allele1_function", "Unknown"),
            allele2_function=call_data.get("allele2_function", "Unknown"),
            allele1_activity=call_data.get("allele1_activity"),
            allele2_activity=call_data.get("allele2_activity"),
            call_quality=call_data.get("call_quality", "definitive"),
            confidence=confidence
        )

        # Translate to phenotype
        result = translate_diplotype_to_phenotype(
            diplotype_call,
            phenotype_maps_path,
            implications_path
        )

        phenotypes.append({
            "gene": result.gene,
            "diplotype": result.diplotype,
            "phenotype": result.phenotype,
            "activity_score": result.activity_score,
            "functional_status": result.functional_status,
            "prescribing_implication": result.prescribing_implication
        })

    output = {
        "sample_id": "${sample_id}",
        "phenotypes": phenotypes
    }

    with open("phenotypes.json", "w") as f:
        json.dump(output, f, indent=2)

    print(f"Translated {len(phenotypes)} phenotypes for sample ${sample_id}")
    """
}
