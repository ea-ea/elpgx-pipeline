#!/usr/bin/env nextflow

/*
 * ElPGx Pipeline - Clinical Pharmacogenetics Analysis
 *
 * Clinical-grade pipeline for pharmacogenetic analysis from WGS/WES data.
 * Predicts drug response phenotypes for clinically relevant pharmacogenes.
 *
 * Author: Elif Arslan
 * Input:  Aligned BAM/CRAM (hg38)
 * Output: Structured JSON report with diplotypes, phenotypes, and confidence
 */

nextflow.enable.dsl = 2

/*
 * Helper functions for provenance collection
 */
def getGitCommit() {
    try {
        def proc = "git rev-parse HEAD".execute()
        proc.waitFor()
        return proc.text.trim()
    } catch (Exception e) {
        return "unknown"
    }
}

def getReferenceMd5(reference) {
    // Try to read from .md5 sidecar file first (standard GATK bundle convention)
    def md5File = file("${reference}.md5")
    if (md5File.exists()) {
        return md5File.text.trim().split()[0]
    }
    // Compute MD5 if no sidecar exists
    try {
        def proc = ["md5sum", reference.toString()].execute()
        proc.waitFor()
        if (proc.exitValue() == 0) {
            return proc.text.trim().split()[0]
        }
    } catch (Exception e) {
        // md5sum not available, try md5 (macOS)
        try {
            def proc2 = ["md5", "-q", reference.toString()].execute()
            proc2.waitFor()
            if (proc2.exitValue() == 0) {
                return proc2.text.trim()
            }
        } catch (Exception e2) {}
    }
    return "unknown"
}

def getInputBamMd5(bam_path) {
    // Try to read from .md5 sidecar file first
    def md5File = file("${bam_path}.md5")
    if (md5File.exists()) {
        return md5File.text.trim().split()[0]
    }
    // Compute MD5 for the input BAM/CRAM
    try {
        def proc = ["md5sum", bam_path.toString()].execute()
        proc.waitFor()
        if (proc.exitValue() == 0) {
            return proc.text.trim().split()[0]
        }
    } catch (Exception e) {
        try {
            def proc2 = ["md5", "-q", bam_path.toString()].execute()
            proc2.waitFor()
            if (proc2.exitValue() == 0) {
                return proc2.text.trim()
            }
        } catch (Exception e2) {}
    }
    return "unknown"
}

// Import modules
include { VALIDATE_INPUT }     from './modules/validate'
include { VARIANT_CALL }       from './modules/variant_call'
include { ALLELE_CALL }        from './modules/allele_call'
include { HLA_CALL }           from './modules/hla_call'
include { PHENOTYPE_TRANSLATE } from './modules/phenotype'
include { GENERATE_REPORT }    from './modules/report'

/*
 * Print pipeline header
 */
log.info """
╔═══════════════════════════════════════════════════════════════════╗
║              ElPGx Pipeline v${workflow.manifest.version}                            ║
║         Clinical Pharmacogenetics Analysis Pipeline               ║
║                     by Elif Arslan                                ║
╚═══════════════════════════════════════════════════════════════════╝

Input BAM/CRAM    : ${params.input}
Reference         : ${params.reference}
Gene Tiers        : ${params.gene_tiers}
HLA Module        : ${params.enable_hla}
Output Directory  : ${params.outdir}
"""

/*
 * Validate required parameters
 */
if (!params.input) {
    error "Input BAM/CRAM file not specified. Use --input <path>"
}

if (!params.reference) {
    error "Reference FASTA not specified. Use --reference <path>"
}

/*
 * Main workflow
 */
workflow {

    // Create input channel
    // Extract sample ID from filename if not in read group
    bam_file = file(params.input)
    sample_id = bam_file.baseName.replaceAll(/\.(bam|cram)$/, '')

    // Check for index file
    if (bam_file.extension == 'bam') {
        bai_file = file("${params.input}.bai")
        if (!bai_file.exists()) {
            bai_file = file(params.input.replaceAll(/\.bam$/, '.bai'))
        }
    } else {
        bai_file = file("${params.input}.crai")
        if (!bai_file.exists()) {
            bai_file = file(params.input.replaceAll(/\.cram$/, '.crai'))
        }
    }

    if (!bai_file.exists()) {
        error "Index file not found for ${params.input}"
    }

    // Create input tuple channel
    input_ch = Channel.of(tuple(sample_id, bam_file, bai_file))

    // Reference files
    reference    = file(params.reference)
    reference_fai = file("${params.reference}.fai")
    reference_dict = file(params.reference.replaceAll(/\.fa(sta)?$/, '.dict'))

    // Data files
    pgx_regions     = file(params.pgx_regions)
    pgx_regions_grch37 = file("${projectDir}/data/pgx_regions_grch37.bed")
    allele_defs     = file(params.pgx_db)
    phenotype_maps  = file(params.phenotype_maps)
    implications    = file(params.implications)
    gene_tiers_file = file(params.gene_tiers_file)
    schema_path     = file(params.schema_path)
    gene_strategies = file("${projectDir}/data/gene_strategies")

    // Provenance information
    git_commit = getGitCommit()
    nextflow_ver = workflow.nextflow.version
    reference_md5 = getReferenceMd5(reference)
    reference_path = params.reference
    input_bam_md5 = getInputBamMd5(bam_file)

    /*
     * Stage 1: Input Validation
     */
    VALIDATE_INPUT(
        input_ch,
        reference,
        reference_fai,
        pgx_regions
    )

    /*
     * Stage 2: Variant Calling
     */
    VARIANT_CALL(
        VALIDATE_INPUT.out.validated_bam,
        reference,
        reference_fai,
        reference_dict,
        pgx_regions,
        pgx_regions_grch37,
        VALIDATE_INPUT.out.validation_json
    )

    /*
     * Stage 3: Allele Calling
     */
    // Pass coordinate map for GRCh37 translation
    grch37_coord_map = file("${projectDir}/data/grch37_coordinate_map.json")

    ALLELE_CALL(
        VARIANT_CALL.out.vcf,
        allele_defs,
        params.gene_tiers,
        VALIDATE_INPUT.out.validation_json,
        grch37_coord_map
    )

    /*
     * Stage 3b: HLA Calling (Optional)
     */
    if (params.enable_hla) {
        HLA_CALL(
            VALIDATE_INPUT.out.validated_bam,
            reference
        )
    }

    /*
     * Stage 4: Phenotype Translation
     */
    PHENOTYPE_TRANSLATE(
        ALLELE_CALL.out.alleles,
        phenotype_maps,
        implications
    )

    /*
     * Stage 5: Report Generation
     */
    // Join alleles and phenotypes by sample_id
    report_input = ALLELE_CALL.out.alleles
        .join(PHENOTYPE_TRANSLATE.out.phenotypes)
        .map { sample_id, alleles, phenotypes ->
            tuple(sample_id, alleles, phenotypes)
        }

    GENERATE_REPORT(
        report_input,
        gene_tiers_file,
        schema_path,
        gene_strategies,
        params.gene_tiers,
        params.enable_hla,
        workflow.manifest.version,
        params.pgx_db_version,
        git_commit,
        nextflow_ver,
        reference_md5,
        reference_path,
        input_bam_md5
    )
}

/*
 * Workflow completion handler
 */
workflow.onComplete {
    log.info """
═══════════════════════════════════════════════════════════════════
ElPGx Pipeline completed at: ${workflow.complete}
Duration                   : ${workflow.duration}
Success                    : ${workflow.success}
Exit status                : ${workflow.exitStatus}
Output directory           : ${params.outdir}
═══════════════════════════════════════════════════════════════════
    """
}

workflow.onError {
    log.error "ElPGx Pipeline execution stopped with error: ${workflow.errorMessage}"
}
