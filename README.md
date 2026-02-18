# ElPGx - Clinical Pharmacogenetics Pipeline

Clinical-grade Pharmacogenetics (PGx) NGS Pipeline for predicting drug response phenotypes from WGS/WES data.

**Version:** 1.1.0
**Author:** Elif Arslan
**License:** MIT
**PGx Database:** CPIC v1.24

## Overview

ElPGx analyzes aligned BAM/CRAM sequencing data to identify pharmacogenomic star alleles, predict drug-response phenotypes, and generate structured clinical reports. It supports both GRCh38 and GRCh37 reference builds with automatic detection.

## Input / Output

### Input
- **Aligned BAM/CRAM file** (coordinate-sorted, indexed)
- **Reference genome:** GRCh38 or GRCh37 FASTA with `.fai` index and `.dict` dictionary in the same directory

### Output
A structured **JSON report** (schema v1.1.0) containing:

1. **Star-allele diplotypes** (e.g., `*1/*17`)
2. **Clinical phenotypes** (e.g., "Rapid Metabolizer")
3. **Activity scores** and functional status
4. **Prescribing implications** (CPIC-based)
5. **Confidence scores** (0.0-1.0) with quality flags
6. **Full provenance** (git commit, tool versions, reference checksums)

<details>
<summary>Example JSON output</summary>

```json
{
  "schema_version": "1.1.0",
  "meta": {
    "sample_id": "SAMPLE_001",
    "pipeline_version": "1.1.0",
    "run_date": "2026-02-18T12:00:00Z",
    "reference_build": "GRCh38",
    "pgx_db_version": "CPIC v1.24",
    "gene_tiers_included": ["group1", "group2"],
    "hla_module_enabled": false,
    "provenance": {
      "git_commit": "14ec7af...",
      "nextflow_version": "25.10.3",
      "variant_caller": "GATK HaplotypeCaller",
      "variant_caller_version": "4.5.0.0"
    }
  },
  "results": [
    {
      "gene": "CYP2C19",
      "gene_tier": "group1",
      "biomarker_type": "metabolic_enzyme",
      "diplotype": "*1/*17",
      "functional_status": "Normal Function / Increased Function",
      "phenotype": "Rapid Metabolizer",
      "activity_score": 1.5,
      "prescribing_implication": "Consider alternative to clopidogrel...",
      "confidence": {
        "score": 0.95,
        "flags": []
      },
      "variants": [...]
    }
  ],
  "warnings": []
}
```

</details>

## Pipeline Architecture

```
BAM/CRAM --> Validate --> Variant Call --> Allele Call --> Phenotype --> Report
               |           (GATK)         (Python)        (CPIC)      (JSON)
               |                              |
          Build detect                   GRCh37 coord
          (GRCh38/37)                    translation
```

| Stage | Tool | Description |
|-------|------|-------------|
| 1. Validation | Python/pysam | Verify BAM integrity, detect reference build (GRCh38/GRCh37) |
| 2. Variant Calling | GATK HaplotypeCaller | Call variants in PGx target regions |
| 3. Allele Calling | Python | Match variants to star alleles (with GRCh37 coordinate translation) |
| 3b. HLA Calling | arcasHLA (optional) | HLA typing for Group 4 genes |
| 4. Phenotype | Python | Translate diplotype to phenotype using CPIC guidelines |
| 5. Report | Python/Pydantic | Generate validated JSON report with provenance |

## Supported Genes (17 genes)

| Group | Genes | Biomarker Type |
|-------|-------|----------------|
| **Group 1** | CYP2C19, CYP2C9, DPYD, NUDT15, SLCO1B1, TPMT, UGT1A1, VKORC1 | Metabolic enzymes, transporter, drug target |
| **Group 2** | CYP2B6, CYP2D6, CYP3A5, G6PD | Metabolic enzymes |
| **Group 3** | CACNA1S, IFNL3, NAT2, RYR1 | Drug targets, metabolic enzyme |
| **Group 4** | HLA-A, HLA-B | HLA markers (requires `--enable_hla`) |

## Quick Start

### Prerequisites

- Python 3.10+
- Reference genome FASTA with `.fai` index and `.dict` dictionary
- For Nextflow: [Nextflow](https://nextflow.io) 21.04+ and Docker

### Option 1: Command-Line Interface (Local)

```bash
# Install from source
git clone https://github.com/ea-ea/elpgx-pipeline.git
cd elpgx
pip install -e .

# Full pipeline from BAM/CRAM (GRCh38 or GRCh37 auto-detected)
elpgx run \
  --bam sample.bam \
  --reference hg38.fa \
  --output report.json \
  --gene-tiers group1,group2,group3

# From pre-called VCF
elpgx analyze \
  --vcf sample.vcf.gz \
  --output report.json

# Validate input file
elpgx validate --bam sample.bam --reference hg38.fa

# Show version
elpgx version
```

### Option 2: Docker

```bash
# Pull the image
docker pull elifffarslan/elpgx-python:1.1.0

# Run pipeline (mount your data directory)
docker run -v /path/to/data:/data elifffarslan/elpgx-python:1.1.0 \
  run \
  --bam /data/sample.bam \
  --reference /data/reference.fa \
  --output /data/report.json \
  --gene-tiers group1,group2,group3
```

### Option 3: Nextflow Pipeline (Recommended for Production)

```bash
# Single sample (GRCh38)
nextflow run main.nf \
  --input sample.bam \
  --reference GRCh38.fa \
  --gene_tiers group1,group2,group3 \
  --outdir results/ \
  -profile docker

# Single sample (GRCh37 - auto-detected)
nextflow run main.nf \
  --input sample_grch37.bam \
  --reference human_g1k_v37.fasta \
  --gene_tiers group1,group2,group3 \
  --outdir results/ \
  -profile docker

# CRAM input (also supported)
nextflow run main.nf \
  --input sample.cram \
  --reference GRCh38.fa \
  --outdir results/ \
  -profile docker
```

### Option 4: Desktop Application (GUI)

```bash
pip install -e ".[gui]"
elpgx-gui
```

## Configuration

### CLI Options

| Option | Description |
|--------|-------------|
| `--bam` / `-b` | Input BAM/CRAM file |
| `--reference` / `-r` | Reference FASTA with `.fai` index |
| `--output` / `-o` | Output JSON report path |
| `--gene-tiers` / `-g` | Gene tiers: `group1`, `group1,group2,group3`, etc. |
| `--enable-hla` | Enable HLA typing module |
| `--data-dir` / `-d` | Custom PGx data directory |

### Nextflow Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | required | BAM/CRAM file path |
| `--reference` | required | Reference FASTA (GRCh38 or GRCh37) |
| `--gene_tiers` | `group1` | Comma-separated gene tiers |
| `--enable_hla` | `false` | Enable HLA typing |
| `--outdir` | `./results` | Output directory |
| `--min_base_quality` | `20` | Minimum base quality for variant calling |
| `--min_mapping_quality` | `20` | Minimum mapping quality |

### Nextflow Profiles

| Profile | Description |
|---------|-------------|
| `docker` | Run with Docker containers (recommended) |
| `singularity` | Run with Singularity containers |
| `local` | Run with local Python installation |
| `awsbatch` | Run on AWS Batch |
| `sevenbridges` | Run on Seven Bridges platform |

## Project Structure

```
elpgx/
├── main.nf                 # Nextflow entry point
├── nextflow.config         # Pipeline configuration
├── modules/                # Nextflow process modules
│   ├── validate.nf         # Input validation & build detection
│   ├── variant_call.nf     # GATK HaplotypeCaller
│   ├── allele_call.nf      # Star-allele calling
│   ├── phenotype.nf        # Phenotype translation
│   ├── report.nf           # Report generation
│   └── hla_call.nf         # HLA typing (optional)
├── elpgx/                  # Python package
│   ├── cli.py              # Command-line interface
│   ├── gui.py              # Desktop GUI (tkinter)
│   ├── models.py           # Pydantic data models
│   ├── validate.py         # Input validation (GRCh38/GRCh37)
│   ├── vcf_parser.py       # VCF parsing & GRCh37 translation
│   ├── allele_caller.py    # Star-allele matching
│   ├── phenotype.py        # Phenotype translation (CPIC)
│   ├── confidence.py       # Confidence scoring
│   └── report.py           # Report generation
├── data/                   # PGx reference data (CPIC v1.24)
│   ├── allele_definitions/ # Star-allele JSON definitions (17 genes)
│   ├── phenotype_maps/     # Diplotype-to-phenotype tables
│   ├── prescribing_implications/
│   ├── gene_strategies/    # Gene-specific calling strategies
│   ├── gene_tiers.json     # Gene tier configuration
│   ├── pgx_regions.bed     # GRCh38 target regions
│   ├── pgx_regions_grch37.bed  # GRCh37 target regions
│   └── grch37_coordinate_map.json
├── docker/                 # Dockerfiles
├── schemas/                # JSON report schema
├── tests/                  # Test suite (202 tests)
└── docs/                   # Documentation
```

## Testing

```bash
# Install with dev dependencies
pip install -e ".[dev]"

# Run all tests
pytest tests/

# Run with coverage
pytest tests/ --cov=elpgx

# Run specific test module
pytest tests/unit/test_allele_caller.py -v
```

## Docker Images

| Image | Description |
|-------|-------------|
| `elifffarslan/elpgx-python:1.1.0` | Main pipeline image (Python + GATK + samtools) |
| `broadinstitute/gatk:4.5.0.0` | GATK variant caller (used by Nextflow) |

## Key Features

- **Dual-build support**: Automatic GRCh38/GRCh37 detection with coordinate translation
- **BAM and CRAM**: Full support for both alignment formats (including alt-contig CRAMs)
- **Clinical-grade confidence**: Multi-factor scoring with depth, quality, and phasing assessment
- **Phenotype blocking**: Strategy-based rules prevent reporting uncertain phenotypes (e.g., CYP2D6)
- **Full provenance**: Git commit, tool versions, reference checksums tracked in every report
- **Schema validation**: Pydantic models ensure report structure integrity
- **Scalable**: Nextflow enables HPC/cloud execution (AWS Batch, Seven Bridges)

## References

- [CPIC](https://cpicpgx.org) - Clinical Pharmacogenetics Implementation Consortium
- [PharmVar](https://www.pharmvar.org) - Pharmacogene Variation Consortium
- [GATK](https://gatk.broadinstitute.org) - Genome Analysis Toolkit

## Usage Notice

This repository is provided solely for technical assignment evaluation purposes.
All rights reserved. No permission is granted to copy, modify, distribute, or
reuse this code outside the scope of the evaluation without explicit written
permission from the author.
