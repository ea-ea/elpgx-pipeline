# ElPGx - Complete Project Documentation

## Clinical Pharmacogenetics Analysis Pipeline

**Version:** 1.1.0
**Author:** Elif Arslan
**Repository:** https://github.com/ea-ea/elpgx
**CPIC Version:** v1.24 (February 2026)
**Reference Build:** GRCh38

---

## Table of Contents

1. [Project Overview](#1-project-overview)
2. [Architecture](#2-architecture)
3. [Installation Methods](#3-installation-methods)
4. [Tools & Interfaces](#4-tools--interfaces)
5. [Python Modules](#5-python-modules)
6. [Nextflow Pipeline](#6-nextflow-pipeline)
7. [Docker Images](#7-docker-images)
8. [Data Files](#8-data-files)
9. [Configuration Files](#9-configuration-files)
10. [Usage Examples](#10-usage-examples)
11. [Output Format](#11-output-format)
12. [Troubleshooting](#12-troubleshooting)

---

## 1. Project Overview

### What is ElPGx?

ElPGx is a **Clinical Pharmacogenetics Analysis Pipeline** that processes whole-genome sequencing (WGS) or whole-exome sequencing (WES) data to predict drug response phenotypes based on genetic variants.

### What it Does

```
BAM/CRAM file (aligned to hg38)
    │
    ▼
┌─────────────────────────────────────────────────────────────┐
│                     ElPGx Pipeline                          │
│                                                             │
│  1. Validate Input  →  Check BAM integrity, hg38 alignment │
│  2. Call Variants   →  GATK HaplotypeCaller in PGx regions │
│  3. Call Alleles    →  Match variants to star alleles      │
│  4. Phenotype       →  Diplotype → Activity → Phenotype    │
│  5. Report          →  JSON with clinical implications     │
│                                                             │
└─────────────────────────────────────────────────────────────┘
    │
    ▼
JSON Report with:
  - Star-allele diplotypes (*1/*17)
  - Phenotypes (Rapid Metabolizer)
  - Prescribing implications
  - Confidence scores
```

### Supported Genes

| Tier | Genes | Description |
|------|-------|-------------|
| **Group 1** | CYP2C19, CYP2C9, VKORC1, SLCO1B1, TPMT, DPYD, UGT1A1, NUDT15 | Core panel - CPIC Level A |
| **Group 2** | CYP2D6 | Complex CNV analysis (future) |
| **Group 3** | CYP3A5, NAT2, IFNL3 | Extended panel |
| **Group 4** | HLA-A, HLA-B | HLA typing (optional) |

---

## 2. Architecture

### Directory Structure

```
PGx_Genomize_ElifArslan/
│
├── elpgx/                      # Python package (core logic)
│   ├── __init__.py             # Package initialization, version
│   ├── models.py               # Pydantic data models
│   ├── validate.py             # Input validation
│   ├── vcf_parser.py           # VCF parsing + GATK calling
│   ├── allele_caller.py        # Star-allele matching
│   ├── phenotype.py            # Phenotype translation
│   ├── confidence.py           # Confidence scoring
│   ├── report.py               # JSON report generation
│   ├── cli.py                  # Command-line interface
│   ├── gui.py                  # Graphical interface
│   └── hla_parser.py           # HLA parsing (placeholder)
│
├── data/                       # Reference data
│   ├── allele_definitions/     # Gene → star allele definitions
│   ├── phenotype_maps/         # Diplotype → phenotype mappings
│   ├── prescribing_implications/  # Clinical guidance text
│   ├── pgx_regions.bed         # Target genomic regions
│   └── gene_tiers.json         # Gene classification
│
├── modules/                    # Nextflow DSL2 modules
│   ├── validate.nf             # Input validation process
│   ├── variant_call.nf         # GATK variant calling
│   ├── allele_call.nf          # Star-allele calling
│   ├── phenotype.nf            # Phenotype translation
│   ├── report.nf               # Report generation
│   └── hla_call.nf             # HLA typing (optional)
│
├── docker/                     # Docker images
│   ├── elpgx-python/Dockerfile # Python analysis image
│   ├── gatk/Dockerfile         # GATK image
│   └── hla/Dockerfile          # HLA typing image
│
├── schemas/                    # JSON schemas
│   └── elpgx_report_schema.json
│
├── tests/                      # Unit tests
│   └── unit/
│
├── scripts/                    # Utility scripts
│   └── build_docker.sh         # Docker build script
│
├── build/                      # PyInstaller specs
│   └── elpgx.spec
│
├── dist/                       # Built applications
│   └── ElPGx.app               # Standalone macOS app
│
├── main.nf                     # Nextflow entry point
├── nextflow.config             # Nextflow configuration
├── docker-compose.yml          # Docker compose file
├── pyproject.toml              # Python project config
└── README.md                   # Project readme
```

---

## 3. Installation Methods

### Method A: Python Package (Development)

```bash
# Clone repository
git clone https://github.com/ea-ea/elpgx-pipeline.git
cd elpgx

# Install in editable mode
pip install -e .

# Run CLI
python -m elpgx.cli run --help
```

### Method B: Standalone App (End Users)

```bash
# Use pre-built app (macOS)
open dist/ElPGx.app

# Or build from source
pip install pyinstaller
pyinstaller build/elpgx.spec --clean
```

### Method C: Nextflow + Docker (Production)

```bash
# Install Docker Desktop
brew install --cask docker

# Install Nextflow
curl -s https://get.nextflow.io | bash

# Build Docker images
./scripts/build_docker.sh --push

# Run pipeline
nextflow run main.nf -profile docker --input sample.bam --reference hg38.fa
```

---

## 4. Tools & Interfaces

### 4.1 Command-Line Interface (CLI)

**Location:** `elpgx/cli.py`

**Why:** For scripting, automation, and command-line users.

**Commands:**

| Command | Purpose | Example |
|---------|---------|---------|
| `run` | Full pipeline from BAM/CRAM | `python -m elpgx.cli run -b sample.bam -r hg38.fa` |
| `analyze` | Analyze pre-called VCF | `python -m elpgx.cli analyze -v variants.vcf` |
| `validate` | Validate input files | `python -m elpgx.cli validate -b sample.bam -r hg38.fa` |
| `version` | Show version info | `python -m elpgx.cli version` |

**Full Pipeline Example:**
```bash
python -m elpgx.cli run \
  -b ~/Downloads/sample.cram \
  -r ~/Downloads/hg38.fa \
  -g group1,group2,group3 \
  --pretty
```

**Options:**
```
-b, --bam PATH          Input BAM/CRAM file (required)
-r, --reference PATH    Reference FASTA (required)
-o, --output PATH       Output JSON (default: derived from input)
-g, --gene-tiers TEXT   Gene tiers: group1,group2,group3 (default: group1)
-s, --sample-id TEXT    Sample ID (default: from BAM header)
--enable-hla            Enable HLA typing module
--pretty                Pretty-print JSON output
```

---

### 4.2 Graphical User Interface (GUI)

**Location:** `elpgx/gui.py`

**Why:** For users who prefer visual interfaces, lab technicians.

**How to Run:**
```bash
# From Python
python -m elpgx.gui

# Or standalone app
open dist/ElPGx.app
```

**Features:**
- File browser for BAM/CRAM and reference selection
- Gene tier checkboxes
- Real-time progress display
- Results viewer
- Export to JSON

---

### 4.3 Standalone Application

**Location:** `dist/ElPGx.app`

**Why:** No Python installation required, for end users.

**How to Build:**
```bash
pip install pyinstaller
pyinstaller build/elpgx.spec --clean --noconfirm
```

**Contents:**
- Python runtime (bundled)
- GATK JAR (406MB, bundled)
- All dependencies
- Data files

---

### 4.4 Nextflow Pipeline

**Location:** `main.nf`

**Why:** Production deployment, cloud platforms (Seven Bridges, AWS), reproducibility.

**How to Run:**
```bash
nextflow run main.nf \
  --input sample.bam \
  --reference hg38.fa \
  --gene_tiers group1,group2,group3 \
  -profile docker
```

**Parameters:**
```
--input             BAM/CRAM file path
--reference         Reference FASTA path
--gene_tiers        Comma-separated tiers (default: group1)
--enable_hla        Enable HLA module (default: false)
--outdir            Output directory (default: ./results)
```

**Profiles:**
```
-profile docker         Use Docker containers
-profile singularity    Use Singularity containers
-profile sevenbridges   Seven Bridges platform
-profile awsbatch       AWS Batch execution
-profile local          Local execution
-profile test           Run with test data
```

---

## 5. Python Modules

### 5.1 models.py

**Purpose:** Define data structures for the entire pipeline.

**Key Classes:**

| Class | Purpose |
|-------|---------|
| `Variant` | Single genetic variant (rsid, chrom, pos, ref, alt, genotype) |
| `GeneResult` | Complete result for one gene (diplotype, phenotype, confidence) |
| `PgxReport` | Full report container with metadata |
| `PipelineMeta` | Pipeline provenance (version, date, settings) |
| `ConfidenceMetrics` | Quality metrics (score, flags) |

**Example:**
```python
from elpgx.models import Variant, GeneResult

variant = Variant(
    rsid="rs12248560",
    chrom="chr10",
    pos=94762706,
    ref="C",
    alt="T",
    genotype="0/1",
    quality=38.5,
    depth=42
)
```

---

### 5.2 validate.py

**Purpose:** Validate BAM/CRAM files before analysis.

**Key Function:**
```python
def validate_bam_file(bam_path: Path, reference_path: Path) -> ValidationResult
```

**What it Checks:**
- File exists and is readable
- BAM/CRAM has valid index (.bai/.crai)
- Reference is GRCh38/hg38
- Chromosome naming convention (chr1 vs 1)

**Example:**
```python
from elpgx.validate import validate_bam_file

result = validate_bam_file(Path("sample.bam"), Path("hg38.fa"))
if result.is_valid:
    print(f"Sample: {result.sample_id}, Build: {result.reference_build}")
else:
    print(f"Errors: {result.errors}")
```

---

### 5.3 vcf_parser.py

**Purpose:** Parse VCF files and call variants from BAM.

**Key Functions:**

| Function | Purpose |
|----------|---------|
| `call_variants_gatk()` | Call variants using GATK HaplotypeCaller |
| `call_variants_bcftools()` | Fallback variant calling with bcftools |
| `parse_vcf()` | Parse VCF file into variant objects |
| `filter_by_positions()` | Filter variants to PGx target positions |

**Example:**
```python
from elpgx.vcf_parser import call_variants_gatk

vcf_path = call_variants_gatk(
    bam_path=Path("sample.bam"),
    reference_path=Path("hg38.fa"),
    regions_bed=Path("pgx_regions.bed"),
    output_vcf=Path("variants.vcf.gz")
)
```

---

### 5.4 allele_caller.py

**Purpose:** Match variants to star alleles and determine diplotypes.

**Key Functions:**

| Function | Purpose |
|----------|---------|
| `load_gene_definitions()` | Load allele definitions from JSON |
| `determine_diplotype()` | Call diplotype for one gene |
| `DiplotypeCall` | Result dataclass with alleles, variants, confidence |

**Example:**
```python
from elpgx.allele_caller import load_gene_definitions, determine_diplotype

gene_defs = load_gene_definitions(Path("data/allele_definitions/CYP2C19.json"))
call = determine_diplotype(variants, gene_defs)
print(f"{call.gene}: {call.diplotype}")  # CYP2C19: *1/*17
```

---

### 5.5 phenotype.py

**Purpose:** Translate diplotypes to phenotypes and prescribing implications.

**Key Function:**
```python
def translate_diplotype_to_phenotype(
    diplotype_call: DiplotypeCall,
    phenotype_maps_dir: Path,
    implications_dir: Path
) -> PhenotypeResult
```

**Phenotype Categories:**
- Poor Metabolizer
- Intermediate Metabolizer
- Normal Metabolizer
- Rapid Metabolizer
- Ultrarapid Metabolizer

**Example:**
```python
from elpgx.phenotype import translate_diplotype_to_phenotype

result = translate_diplotype_to_phenotype(
    diplotype_call,
    Path("data/phenotype_maps"),
    Path("data/prescribing_implications")
)
print(f"Phenotype: {result.phenotype}")
print(f"Implication: {result.prescribing_implication}")
```

---

### 5.6 confidence.py

**Purpose:** Calculate confidence scores for variant calls.

**Formula:**
```
confidence = 0.40 * depth_score + 0.35 * quality_score + 0.25 * completeness_score
```

**Flags:**
| Flag | Trigger |
|------|---------|
| `LOW_DEPTH` | Any position < 10× depth |
| `LOW_QUALITY` | Any variant QUAL < 20 |
| `INCOMPLETE_ALLELE` | Missing defining positions |
| `AMBIGUOUS_CALL` | Multiple possible diplotypes |

---

### 5.7 report.py

**Purpose:** Generate final JSON report.

**Key Function:**
```python
def generate_report(
    gene_results: list[GeneResult],
    sample_id: str,
    reference_build: str = "GRCh38"
) -> PgxReport
```

---

## 6. Nextflow Pipeline

### 6.1 main.nf

**Purpose:** Pipeline entry point, orchestrates all modules.

**Workflow:**
```
VALIDATE_INPUT → VARIANT_CALL → ALLELE_CALL → PHENOTYPE_TRANSLATE → GENERATE_REPORT
                                     │
                                     └── HLA_CALL (optional)
```

---

### 6.2 Module: validate.nf

**Process:** `VALIDATE_INPUT`

**Purpose:** Validate BAM/CRAM files.

**Container:** `elifffarslan/elpgx-python:1.1.0`

**Input:**
- BAM/CRAM file with index
- Reference FASTA

**Output:**
- Validated BAM tuple
- validation_result.json

---

### 6.3 Module: variant_call.nf

**Process:** `VARIANT_CALL`

**Purpose:** Call variants in PGx regions using GATK.

**Container:** `broadinstitute/gatk:4.5.0.0`

**Command:**
```bash
gatk HaplotypeCaller \
    -R reference.fa \
    -I sample.bam \
    -L pgx_regions.bed \
    -O sample.vcf.gz
```

**Output:**
- VCF file with variants
- VCF index

---

### 6.4 Module: allele_call.nf

**Process:** `ALLELE_CALL`

**Purpose:** Determine star-allele diplotypes from variants.

**Container:** `elifffarslan/elpgx-python:1.1.0`

**Output:**
- alleles.json with diplotype calls

---

### 6.5 Module: phenotype.nf

**Process:** `PHENOTYPE_TRANSLATE`

**Purpose:** Convert diplotypes to phenotypes.

**Container:** `elifffarslan/elpgx-python:1.1.0`

**Output:**
- phenotypes.json

---

### 6.6 Module: report.nf

**Process:** `GENERATE_REPORT`

**Purpose:** Generate final clinical report.

**Container:** `elifffarslan/elpgx-python:1.1.0`

**Output:**
- {sample_id}_pgx_report.json

---

### 6.7 Module: hla_call.nf

**Process:** `HLA_CALL`

**Purpose:** HLA typing (optional, placeholder).

**Container:** `elifffarslan/elpgx-hla:1.1.0`

**Activation:** Only when `--enable_hla` is set.

---

## 7. Docker Images

### 7.1 elpgx-python

**Dockerfile:** `docker/elpgx-python/Dockerfile`

**Purpose:** Python analysis environment.

**Contents:**
- Python 3.11
- pysam, cyvcf2, pydantic, jsonschema
- elpgx package
- Data files

**Build:**
```bash
docker build -t elifffarslan/elpgx-python:1.1.0 -f docker/elpgx-python/Dockerfile .
```

---

### 7.2 elpgx-gatk

**Dockerfile:** `docker/gatk/Dockerfile`

**Purpose:** GATK variant calling.

**Base:** `broadinstitute/gatk:4.5.0.0`

**Note:** Can use base image directly in production.

---

### 7.3 elpgx-hla

**Dockerfile:** `docker/hla/Dockerfile`

**Purpose:** HLA typing (placeholder for arcasHLA/OptiType).

**Status:** Placeholder - implement with arcasHLA for production.

---

## 8. Data Files

### Data Source & Provenance

All pharmacogenomics data in this pipeline is sourced from **official clinical databases**:

| Source | Description | URL |
|--------|-------------|-----|
| **CPIC** | Clinical Pharmacogenetics Implementation Consortium | https://cpicpgx.org/ |
| **PharmVar** | Pharmacogene Variation Consortium | https://www.pharmvar.org/ |
| **PharmGKB** | Pharmacogenomics Knowledge Base | https://www.pharmgkb.org/ |

**Current Data Version:** CPIC v1.24 (February 2026)
**Reference Build:** GRCh38
**License:** CC BY-SA 4.0

To update data, run:
```bash
python3 scripts/download_cpic_complete.py
python3 scripts/add_defining_variants.py
```

### 8.1 Allele Definitions

**Location:** `data/allele_definitions/`

**Format:** JSON per gene with GRCh38 coordinates

**Example (CYP2C19.json):**
```json
{
  "gene": "CYP2C19",
  "chromosome": "chr10",
  "reference_build": "GRCh38",
  "cpic_version": "v1.24",
  "reference_allele": "*1",
  "sequence_locations": [
    {"rsid": "rs4244285", "chrom": "chr10", "pos": 94781859, "ref": "G", "alt": "A"}
  ],
  "alleles": [
    {
      "allele": "*2",
      "function": "No function",
      "activity_value": 0.0,
      "evidence_strength": "Definitive",
      "defining_variants": [
        {"rsid": "rs4244285", "chrom": "chr10", "pos": 94781859, "ref": "G", "alt": "A"}
      ]
    }
  ]
}
```

**Available Genes (18 total):**
- Group 1: CYP2C19, CYP2C9, VKORC1, SLCO1B1, TPMT, DPYD, UGT1A1, NUDT15
- Group 2: CYP2D6, CYP3A5, CYP2B6, G6PD
- Group 3: NAT2, IFNL3, CACNA1S, RYR1
- Group 4: HLA-A, HLA-B

---

### 8.2 Phenotype Maps

**Location:** `data/phenotype_maps/`

**Purpose:** Map diplotypes to phenotypes.

**Example (CYP2C19_phenotype.json):**
```json
{
  "gene": "CYP2C19",
  "diplotype_phenotype_map": {
    "*1/*1": {
      "phenotype": "Normal Metabolizer",
      "activity_score": 2.0
    },
    "*1/*17": {
      "phenotype": "Rapid Metabolizer",
      "activity_score": 2.5
    }
  }
}
```

---

### 8.3 Prescribing Implications

**Location:** `data/prescribing_implications/`

**Purpose:** Clinical guidance based on phenotype.

**Example (CYP2C19_implications.json):**
```json
{
  "gene": "CYP2C19",
  "implications": [
    {
      "phenotype": "Poor Metabolizer",
      "implication": "Clopidogrel: Consider alternative antiplatelet therapy. PPIs: Standard dosing likely effective."
    }
  ]
}
```

---

### 8.4 PGx Regions BED

**Location:** `data/pgx_regions.bed`

**Purpose:** Target genomic regions for variant calling.

**Format:**
```
chr1    97070797    97921042    DPYD
chr10   94757546    94855547    CYP2C19
chr10   94938657    94990091    CYP2C9
...
```

---

### 8.5 Gene Tiers Configuration

**Location:** `data/gene_tiers.json`

**Purpose:** Gene classification and metadata.

**Format:**
```json
{
  "genes": {
    "CYP2C19": {
      "tier": "group1",
      "biomarker_type": "metabolic_enzyme",
      "chromosome": "chr10",
      "cpic_level": "A"
    }
  }
}
```

---

## 9. Configuration Files

### 9.1 nextflow.config

**Purpose:** Nextflow pipeline configuration.

**Key Sections:**

```groovy
// Parameters
params {
    input = null              // BAM/CRAM path
    reference = null          // Reference FASTA
    gene_tiers = 'group1'     // Tiers to analyze
    enable_hla = false        // HLA module
    outdir = './results'      // Output directory
}

// Process containers
process {
    withName: 'VARIANT_CALL' {
        container = 'broadinstitute/gatk:4.5.0.0'
        cpus = 4
        memory = 8.GB
    }
}

// Profiles
profiles {
    docker { docker.enabled = true }
    sevenbridges { process.executor = 'sbg' }
    awsbatch { process.executor = 'awsbatch' }
}
```

---

### 9.2 pyproject.toml

**Purpose:** Python package configuration.

**Contents:**
- Package metadata
- Dependencies
- Entry points (CLI commands)
- Build settings

---

### 9.3 docker-compose.yml

**Purpose:** Build all Docker images.

**Usage:**
```bash
docker-compose build        # Build all images
docker-compose push         # Push to registry
```

---

## 10. Usage Examples

### Example 1: Quick Analysis (CLI)

```bash
# Analyze with Group 1 genes only
python -m elpgx.cli run \
  -b sample.cram \
  -r hg38.fa \
  -g group1
```

### Example 2: Full Analysis (CLI)

```bash
# Analyze with Groups 1, 2, 3
python -m elpgx.cli run \
  -b sample.cram \
  -r hg38.fa \
  -g group1,group2,group3 \
  -o results/sample_report.json \
  --pretty
```

### Example 3: Production (Nextflow)

```bash
nextflow run main.nf \
  --input s3://bucket/sample.bam \
  --reference s3://bucket/hg38.fa \
  --gene_tiers group1,group2,group3 \
  --outdir s3://bucket/results \
  -profile docker,awsbatch
```

### Example 4: Seven Bridges Platform

```bash
nextflow run main.nf \
  --input /sbgenomics/project-files/sample.bam \
  --reference /sbgenomics/project-files/hg38.fa \
  -profile sevenbridges
```

---

## 11. Output Format

### JSON Report Structure

```json
{
  "meta": {
    "sample_id": "HG00100",
    "pipeline_version": "1.0.3",
    "run_date": "2026-02-13T15:30:00Z",
    "reference_build": "GRCh38",
    "pgx_db_version": "CPIC v1.22",
    "gene_tiers_included": ["group1", "group2", "group3"],
    "hla_module_enabled": false
  },
  "results": [
    {
      "gene": "CYP2C19",
      "gene_tier": "group1",
      "biomarker_type": "metabolic_enzyme",
      "diplotype": "*1/*17",
      "functional_status": "Normal Function / Increased Function",
      "phenotype": "Rapid Metabolizer",
      "activity_score": 2.5,
      "prescribing_implication": "Clopidogrel: Standard dosing. PPIs: May need higher dose.",
      "confidence": {
        "score": 0.95,
        "flags": []
      },
      "variants": [
        {
          "rsid": "rs12248560",
          "chrom": "chr10",
          "pos": 94762706,
          "ref": "C",
          "alt": "T",
          "genotype": "0/1",
          "quality": 38.5,
          "depth": 42
        }
      ]
    }
  ],
  "warnings": [],
  "provenance": {
    "pipeline_name": "elpgx-pipeline",
    "pipeline_version": "1.0.3"
  }
}
```

---

## 12. Troubleshooting

### Issue: "No module named 'elpgx'"

**Solution:** Install the package:
```bash
pip install -e .
```

### Issue: "GATK not found"

**Solution:** GATK is bundled in the standalone app. For CLI:
```bash
conda install -c bioconda gatk4
```

### Issue: "Reference does not appear to be hg38"

**Solution:** Use GRCh38 reference with chr1-chr22 naming:
```bash
# Download correct reference
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
```

### Issue: "Docker not found"

**Solution:** Install Docker Desktop:
```bash
brew install --cask docker
```

### Issue: "Java not found"

**Solution:** Install Java or use Docker:
```bash
brew install openjdk@17
```

---

## Summary

| Component | File | Purpose |
|-----------|------|---------|
| **CLI** | `elpgx/cli.py` | Command-line interface |
| **GUI** | `elpgx/gui.py` | Graphical interface |
| **Validation** | `elpgx/validate.py` | Input file validation |
| **Variant Calling** | `elpgx/vcf_parser.py` | GATK/bcftools integration |
| **Allele Calling** | `elpgx/allele_caller.py` | Star-allele matching |
| **Phenotype** | `elpgx/phenotype.py` | Clinical phenotype translation |
| **Confidence** | `elpgx/confidence.py` | Quality scoring |
| **Report** | `elpgx/report.py` | JSON report generation |
| **Nextflow** | `main.nf` | Pipeline orchestration |
| **Docker** | `docker/` | Container images |
| **Data** | `data/` | Allele definitions, phenotype maps |

---

*Document generated: February 2026*
*ElPGx Version: 1.0.3*
