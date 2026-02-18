# Sample Test Data

End-to-end test data for the ElPGx pipeline.

## Test Datasets

### 1. HG00100 (CRAM, GRCh38)
1000 Genomes Project exome sample aligned to GRCh38.

| File | Size | Source |
|------|------|--------|
| `HG00100.alt_bwamem_GRCh38DH.20150826.GBR.exome.cram` | ~5.8 GB | [1000 Genomes FTP](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/GBR/HG00100/exome_alignment/) |
| `HG00100.alt_bwamem_GRCh38DH.20150826.GBR.exome.cram.crai` | ~412 KB | Same location |
| `GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta` | ~2.9 GB | [NCBI](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/) |
| `GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.fai` | ~8 KB | Same location |

### 2. HG002 (BAM, GRCh37/hg19)
Genome in a Bottle (GIAB) reference sample aligned to GRCh37.
Tests the pipeline's GRCh37 coordinate translation capability.

| File | Size | Source |
|------|------|--------|
| `151002_*.posiSrt.markDup.bam` | ~10 GB | GIAB / Diagnostic sequencing |
| `151002_*.posiSrt.markDup.bai` | ~7 MB | Same location |
| `human_g1k_v37.fasta` | ~2.9 GB | [1000 Genomes FTP](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz) |
| `human_g1k_v37.fasta.fai` | ~3 KB | Generated with `samtools faidx` |
| `human_g1k_v37.dict` | ~14 KB | Generated with `gatk CreateSequenceDictionary` |

## Setup

```bash
# --- GRCh38 (HG00100 CRAM) ---
# 1. Download reference
curl -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta
samtools faidx GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta
gatk CreateSequenceDictionary -R GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta

# 2. Download CRAM
curl -O https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/GBR/HG00100/exome_alignment/HG00100.alt_bwamem_GRCh38DH.20150826.GBR.exome.cram
curl -O https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/GBR/HG00100/exome_alignment/HG00100.alt_bwamem_GRCh38DH.20150826.GBR.exome.cram.crai

# --- GRCh37 (HG002 BAM) ---
# 1. Download GRCh37 reference
curl -O https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
gunzip human_g1k_v37.fasta.gz
samtools faidx human_g1k_v37.fasta
gatk CreateSequenceDictionary -R human_g1k_v37.fasta
```

## Running E2E Tests

```bash
# Run only e2e tests
python -m pytest tests/e2e/ -v

# Run CRAM tests only (GRCh38)
python -m pytest tests/e2e/test_full_pipeline.py -v

# Run BAM tests only (GRCh37)
python -m pytest tests/e2e/test_bam_pipeline.py -v

# Run all tests (unit + e2e)
python -m pytest tests/ -v
```

Tests are automatically skipped when sample data is not present.

## Expected Results

### HG00100 (GRCh38 CRAM)

| Gene | Diplotype | Phenotype | Confidence |
|------|-----------|-----------|------------|
| CYP2C19 | *1/*1 | Normal Metabolizer | 90% |
| CYP2D6 | *13+*2/*13+*2 | Normal Metabolizer | 70% |
| CYP2B6 | *1/*6 | Intermediate Metabolizer | 95% |
| NAT2 | *6/*6 | Poor Metabolizer | 70% |
| SLCO1B1 | *1/*1 | Normal Function | 90% |
| All others | *1/*1 | Normal Metabolizer | 90% |

Population: GBR (British, European ancestry)

### HG002 (GRCh37 BAM)

| Gene | Diplotype | Phenotype | Confidence |
|------|-----------|-----------|------------|
| NAT2 | *6/*60 | Indeterminate | 70% |
| All others | *1/*1 | Normal Metabolizer | 90% |

Population: Ashkenazi Jewish ancestry (GIAB reference)
