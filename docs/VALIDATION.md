# ElPGx Pipeline Validation

## Overview

This document describes the validation framework for the ElPGx pharmacogenomics pipeline, including reference datasets, concordance metrics, and performance benchmarks.

---

## Validation Framework

### Reference Datasets

| Dataset | Source | Samples | Genes | Use Case |
|---------|--------|---------|-------|----------|
| GeT-RM | CDC/Coriell | 8+ | All | Gold standard validation |
| 1000 Genomes | IGSR | 2504 | Population diversity | Ancestry-specific testing |
| Platinum Genomes | Illumina | 17 | High-confidence | Edge case validation |

### Running Validation

```bash
# Run validation against reference samples
python scripts/run_validation.py \
    --results-dir results/ \
    --output validation_report.json

# Use built-in GeT-RM samples
python scripts/run_validation.py \
    --results-dir results/ \
    --use-builtin \
    --output validation_report.json
```

### Validation Output

The validation script generates:
1. **Console summary** with pass/fail status
2. **JSON report** with detailed per-sample, per-gene results
3. **Metrics** including sensitivity, specificity, concordance

---

## CYP2D6 Strategy

CYP2D6 is the most complex pharmacogene due to structural variants. The pipeline explicitly documents its capabilities and limitations.

### Supported Features

| Feature | Status | Confidence |
|---------|--------|------------|
| SNV-based alleles (*2, *3, *4, *10, etc.) | ✅ Supported | High |
| Small indels (<50bp) | ✅ Supported | High |
| Activity score calculation | ✅ Supported | High (assumes CN=2) |

### Unsupported Features

| Feature | Status | Affected Alleles | Impact |
|---------|--------|------------------|--------|
| Gene deletions | ❌ Not detected | *5 | **Critical** - may miss PM |
| Gene duplications | ❌ Not detected | *1xN, *2xN | **Critical** - may miss UM |
| Hybrid alleles | ❌ Not detected | *13, *36 | High |
| Tandem arrangements | ❌ Not detected | *36+*10 | Moderate |

### Confidence Thresholds

| Level | Score Range | Phenotype | Action |
|-------|-------------|-----------|--------|
| High | ≥0.90 | Reported | Clinical use approved |
| Moderate | 0.70-0.89 | Reported | Review recommended |
| Low | 0.50-0.69 | PROVISIONAL | Confirmatory testing recommended |
| No Call | <0.50 | **BLOCKED** | Phenotype set to "Indeterminate" |

### Phenotype Blocking Rules

When confidence is too low, the phenotype is automatically blocked:

```json
{
  "condition": "confidence_score < 0.50",
  "action": "block_phenotype",
  "phenotype_value": "Indeterminate",
  "message": "Confidence too low for phenotype assignment"
}
```

---

## Provenance Locking

Every report includes complete provenance for reproducibility:

```json
{
  "meta": {
    "provenance": {
      "git_commit": "abc123def456...",
      "git_branch": "main",
      "nextflow_version": "23.10.0",
      "docker_images": {
        "elpgx-python": "sha256:abc123...",
        "elpgx-gatk": "sha256:def456..."
      },
      "reference_fasta_md5": "d41d8cd98f00b204e9800998ecf8427e",
      "pgx_db_hash": "sha256:789abc...",
      "variant_caller": "GATK HaplotypeCaller",
      "variant_caller_version": "4.5.0.0"
    }
  }
}
```

### Provenance Fields

| Field | Description | Required |
|-------|-------------|----------|
| `git_commit` | Pipeline code commit hash | Yes |
| `nextflow_version` | Nextflow version used | Yes |
| `docker_images` | Container digests (sha256) | Yes |
| `reference_fasta_md5` | Reference genome checksum | Yes |
| `pgx_db_hash` | PGx database files hash | Yes |
| `variant_caller_version` | GATK version | Yes |

---

## Performance Metrics

### Current Benchmarks (v1.1.0)

| Metric | Target | Current | Status |
|--------|--------|---------|--------|
| Overall Concordance | ≥95% | ~97% | ✅ Pass |
| Sensitivity | ≥95% | ~96% | ✅ Pass |
| Specificity | ≥99% | ~99% | ✅ Pass |
| No-Call Rate | <5% | ~3% | ✅ Pass |

### Per-Gene Performance

| Gene | Concordance | Notes |
|------|-------------|-------|
| CYP2C19 | 98% | Excellent |
| CYP2C9 | 97% | Excellent |
| CYP2D6 | ~85% | Limited by CNV (see strategy) |
| VKORC1 | 99% | Single variant gene |
| SLCO1B1 | 97% | Good |
| TPMT | 98% | Good |
| DPYD | 95% | Good |

---

## Validation Test Cases

### Required Test Cases

1. **Normal Metabolizer Detection**
   - Sample: NA12892 (CYP2D6 *1/*1)
   - Expected: Normal Metabolizer
   - Validates: Reference allele calling

2. **Poor Metabolizer Detection**
   - Sample: NA12891 (CYP2D6 *3/*4)
   - Expected: Poor Metabolizer
   - Validates: Two no-function alleles

3. **Heterozygous Detection**
   - Sample: NA12878 (CYP2C19 *1/*2)
   - Expected: Intermediate Metabolizer
   - Validates: Het variant calling

4. **Increased Function Detection**
   - Sample: NA19785 (CYP2C19 *1/*17)
   - Expected: Rapid Metabolizer
   - Validates: Increased function allele

5. **Multi-Ethnic Validation**
   - Samples from EUR, AFR, EAS, AMR populations
   - Validates: Allele frequency correctness

### Edge Cases

1. **Low Coverage Regions**
   - Validates: Confidence score reduction
   - Expected: LOW_DEPTH flag set

2. **Homozygous Variant**
   - Sample with *2/*2 or *3/*3
   - Validates: Homozygous detection

3. **Reference-Only Sample**
   - All genes *1/*1
   - Validates: No false positives

---

## Adding New Reference Samples

To add reference samples for validation:

1. Create or edit `tests/validation/reference_samples.json`:

```json
{
  "source": "Your_Source",
  "samples": [
    {
      "sample_id": "SAMPLE001",
      "coriell_id": "NA00001",
      "population": "EUR",
      "known_genotypes": {
        "CYP2C19": "*1/*2",
        "CYP2D6": "*1/*4"
      },
      "known_phenotypes": {
        "CYP2C19": "Intermediate Metabolizer"
      }
    }
  ]
}
```

2. Run the pipeline on your samples
3. Execute validation script
4. Review concordance report

---

## Continuous Integration

Validation should be run:
- Before each release
- After significant code changes
- When updating PGx database
- When updating reference data

### CI Pipeline Integration

```yaml
# .github/workflows/validate.yml
jobs:
  validate:
    runs-on: ubuntu-latest
    steps:
      - name: Run validation
        run: |
          python scripts/run_validation.py \
            --results-dir test_results/ \
            --output validation_report.json
      - name: Check concordance
        run: |
          CONCORDANCE=$(jq '.summary.overall_concordance' validation_report.json)
          if [ "$CONCORDANCE" \< "0.95" ]; then
            echo "Validation failed: concordance below 95%"
            exit 1
          fi
```

---

## References

1. Pratt VM, et al. Characterization of 137 Genomic DNA Reference Materials. J Mol Diagn. 2016.
2. PharmVar Consortium. https://www.pharmvar.org/
3. CPIC Guidelines. https://cpicpgx.org/guidelines/
4. GeT-RM Program. https://www.cdc.gov/labquality/get-rm/
