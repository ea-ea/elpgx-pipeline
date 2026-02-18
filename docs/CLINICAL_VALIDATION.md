# ElPGx Clinical-Grade Validation System

**Version:** 1.1.0
**Last Updated:** 2026-02-14
**Compliance Target:** FDA 21 CFR Part 11, CAP/CLIA Laboratory Standards
**Author:** Elif Arslan

---

## Table of Contents

1. [Overview](#overview)
2. [Validation Categories](#validation-categories)
3. [Severity Levels](#severity-levels)
4. [Critical Gene Requirements](#critical-gene-requirements)
5. [Validation Checks Reference](#validation-checks-reference)
6. [Running Validation](#running-validation)
7. [Validation Report Format](#validation-report-format)
8. [Remediation Guide](#remediation-guide)
9. [Regulatory Compliance](#regulatory-compliance)
10. [Audit Trail](#audit-trail)

---

## Overview

The ElPGx Clinical Validation System ensures that the pharmacogenomics pipeline meets the rigorous standards required for clinical laboratory use. This system performs comprehensive checks before any analysis to verify data integrity, completeness, and safety rules.

### Purpose

- **Patient Safety**: Ensure accurate phenotype predictions for drug dosing
- **Data Integrity**: Verify reference data has not been corrupted or tampered with
- **Reproducibility**: Track all versions and checksums for audit trails
- **Compliance**: Meet FDA, CAP, and CLIA requirements for clinical software

### When Validation Runs

| Entry Point | Validation Mode | Behavior on Failure |
|-------------|-----------------|---------------------|
| CLI startup | Non-strict | Warns on errors, blocks on critical |
| GUI launch | Non-strict | Shows dialog, blocks on critical |
| Nextflow pipeline | Strict | Fails workflow on any error |
| Manual run | Configurable | Returns detailed report |

---

## Validation Categories

### CV-100: Data Integrity

**Purpose:** Verify that all reference data files are intact and from known sources.

| Code | Check | Severity |
|------|-------|----------|
| CV-101 | Calculate checksums for allele definitions | INFO |
| CV-102 | Allele definitions directory exists | CRITICAL |
| CV-103 | Phenotype maps directory exists | CRITICAL |
| CV-104 | CPIC version matches expected | WARNING |
| CV-105 | CPIC version verified | INFO |
| CV-106 | PROVENANCE.json valid JSON | CRITICAL |
| CV-107 | PROVENANCE.json exists | ERROR |

**Checksums Calculated:**
- Combined SHA256 of all allele definition files
- Combined SHA256 of all phenotype map files
- Individual checksum of PROVENANCE.json
- Individual checksum of schema file

---

### CV-200: Critical Gene Completeness

**Purpose:** Ensure genes with high clinical impact have complete allele definitions.

**Why This Matters:**
Missing allele definitions can lead to incorrect phenotype predictions, potentially causing:
- Drug toxicity (DPYD → 5-FU)
- Therapeutic failure (CYP2C19 → Clopidogrel)
- Severe adverse reactions (TPMT → Azathioprine)

| Code | Check | Severity |
|------|-------|----------|
| CV-201 | Critical gene file exists | CRITICAL |
| CV-202 | Minimum allele count met | ERROR |
| CV-203 | Required alleles present | ERROR |
| CV-204 | All required alleles verified | INFO |
| CV-205 | Gene definition valid JSON | CRITICAL |

---

### CV-300: Phenotype Map Completeness

**Purpose:** Verify that diplotype-to-phenotype translations cover all expected categories.

| Code | Check | Severity |
|------|-------|----------|
| CV-301 | Phenotype map file exists | ERROR |
| CV-302 | All phenotype categories covered | WARNING |
| CV-303 | Phenotype coverage verified | INFO |
| CV-304 | Phenotype map valid JSON | ERROR |

**Required Phenotype Categories:**

| Gene | Required Phenotypes |
|------|---------------------|
| CYP2C19 | Ultrarapid, Normal, Intermediate, Poor Metabolizer |
| CYP2C9 | Normal, Intermediate, Poor Metabolizer |
| CYP2D6 | Ultrarapid, Normal, Intermediate, Poor Metabolizer |
| DPYD | Normal, Intermediate, Poor Metabolizer |
| TPMT | Normal, Intermediate, Poor Metabolizer |

---

### CV-400: Clinical Safety Rules

**Purpose:** Verify that safety mechanisms are properly configured.

| Code | Check | Severity |
|------|-------|----------|
| CV-401 | CYP2D6 no_call blocks phenotype | CRITICAL |
| CV-402 | Phenotype blocking verified | INFO |
| CV-403 | CNV limitation warning configured | WARNING |
| CV-404 | Strategy file valid JSON | ERROR |
| CV-405 | Gene strategies directory exists | WARNING |

**Safety Rules Verified:**

1. **Phenotype Blocking**: Low-confidence calls must NOT report phenotypes
2. **CNV Warnings**: CYP2D6 must warn about CNV detection limitations
3. **Confidence Thresholds**: Minimum thresholds must be enforced

---

### CV-500: Schema Compliance

**Purpose:** Ensure output schema meets clinical reporting requirements.

| Code | Check | Severity |
|------|-------|----------|
| CV-501 | Schema file exists | CRITICAL |
| CV-502 | Required top-level fields | ERROR |
| CV-503 | Required meta fields | ERROR |
| CV-504 | Provenance field exists | ERROR |
| CV-505 | Provenance verified | INFO |
| CV-506 | Confidence required in results | ERROR |
| CV-507 | Schema valid JSON | CRITICAL |

**Required Schema Fields:**

```
Top-level: schema_version, meta, results
Meta: sample_id, pipeline_version, run_date, reference_build, pgx_db_version
Results: gene, diplotype, phenotype, confidence
```

---

### CV-600: Version Tracking

**Purpose:** Ensure data provenance through version files.

| Code | Check | Severity |
|------|-------|----------|
| CV-601 | VERSION file present | INFO |
| CV-602 | VERSION file missing | ERROR |

**Required VERSION Files:**

```
data/VERSION
data/allele_definitions/VERSION
data/phenotype_maps/VERSION
data/prescribing_implications/VERSION
schemas/VERSION
```

---

## Severity Levels

| Level | Symbol | Description | Pipeline Behavior |
|-------|--------|-------------|-------------------|
| **CRITICAL** | 🔴 | Blocks all clinical use | Immediate halt |
| **ERROR** | 🟠 | Must resolve for production | Halt in strict mode |
| **WARNING** | 🟡 | Should be reviewed | Log and continue |
| **INFO** | 🔵 | Informational only | Log only |

### Strict vs Non-Strict Mode

| Mode | CRITICAL | ERROR | WARNING | INFO |
|------|----------|-------|---------|------|
| Strict | FAIL | FAIL | PASS | PASS |
| Non-Strict | FAIL | PASS | PASS | PASS |

**Recommendation:** Use strict mode for all clinical production runs.

---

## Critical Gene Requirements

### Gene: DPYD

| Requirement | Value | Clinical Impact |
|-------------|-------|-----------------|
| Minimum Alleles | 5 | 5-FU/Capecitabine toxicity |
| Required Alleles | *1, *2A, *13 | Fatal toxicity risk |
| Safety Level | CRITICAL | |

**Clinical Context:** DPYD deficiency can cause fatal fluoropyrimidine toxicity. Missing the *2A allele (most common variant) could result in a patient receiving full-dose chemotherapy when dose reduction is required.

### Gene: TPMT

| Requirement | Value | Clinical Impact |
|-------------|-------|-----------------|
| Minimum Alleles | 4 | Thiopurine myelosuppression |
| Required Alleles | *1, *2, *3A, *3C | Life-threatening toxicity |
| Safety Level | CRITICAL | |

**Clinical Context:** TPMT deficiency causes severe, potentially fatal myelosuppression with standard azathioprine/6-MP doses.

### Gene: CYP2C19

| Requirement | Value | Clinical Impact |
|-------------|-------|-----------------|
| Minimum Alleles | 8 | Clopidogrel efficacy |
| Required Alleles | *1, *2, *3, *17 | Cardiovascular events |
| Safety Level | CRITICAL | |

**Clinical Context:** Poor metabolizers have increased risk of cardiovascular events on clopidogrel. Ultrarapid metabolizers may have increased bleeding risk.

### Gene: CYP2C9

| Requirement | Value | Clinical Impact |
|-------------|-------|-----------------|
| Minimum Alleles | 6 | Warfarin bleeding risk |
| Required Alleles | *1, *2, *3 | Hemorrhage risk |
| Safety Level | CRITICAL | |

**Clinical Context:** CYP2C9 variants significantly affect warfarin dose requirements. Incorrect dosing can cause life-threatening bleeding.

### Gene: NUDT15

| Requirement | Value | Clinical Impact |
|-------------|-------|-----------------|
| Minimum Alleles | 3 | Thiopurine toxicity (Asian) |
| Required Alleles | *1, *3 | Myelosuppression |
| Safety Level | CRITICAL | |

**Clinical Context:** NUDT15 variants are particularly important in Asian populations for thiopurine dosing.

### Gene: UGT1A1

| Requirement | Value | Clinical Impact |
|-------------|-------|-----------------|
| Minimum Alleles | 3 | Irinotecan neutropenia |
| Required Alleles | *1, *28, *80 | Severe neutropenia |
| Safety Level | CRITICAL | |

**Clinical Context:** UGT1A1*28 homozygotes have significantly increased risk of severe neutropenia with irinotecan.

---

## Running Validation

### Command Line

```bash
# Run validation standalone
python -m elpgx.clinical_validation

# Run with specific project root
python -c "
from pathlib import Path
from elpgx.clinical_validation import run_clinical_validation, print_validation_report

report = run_clinical_validation(Path('.'), strict=True, raise_on_failure=False)
print_validation_report(report)
"
```

### Programmatic Usage

```python
from pathlib import Path
from elpgx.clinical_validation import (
    run_clinical_validation,
    ClinicalValidationError,
    ValidationSeverity,
)

# Run validation
report = run_clinical_validation(
    project_root=Path("/path/to/elpgx"),
    strict=True,           # Treat errors as failures
    raise_on_failure=True  # Raise exception on failure
)

# Check results
if report.passed:
    print("Validation PASSED")
else:
    print("Validation FAILED")
    for finding in report.findings:
        if finding.severity in [ValidationSeverity.CRITICAL, ValidationSeverity.ERROR]:
            print(f"[{finding.code}] {finding.message}")
```

### In Nextflow Pipeline

The validation runs automatically at pipeline start:

```groovy
workflow {
    // Validation runs here before any processing
    VALIDATE_CLINICAL()

    // Pipeline continues only if validation passes
    VARIANT_CALL(validated_input)
    // ...
}
```

---

## Validation Report Format

### JSON Schema

```json
{
  "validation_id": "CV-20260214153042-a1b2c3d4",
  "timestamp": "2026-02-14T15:30:42.123456+00:00",
  "pipeline_version": "1.1.0",
  "validation_mode": "full",
  "passed": true,
  "findings": [
    {
      "code": "CV-101",
      "severity": "INFO",
      "category": "Data Integrity",
      "message": "Calculated checksums for 8 allele definition files",
      "details": "Combined hash: 969982cf1562accf...",
      "remediation": null
    }
  ],
  "checksums": {
    "allele_definitions": "sha256:abcd1234...",
    "phenotype_maps": "sha256:efgh5678...",
    "schema": "sha256:ijkl9012..."
  },
  "metadata": {
    "strict_mode": true,
    "critical_count": 0,
    "error_count": 0,
    "warning_count": 2,
    "info_count": 15
  }
}
```

### Console Output

```
======================================================================
ElPGx CLINICAL VALIDATION REPORT
======================================================================
Validation ID: CV-20260214153042-a1b2c3d4
Timestamp:     2026-02-14T15:30:42.123456+00:00
Version:       1.1.0
Mode:          full
Status:        PASSED
----------------------------------------------------------------------

INFO (15):
  [CV-101] Calculated checksums for 8 allele definition files
  [CV-105] CPIC version verified: v1.24
  [CV-204] DPYD: All 3 required alleles present
  ...

WARNING (2):
  [CV-302] NAT2: Missing phenotype mappings: Rapid Acetylator
           Fix: Add missing diplotype-to-phenotype mappings

----------------------------------------------------------------------
Summary: 0 critical, 0 errors, 2 warnings
======================================================================
```

---

## Remediation Guide

### CRITICAL Issues

#### CV-102: Allele definitions directory not found

```bash
# Verify data directory structure
ls -la data/allele_definitions/

# If missing, restore from backup or re-download
python scripts/download_cpic_data.py
```

#### CV-201: Critical gene definition file missing

```bash
# Check which gene is missing
ls data/allele_definitions/*.json

# Download specific gene from CPIC
curl -o data/allele_definitions/DPYD.json \
  "https://api.cpicpgx.org/v1/gene/DPYD"
```

#### CV-401: CYP2D6 phenotype blocking not configured

```bash
# Edit strategy file
vim data/gene_strategies/CYP2D6_STRATEGY.json

# Ensure this setting exists:
{
  "confidence_thresholds": {
    "no_call": {
      "phenotype_reported": false,
      "action": "block_phenotype"
    }
  }
}
```

### ERROR Issues

#### CV-202: Insufficient alleles for critical gene

1. Check current allele count in gene definition
2. Compare against CPIC database
3. Add missing allele definitions
4. Update VERSION file

#### CV-602: Missing VERSION file

```bash
# Create VERSION file
echo "1.1.0" > data/allele_definitions/VERSION
echo "cpic_version: v1.24" >> data/allele_definitions/VERSION
echo "last_updated: $(date +%Y-%m-%d)" >> data/allele_definitions/VERSION
```

### WARNING Issues

#### CV-104: CPIC version mismatch

1. Verify data source version
2. Update PROVENANCE.json if intentional
3. Or update data to match expected version

#### CV-302: Missing phenotype mappings

1. Review phenotype map file
2. Add missing diplotype entries
3. Consult CPIC guidelines for mappings

---

## Regulatory Compliance

### FDA 21 CFR Part 11

| Requirement | Implementation |
|-------------|----------------|
| Electronic signatures | Validation ID + timestamp |
| Audit trails | Full finding history in report |
| Data integrity | SHA256 checksums |
| Version control | VERSION files in all directories |

### CAP/CLIA Laboratory Standards

| Standard | Implementation |
|----------|----------------|
| Validation documentation | Validation report JSON |
| Quality control | Pre-run checks |
| Proficiency testing | Test suite with known samples |
| Error detection | Multi-level severity system |

### ISO 15189 (Medical Laboratories)

| Clause | Implementation |
|--------|----------------|
| 5.5 Examination processes | Validation before analysis |
| 5.6 Quality assurance | Checksum verification |
| 5.9 Release of results | Block on validation failure |

---

## Audit Trail

### What Is Logged

Every validation run creates an audit record containing:

1. **Validation ID**: Unique identifier (CV-YYYYMMDDHHMMSS-RANDOM)
2. **Timestamp**: ISO8601 UTC timestamp
3. **Pipeline Version**: Software version that ran validation
4. **All Findings**: Complete list with severity and remediation
5. **Checksums**: Cryptographic hashes of all data files
6. **Pass/Fail Status**: Final determination

### Retention Requirements

| Environment | Retention Period | Storage |
|-------------|------------------|---------|
| Development | 30 days | Local filesystem |
| Testing | 90 days | Version control |
| Production | 7 years | Secure archive |

### Accessing Audit Records

```python
# Save validation report for audit
import json
from datetime import datetime

report = run_clinical_validation(project_root)

audit_file = f"audits/validation_{report.validation_id}.json"
with open(audit_file, "w") as f:
    json.dump(report.model_dump(), f, indent=2, default=str)
```

---

## Appendix A: Complete Validation Code Reference

### Finding Codes

| Range | Category |
|-------|----------|
| CV-100 - CV-199 | Data Integrity |
| CV-200 - CV-299 | Critical Gene Completeness |
| CV-300 - CV-399 | Phenotype Map Completeness |
| CV-400 - CV-499 | Clinical Safety Rules |
| CV-500 - CV-599 | Schema Compliance |
| CV-600 - CV-699 | Version Tracking |

### Expected CPIC Version

```python
EXPECTED_CPIC_VERSION = "v1.24"
```

### Confidence Thresholds

```python
CLINICAL_CONFIDENCE_THRESHOLDS = {
    "reportable_minimum": 0.70,      # Below this = block phenotype
    "high_confidence": 0.90,         # High confidence threshold
    "cnv_gene_warning": 0.85,        # CYP2D6 without CNV warning
}
```

---

## Appendix B: Test Suite

The clinical validation system includes comprehensive tests:

```bash
# Run all clinical validation tests
pytest tests/unit/test_clinical_validation.py -v

# Run specific test category
pytest tests/unit/test_clinical_validation.py::TestCriticalGeneCompleteness -v

# Run with coverage
pytest tests/unit/test_clinical_validation.py --cov=elpgx.clinical_validation
```

### Test Categories

- `TestDataIntegrity`: 6 tests
- `TestCriticalGeneCompleteness`: 18 tests (3 per critical gene)
- `TestPhenotypeMapCompleteness`: 5 tests
- `TestClinicalSafetyRules`: 4 tests
- `TestSchemaCompliance`: 4 tests
- `TestVersionTracking`: 5 tests
- `TestFullValidation`: 3 tests
- `TestHashFunctions`: 3 tests

**Total: 48 tests**

---

## Document History

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0.0 | 2026-02-14 | Elif Arslan | Initial clinical validation system |
| 1.1.0 | 2026-02-14 | Elif Arslan | Added GUI integration, expanded tests |

---

*This document is part of the ElPGx Clinical Pharmacogenetics Pipeline quality management system.*
