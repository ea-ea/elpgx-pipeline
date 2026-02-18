# PGx Pipeline: Strategy & Feasibility Presentation

## Pharmacogenetics Integration for Clinical NGS Platform

**20-Minute Product Team Presentation**

---

# Agenda

1. Current State of Pharmacogenomics (5 min)
2. Architectural Defense (7 min)
3. Feasibility Analysis (5 min)
4. Product & UI/UX Vision (3 min)

---

# Part 1: Current State of the Art

---

## 1.1 Industry Approaches to Pharmacogenomics

### Leading Approaches

| Approach | Examples | Coverage |
|----------|----------|----------|
| **Targeted PGx Panels** | Genotyping arrays | 15-30 genes, SNPs only |
| **WES/WGS-derived PGx** | PharmCAT, Aldy | Broader coverage, variable quality |
| **Clinical Decision Support** | EHR-integrated alerts | Real-time prescribing guidance |

### Standard Reporting Elements
- Star-allele diplotypes (e.g., CYP2C19 *1/*17)
- Metabolizer phenotypes (Poor → Ultrarapid)
- CPIC/DPWG guideline-based interpretations
- Drug-specific dosing recommendations

---

## 1.2 PGx in Rare Disease Patients

### Current Adoption Status

| Status | Reality |
|--------|---------|
| Standard practice? | **No** - but rapidly emerging |
| Clinical value? | **High** - innovative and actionable |
| Implementation | Translational/precision medicine programs |

### Why Rare Disease Patients Need PGx

- **High medication burden** → More drug interactions
- **Polypharmacy** → Higher ADR risk
- **WES/WGS already available** → Reuse existing data
- **No additional testing required** → Cost-effective

> **Opportunity**: Convert diagnostic sequencing into therapeutic guidance

---

## 1.3 Sequencing Modalities Comparison

| Modality | Pros | Cons | PGx Suitability |
|----------|------|------|-----------------|
| **WES** | Lower cost, widely available | Uneven coverage, misses regulatory regions | Medium |
| **WGS** | Uniform coverage, CNV detection | Higher cost, more data | High |
| **Short-read** | Industry standard, validated | Complex loci difficult (CYP2D6) | High (most genes) |
| **Long-read** | Superior for structural variants | Limited availability, cost | Future for complex genes |

### Industry Trajectory
- **Near-term**: Short-read WGS as clinical standard
- **Future**: Long-read for complex PGx genes (CYP2D6)

---

# Part 2: Architectural Defense

---

## 2.1 Pipeline Architecture Overview

```
BAM/CRAM (hg38)
     │
     ▼
┌──────────────┐
│   VALIDATE   │  ← Verify hg38, check index, extract sample ID
└──────┬───────┘
       │
       ▼
┌──────────────┐
│   VARIANT    │  ← GATK HaplotypeCaller (PGx regions only)
│   CALLING    │
└──────┬───────┘
       │
       ▼
┌──────────────┐
│   ALLELE     │  ← Star-allele matching engine
│   CALLING    │
└──────┬───────┘
       │
       ▼
┌──────────────┐
│  PHENOTYPE   │  ← CPIC-based translation
│ TRANSLATION  │
└──────┬───────┘
       │
       ▼
┌──────────────┐
│   REPORT     │  ← Schema-validated JSON
│ GENERATION   │
└──────────────┘
```

---

## 2.2 Technology Stack Rationale

| Component | Choice | Why |
|-----------|--------|-----|
| **Orchestration** | Nextflow DSL2 | Industry standard, reproducible, scalable |
| **Variant Calling** | GATK HaplotypeCaller | Clinical-grade, well-validated |
| **Logic Layer** | Python 3.11+ | Type safety, Pydantic models, testability |
| **Containers** | Docker | Reproducibility, version control |
| **Data Models** | Pydantic v2 | Validation, serialization, type hints |
| **Reference Data** | CPIC/PharmVar | Regulatory-compliant, peer-reviewed |

### Key Design Principles
1. **Separation of concerns**: Heavy computation ≠ interpretation logic
2. **Deterministic outputs**: Pinned versions, no `:latest` tags
3. **Auditable**: Full provenance in output JSON
4. **Extensible**: Module contracts enable staged rollout

---

## 2.3 Algorithm Choices

### Allele Calling Strategy

**Mode A: VCF-based (Group 1 Genes)**
```
VCF → Match variants → Identify haplotypes → Resolve diplotype
```
- Works for: CYP2C19, CYP2C9, VKORC1, SLCO1B1, TPMT, DPYD, UGT1A1, NUDT15
- Confidence: High (SNP-based)

**Mode B: BAM-based (Future - Group 2)**
```
BAM → CNV detection → Hybrid detection → Complex diplotype
```
- Required for: CYP2D6
- Tools: Cyrius, Stargazer, Aldy

### Phenotype Translation
```
Diplotype → Activity Score → Phenotype → Prescribing Implication
     ↓            ↓              ↓               ↓
  *1/*17        2.5         Rapid         "Consider alternative"
```

---

## 2.4 Handling Complex PGx Genes

### CYP2D6 Challenges

| Challenge | Impact | Our Strategy |
|-----------|--------|--------------|
| **Gene deletions (*5)** | Missed calls | BAM-based detection (v1.1) |
| **Duplications (*1xN)** | Copy number | Dedicated CNV module |
| **CYP2D7 hybrids** | False positives | Specialized caller integration |
| **Pseudogene homology** | Mismapping | Long-read future pathway |

### Confidence Scoring Formula

```
Confidence = 0.40 × Depth + 0.35 × Quality + 0.25 × Completeness
```

### Quality Flags
- `LOW_DEPTH` - Any position < 10x
- `LOW_QUALITY` - QUAL < 20
- `INCOMPLETE_ALLELE` - Missing defining variants
- `AMBIGUOUS_CALL` - Multiple equally likely diplotypes

---

# Part 3: Feasibility Analysis

---

## 3.1 Technical Limitations

### NGS Limitations for PGx

| Limitation | Affected Genes | Mitigation |
|------------|----------------|------------|
| WES coverage gaps | All | Use WGS, flag low-coverage |
| Short-read CNV detection | CYP2D6 | Dedicated tools, defer to v1.1 |
| Homologous region mapping | CYP2D6, HLA | Specialized callers |
| Regulatory region variants | CYP2C19 *17 | WGS preferred |

### Group 1 Genes: High Feasibility
- All 8 genes are SNP-based
- Well-defined star alleles
- CPIC Level A evidence
- Standard variant calling sufficient

### Group 2+ Genes: Requires Specialized Handling
- CYP2D6: CNV calling module
- HLA: Dedicated HLA typers

---

## 3.2 Clinical Validation Strategy

### Phase 1: Pre-Release Validation

| Validation Type | Method | Target |
|-----------------|--------|--------|
| **Reference samples** | GeT-RM, Coriell | Known diplotypes |
| **Orthogonal comparison** | vs. targeted PGx panels | Concordance rate |
| **CPIC concordance** | Expected phenotypes | 100% for Group 1 |

### Phase 2: Gradual Rollout

```
Week 1-4:   Group 1 genes only (8 genes)
Week 5-8:   Add CYP2D6 with confidence flags
Week 9+:    HLA module (opt-in)
```

### Documentation Requirements
- Performance characteristics per gene
- Known limitations and edge cases
- Confidence thresholds for clinical use
- Recommendations for confirmatory testing

---

## 3.3 Accuracy Expectations

### Group 1 Genes (v1.0)

| Gene | Expected Concordance | Confidence |
|------|---------------------|------------|
| CYP2C19 | >99% | High |
| CYP2C9 | >99% | High |
| VKORC1 | >99% | High |
| SLCO1B1 | >99% | High |
| TPMT | >98% | High |
| DPYD | >95% | Medium (rare variants) |
| UGT1A1 | >95% | Medium (*28 is TA-repeat) |
| NUDT15 | >98% | High |

---

# Part 4: Product & UI/UX Vision

---

## 4.1 Clinical Use Case

### Scenario: Clopidogrel Prescription

```
┌─────────────────────────────────────────────────────────┐
│  PRESCRIBING: Clopidogrel 75mg                          │
├─────────────────────────────────────────────────────────┤
│  ⚠️  PGx ALERT: CYP2C19 *1/*2 (Intermediate Metabolizer)│
│                                                          │
│  Recommendation: Consider prasugrel or ticagrelor        │
│  Guideline: CPIC 2022                                    │
│                                                          │
│  [View Details]  [Override with Reason]  [Accept]        │
└─────────────────────────────────────────────────────────┘
```

### Integration Points
- EHR prescribing workflow
- Clinical decision support alerts
- Patient summary dashboard

---

## 4.2 Communicating Uncertainty

### Design Principles

1. **Don't overwhelm** - Hide raw genomic data by default
2. **Surface uncertainty clearly** - Confidence indicators
3. **Enable drill-down** - Expandable evidence sections

### UI Confidence Indicators

| Confidence | Display | Meaning |
|------------|---------|---------|
| High (>0.9) | ✅ Green | Report as-is |
| Moderate (0.7-0.9) | 🟡 Yellow | Clinical review suggested |
| Low (<0.7) | 🔴 Red | Confirmatory testing recommended |

### Example Report Card

```
┌─────────────────────────────────────────────┐
│  CYP2C19                                    │
│  ─────────────────────────────────────────  │
│  Diplotype:    *1/*17                       │
│  Phenotype:    Rapid Metabolizer     ✅      │
│  Confidence:   95%                          │
│                                             │
│  ▸ View variants (1 found)                  │
│  ▸ View guideline (CPIC 2022)               │
│  ▸ Drug implications                        │
└─────────────────────────────────────────────┘
```

---

## 4.3 Report Layering

### Audience-Appropriate Views

| Audience | View | Content |
|----------|------|---------|
| **Physician** | Summary | Phenotype + drug guidance |
| **Pharmacist** | Detailed | Activity scores + dosing |
| **Lab Director** | Technical | Variants + QC metrics |
| **EHR System** | Machine-readable | Full JSON |

---

# Summary

## What We've Demonstrated

| Area | Evidence |
|------|----------|
| **Industry knowledge** | Current PGx landscape understanding |
| **Technical strategy** | Defensible architecture with clear rationale |
| **Feasibility awareness** | Honest assessment of limitations |
| **Product vision** | Clinician-centered UI/UX concepts |

## Key Differentiators

1. **Reuses diagnostic WGS** - No additional testing
2. **CPIC-compliant** - Regulatory defensibility
3. **Confidence scoring** - Transparent uncertainty
4. **Modular design** - Staged rollout capability
5. **Clinical-grade** - Schema validation, provenance

---

# Questions?

## Appendix Available
- Detailed technical rationale document
- Full pipeline specification
- Test strategy documentation
- Data model schemas
