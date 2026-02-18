# PGx Pipeline: Technical Rationale Document

## Comprehensive Architecture and Technology Decisions

This document provides detailed justification for every major technical decision in the PGx Pipeline, including algorithm selection, language choices, package dependencies, and architectural patterns.

---

## Table of Contents

1. [Language Selection](#1-language-selection)
2. [Workflow Orchestration](#2-workflow-orchestration)
3. [Variant Calling Tools](#3-variant-calling-tools)
4. [Python Packages](#4-python-packages)
5. [Allele Calling Algorithm](#5-allele-calling-algorithm)
6. [Phenotype Translation](#6-phenotype-translation)
7. [Confidence Scoring](#7-confidence-scoring)
8. [Data Architecture](#8-data-architecture)
9. [Containerization Strategy](#9-containerization-strategy)
10. [Testing Framework](#10-testing-framework)

---

## 1. Language Selection

### Python 3.11+ for Logic Layer

**Decision**: Use Python 3.11+ for all interpretation logic (validation, allele calling, phenotype translation, reporting).

**Rationale**:

| Factor | Python | Alternatives Considered |
|--------|--------|------------------------|
| **Type Safety** | Native type hints, Pydantic | R (limited), Perl (none) |
| **Bioinformatics Ecosystem** | pysam, cyvcf2, BioPython | R/Bioconductor (different paradigm) |
| **Clinical Software Standards** | Widely accepted in clinical bioinformatics | |
| **Testability** | pytest, excellent tooling | |
| **Maintainability** | Readable, large talent pool | |

**Why Python 3.11+?**
- `match/case` statements for cleaner allele matching
- Native `|` union types (e.g., `str | None`)
- Performance improvements (10-60% faster than 3.10)
- Long-term support (until 2027)

**Rejected Alternatives**:
- **R**: Excellent for statistical analysis but less suitable for production clinical software
- **Java**: Verbosity overhead for relatively simple logic
- **Go**: Limited bioinformatics library ecosystem
- **Rust**: Learning curve, smaller bioinformatics community

---

## 2. Workflow Orchestration

### Nextflow DSL2

**Decision**: Use Nextflow with DSL2 syntax for pipeline orchestration.

**Rationale**:

| Criterion | Nextflow | Snakemake | CWL | WDL |
|-----------|----------|-----------|-----|-----|
| **Industry Adoption** | High (clinical labs) | High (academic) | Medium | Medium |
| **Docker Support** | Native, seamless | Good | Good | Good |
| **Scalability** | Excellent | Good | Good | Excellent |
| **Learning Curve** | Moderate | Low | High | Moderate |
| **Resume Capability** | Built-in | Built-in | Limited | Good |
| **Error Handling** | Process-level | Rule-level | Complex | Good |

**Key Nextflow Features Used**:

```groovy
// DSL2 modular imports
include { VARIANT_CALL } from './modules/variant_call'

// Process-level container specification
process ALLELE_CALL {
    container 'pgx-python:1.0.0'
}

// Built-in reporting
report.enabled = true
timeline.enabled = true
trace.enabled = true
```

**Why Not Snakemake?**
- Nextflow has stronger clinical genomics adoption
- Better Docker integration out-of-the-box
- Process isolation matches our module contract design
- Native support for cloud execution (AWS, GCP, Azure)

---

## 3. Variant Calling Tools

### GATK HaplotypeCaller

**Decision**: Use GATK HaplotypeCaller for variant calling in PGx regions.

**Rationale**:

| Tool | Accuracy | Speed | Clinical Validation | PGx Suitability |
|------|----------|-------|---------------------|-----------------|
| **GATK HC** | Excellent | Moderate | Gold standard | Excellent |
| bcftools mpileup | Good | Fast | Common | Good for SNPs |
| DeepVariant | Excellent | Slow | Emerging | Untested for PGx |
| Strelka2 | Good | Fast | Somatic focus | Not designed for germline PGx |

**Why GATK HaplotypeCaller?**

1. **Clinical Familiarity**: Most clinical labs already use GATK
2. **Haplotype-aware**: Reconstructs haplotypes during calling (critical for PGx)
3. **Extensive Validation**: Used in major sequencing initiatives (GIAB, GeT-RM)
4. **gVCF Support**: Enables joint calling if batch processing needed
5. **Active Maintenance**: Broad Institute continues development

**Configuration Choices**:

```bash
gatk HaplotypeCaller \
    -R hg38.fa \
    -I sample.bam \
    -L pgx_regions.bed \           # Target only PGx regions
    --min-base-quality-score 20 \  # Filter low-quality bases
    --minimum-mapping-quality 20   # Exclude poorly mapped reads
```

**Why Not bcftools?**
- GATK provides haplotype reconstruction which aids diplotype inference
- Better handling of complex indels
- More extensive validation in clinical settings

---

## 4. Python Packages

### 4.1 Pydantic v2

**Decision**: Use Pydantic v2 for all data models.

**Rationale**:

```python
# Type-safe, self-documenting models
class GeneResult(BaseModel):
    gene: str = Field(..., description="HGNC gene symbol")
    diplotype: str
    phenotype: str
    confidence: ConfidenceMetrics
```

**Benefits**:
- **Validation at runtime**: Catches data errors early
- **JSON serialization**: Built-in `model_dump_json()`
- **Schema generation**: Auto-generates JSON Schema
- **IDE support**: Full autocomplete and type checking
- **Performance**: Pydantic v2 is 5-50x faster than v1

**Why Pydantic over dataclasses?**
- Built-in validation
- JSON serialization without extra code
- Field descriptions for documentation
- Union types handled correctly

### 4.2 pysam

**Decision**: Use pysam for BAM/CRAM parsing.

**Rationale**:

| Feature | pysam | pybam | Raw samtools |
|---------|-------|-------|--------------|
| **Pythonic API** | ✅ | ✅ | ❌ |
| **CRAM Support** | ✅ | ❌ | ✅ |
| **Header Access** | ✅ | Limited | Manual |
| **Pileup** | ✅ | ❌ | Manual |
| **Maintenance** | Active | Stale | N/A |

**Usage in Pipeline**:
```python
with pysam.AlignmentFile(bam_path, "rb") as bam:
    # Extract sample ID from read group
    sample_id = bam.header.get("RG", [{}])[0].get("SM")

    # Calculate coverage over PGx regions
    for pileup_col in bam.pileup(chrom, start, end):
        depths.append(pileup_col.n)
```

### 4.3 cyvcf2

**Decision**: Use cyvcf2 as primary VCF parser (with pysam fallback).

**Rationale**:

| Feature | cyvcf2 | PyVCF | pysam.VariantFile |
|---------|--------|-------|-------------------|
| **Speed** | Fastest | Slow | Fast |
| **Memory** | Low | High | Medium |
| **API** | Clean | Dated | Good |
| **Genotype Access** | Easy | Complex | Easy |

**Performance Comparison** (10,000 variants):
- cyvcf2: ~0.5 seconds
- pysam: ~1.2 seconds
- PyVCF: ~5+ seconds

**Code Pattern**:
```python
# cyvcf2 is primary, pysam is fallback
if cyvcf2_available:
    return parse_vcf_cyvcf2(vcf_path)
else:
    return parse_vcf_pysam(vcf_path)
```

### 4.4 jsonschema

**Decision**: Use jsonschema for output validation.

**Rationale**:
- **Standards-compliant**: JSON Schema Draft 2020-12
- **Clinical requirement**: Validates outputs before release
- **Extensible**: Schema can evolve with version control
- **Well-maintained**: Active development, good documentation

---

## 5. Allele Calling Algorithm

### VCF-based Star-Allele Matching

**Decision**: Match observed variants against known star-allele definitions.

**Algorithm**:

```
1. Load allele definitions (star allele → defining variants)
2. Parse VCF for sample
3. For each gene:
   a. Find all observed variants at defining positions
   b. Match observed variants to allele definitions
   c. Count copies of each matched allele (from genotype)
   d. Assign diplotype (allele1/allele2)
   e. Default to reference allele (*1) if no variants found
4. Calculate confidence based on depth/quality
```

**Why This Approach?**

| Approach | Pros | Cons | Our Choice |
|----------|------|------|------------|
| **VCF-based matching** | Simple, transparent, auditable | Misses CNVs | ✅ Mode A (Group 1) |
| **Probabilistic (Beagle/PHASE)** | Handles phase | Complexity, speed | ❌ Not needed for most PGx |
| **BAM-based (Aldy/Cyrius)** | Handles CNVs | Black-box, complex | ✅ Mode B (CYP2D6, future) |

**Handling Ambiguity**:

```python
# If multiple diplotypes equally likely
if len(possible_diplotypes) > 1:
    call_quality = "ambiguous"
    flags.append(ConfidenceFlag.AMBIGUOUS_CALL)
```

**Reference Allele Assignment**:
```python
# No variants found = homozygous reference
if not matched_alleles:
    diplotype = f"{ref_allele}/{ref_allele}"  # e.g., "*1/*1"
```

---

## 6. Phenotype Translation

### CPIC-Based Mapping

**Decision**: Use CPIC (Clinical Pharmacogenetics Implementation Consortium) guidelines as the authoritative source.

**Rationale**:

| Source | Clinical Adoption | Update Frequency | Regulatory Standing |
|--------|-------------------|------------------|---------------------|
| **CPIC** | Gold standard | Regular | FDA-recognized |
| DPWG | Common in Europe | Regular | EU standard |
| PharmGKB | Research reference | Frequent | Informational |
| Custom | - | - | Not defensible |

**Translation Chain**:

```
Diplotype (*1/*17)
    ↓
Activity Score (1.0 + 1.5 = 2.5)
    ↓
Phenotype (Rapid Metabolizer)
    ↓
Prescribing Implication (CPIC guideline text)
```

**Activity Score System** (CYP2C19 example):

| Allele | Function | Activity Value |
|--------|----------|----------------|
| *1 | Normal | 1.0 |
| *2 | No function | 0.0 |
| *17 | Increased | 1.5 |

**Why Activity Scores?**
- Consistent across genes
- Handles novel combinations mathematically
- CPIC-standardized approach
- Auditable calculation

---

## 7. Confidence Scoring

### Composite Score Algorithm

**Decision**: Calculate confidence as weighted composite of depth, quality, and completeness.

**Formula**:
```
Confidence = 0.40 × DepthScore + 0.35 × QualityScore + 0.25 × CompletenessScore
```

**Component Definitions**:

| Component | Calculation | Rationale |
|-----------|-------------|-----------|
| **DepthScore** | min(1.0, mean_depth / 30) | Higher depth = more reliable |
| **QualityScore** | mean(min(1.0, QUAL/100)) | Higher QUAL = better call |
| **CompletenessScore** | variants_found / variants_expected | Complete evidence better |

**Weight Justification**:
- **Depth (0.40)**: Most critical for PGx accuracy
- **Quality (0.35)**: Variant call reliability
- **Completeness (0.25)**: Evidence completeness

**Quality Flags**:

| Flag | Trigger | Clinical Action |
|------|---------|-----------------|
| LOW_DEPTH | Any position < 10x | Review recommended |
| LOW_QUALITY | Any QUAL < 20 | Verify variants |
| INCOMPLETE_ALLELE | Missing defining variants | Confirmatory testing |
| AMBIGUOUS_CALL | Multiple equal diplotypes | Manual review |

---

## 8. Data Architecture

### Layered Interpretation Chain

**Decision**: Structure data to preserve every layer of the interpretation chain.

```
Raw Variants → Star Alleles → Diplotype → Function → Phenotype → Implication
```

**Why Preserve All Layers?**

1. **Auditing**: Trace any result back to raw evidence
2. **Regulatory**: CAP/CLIA require explainable results
3. **Extensibility**: Different consumers need different layers
4. **Debugging**: Identify where interpretation diverges

**Directory Structure Rationale**:

```
data/
├── allele_definitions/     # Star allele → variant mappings
│   └── CYP2C19.json       # One file per gene (modular updates)
├── phenotype_maps/         # Diplotype → phenotype
│   └── CYP2C19_phenotype.json
├── prescribing_implications/  # Phenotype → clinical text
│   └── CYP2C19_implications.json
└── gene_tiers.json        # Gene classification config
```

**Why Separate Files?**
- Independent versioning per gene
- Easier CPIC guideline updates
- Clear separation of concerns
- Simpler testing

---

## 9. Containerization Strategy

### Docker with Pinned Versions

**Decision**: One Docker image per logical component, with version pinning.

**Images**:

| Image | Purpose | Base |
|-------|---------|------|
| `pgx-python:1.0.0` | Python logic | python:3.11-slim |
| `broadinstitute/gatk:4.5.0.0` | Variant calling | (official) |
| `pgx-hla:1.0.0` | HLA typing | python:3.11-slim |

**Why Separate Images?**
- Independent updates
- Smaller image sizes
- Clear dependency boundaries
- Easier security auditing

**Version Pinning Rules**:

```groovy
// BAD: Never use :latest
container = 'broadinstitute/gatk:latest'

// GOOD: Pin specific version
container = 'broadinstitute/gatk:4.5.0.0'

// BEST: Pin with digest (immutable)
container = 'broadinstitute/gatk:4.5.0.0@sha256:abc123...'
```

**Rationale for Pinning**:
- **Reproducibility**: Same results 6 months later
- **Debugging**: Know exact software state
- **Regulatory**: Required for clinical validation
- **Security**: Audit specific versions

---

## 10. Testing Framework

### Test Agent Architecture

**Decision**: Implement specialized test agents for different testing dimensions.

| Agent | Purpose | When to Run |
|-------|---------|-------------|
| **pipeline_evaluator** | End-to-end correctness | After pipeline changes |
| **api_tester** | Contract compliance | After schema/model changes |
| **workflow_optimizer** | Performance analysis | Before releases |
| **performance_benchmarker** | Regression detection | CI/CD |
| **test_results_analyzer** | Unified reporting | Always |

**Why Test Agents vs Traditional Tests?**

| Traditional pytest | Test Agents |
|-------------------|-------------|
| Unit-focused | Holistic validation |
| Fast, isolated | Integration-aware |
| Developer-centric | Quality-centric |
| Code coverage | Correctness coverage |

**Both are used**: pytest for unit tests, agents for system validation.

---

## Summary: Key Architectural Decisions

| Decision | Choice | Primary Rationale |
|----------|--------|-------------------|
| **Language** | Python 3.11+ | Type safety, ecosystem |
| **Orchestration** | Nextflow DSL2 | Clinical adoption, Docker support |
| **Variant Calling** | GATK HC | Clinical validation, haplotype-aware |
| **Data Models** | Pydantic v2 | Validation, serialization |
| **VCF Parsing** | cyvcf2 | Performance |
| **Allele Calling** | VCF-based matching | Transparency, auditability |
| **Phenotype Source** | CPIC | Regulatory compliance |
| **Confidence** | Composite scoring | Multi-factor assessment |
| **Data Structure** | Layered chain | Auditability |
| **Containers** | Docker, pinned | Reproducibility |

---

## References

1. CPIC Guidelines: https://cpicpgx.org
2. PharmVar: https://www.pharmvar.org
3. GATK Best Practices: https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows
4. Nextflow Documentation: https://www.nextflow.io/docs/latest/
5. Pydantic Documentation: https://docs.pydantic.dev/latest/
6. pysam Documentation: https://pysam.readthedocs.io/
7. cyvcf2 Documentation: https://brentp.github.io/cyvcf2/
