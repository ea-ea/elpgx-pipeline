# ElPGx Pipeline Safeguards Implementation Status

**Version:** 1.1.0
**Last Updated:** 2026-02-14
**Author:** Elif Arslan

This document provides a comprehensive assessment of clinical safeguards in the ElPGx pipeline, detailing what is implemented, what is partially implemented, and what requires future development.

---

## Executive Summary

| Safeguard | Status | Priority to Fix |
|-----------|--------|-----------------|
| Defining-Variant Completeness Check | ⚠️ Partial | HIGH |
| Coverage Analysis | ⚠️ Partial | HIGH |
| CNV Flagging | ❌ Not Implemented | CRITICAL |
| Confidence Downgrading | ✅ Implemented | - |
| Explicit Limitations in Output | ⚠️ Partial | MEDIUM |

---

## 1. Defining-Variant Completeness Check

### Requirement
If a required defining variant is missing or has low coverage, the phenotype should be downgraded or flagged.

### Current Implementation Status: ⚠️ PARTIAL

### What IS Implemented

#### Completeness Score Calculation
**File:** `elpgx/confidence.py`
**Lines:** 88-107

```python
def calculate_completeness_score(
    genotyped_positions: int,
    total_defining_positions: int,
) -> float:
    """
    Calculate completeness score component.

    Ratio of genotyped defining positions to total defining positions.
    """
    if total_defining_positions == 0:
        return 1.0  # No defining variants = reference allele, fully complete

    return genotyped_positions / total_defining_positions
```

#### Incomplete Allele Flag
**File:** `elpgx/confidence.py`
**Lines:** 142-144

```python
# Check completeness
if completeness_score < 1.0:
    flags.append(ConfidenceFlag.INCOMPLETE_ALLELE)
```

#### Allele Match Completeness Tracking
**File:** `elpgx/allele_caller.py`
**Lines:** 158-167

```python
is_complete = len(matched_variants) == total_defining and total_defining > 0

return AlleleMatch(
    allele=allele_def.allele,
    function=allele_def.function,
    activity_value=allele_def.activity_value,
    matched_variants=matched_variants,
    total_defining_variants=total_defining,
    is_complete_match=is_complete,  # Tracks if all defining variants found
)
```

### What is NOT Implemented

1. **Coverage verification at non-variant positions**
   - Current: Only counts variants that were called by GATK
   - Missing: Does not query BAM to verify coverage at positions where NO variant was called
   - Impact: If a defining position has 0 reads, GATK produces no call, and pipeline incorrectly assumes reference

2. **Per-position coverage requirement**
   - Missing: No minimum depth requirement at each defining position
   - Impact: A position with 1 read is treated same as position with 100 reads

3. **Phenotype downgrade based on completeness**
   - Current: Only sets `INCOMPLETE_ALLELE` flag
   - Missing: No automatic phenotype downgrade when completeness < threshold

### Required Implementation

```python
# Pseudocode for full implementation
def check_defining_position_coverage(bam_path, defining_positions, min_depth=10):
    """
    Query BAM file to verify coverage at ALL defining positions,
    not just positions where variants were called.
    """
    uncovered_positions = []
    low_coverage_positions = []

    for chrom, pos in defining_positions:
        depth = get_depth_at_position(bam_path, chrom, pos)
        if depth == 0:
            uncovered_positions.append((chrom, pos))
        elif depth < min_depth:
            low_coverage_positions.append((chrom, pos, depth))

    return uncovered_positions, low_coverage_positions
```

---

## 2. Coverage Analysis

### Requirement
Read depth and uniformity across the gene region should be evaluated.

### Current Implementation Status: ⚠️ PARTIAL

### What IS Implemented

#### Depth Score Calculation
**File:** `elpgx/confidence.py`
**Lines:** 33-57

```python
def calculate_depth_score(
    variants: list[Variant],
    target_depth: int = 30,
) -> tuple[float, float, int]:
    """
    Calculate depth score component.
    """
    if not variants:
        return 0.0, 0.0, 0

    depths = [v.depth for v in variants]
    mean_depth = sum(depths) / len(depths)
    min_depth = min(depths)

    # Score: ratio of mean depth to target, capped at 1.0
    depth_score = min(1.0, mean_depth / target_depth)

    return depth_score, mean_depth, min_depth
```

#### Low Depth Flag
**File:** `elpgx/confidence.py`
**Lines:** 130-134

```python
# Check depth
for v in variants:
    if v.depth < config.min_depth_threshold:  # default: 10
        flags.append(ConfidenceFlag.LOW_DEPTH)
        break
```

#### Confidence Reduction for Low Depth
**File:** `elpgx/allele_caller.py`
**Lines:** 357-362

```python
# Calculate average depth for confidence assessment
if all_variants:
    avg_depth = sum(v.depth or 0 for v in all_variants) / len(all_variants)
    if avg_depth < 10:
        confidence_score = min(confidence_score, 0.60)
    elif avg_depth < 20:
        confidence_score = min(confidence_score, 0.80)
```

#### Depth Statistics in Output
**File:** `elpgx/models.py`
**Lines:** 72-76

```python
class VariantDepthStats(BaseModel):
    """Depth statistics across all variants for a gene."""
    mean_depth: float = Field(..., description="Mean read depth across defining positions")
    min_depth: int = Field(..., description="Minimum read depth at any defining position")
```

### What is NOT Implemented

1. **Gene region uniformity check**
   - Missing: No coefficient of variation or uniformity metric
   - Impact: High average depth with one position at 0x is not detected

2. **Full gene region coverage analysis**
   - Current: Only checks depth at variant positions
   - Missing: Does not assess coverage across entire gene region
   - Impact: Deletions or poorly covered regions go undetected

3. **Coverage thresholds per gene**
   - Missing: All genes use same depth threshold (10x)
   - Impact: Complex genes like CYP2D6 may need higher thresholds

### Required Implementation

```python
# Pseudocode for full implementation
def analyze_gene_coverage(bam_path, gene_region, target_depth=30):
    """
    Analyze coverage across entire gene region.
    """
    chrom, start, end = gene_region
    depths = get_depth_across_region(bam_path, chrom, start, end)

    metrics = {
        "mean_depth": np.mean(depths),
        "min_depth": np.min(depths),
        "max_depth": np.max(depths),
        "median_depth": np.median(depths),
        "std_depth": np.std(depths),
        "cv": np.std(depths) / np.mean(depths),  # Coefficient of variation
        "pct_above_10x": sum(d >= 10 for d in depths) / len(depths),
        "pct_above_30x": sum(d >= 30 for d in depths) / len(depths),
        "zero_coverage_positions": sum(d == 0 for d in depths),
    }

    # Uniformity check
    if metrics["cv"] > 0.5:  # High variability
        flags.append("UNEVEN_COVERAGE")

    return metrics
```

---

## 3. CNV Flagging

### Requirement
Copy number inference based on coverage ratio heuristics. If ambiguous, flag for review.

### Current Implementation Status: ❌ NOT IMPLEMENTED

### What IS Implemented

**Nothing.** There is no CNV detection in the current production code.

### What Exists But is NOT Integrated

#### Future CNV Module (Not Active)
**File:** `elpgx/confidence_v2.py`
**Status:** Created but NOT integrated into pipeline

```python
@dataclass
class CNVEvidence:
    """Copy number variation evidence for a gene."""
    gene: str
    estimated_copy_number: float
    confidence: float  # 0-1 confidence in CNV call
    method: str  # "read_depth", "paired_end", "split_read", "none"
    is_deletion: bool = False
    is_duplication: bool = False
```

#### Documentation of CNV Plans
**File:** `docs/FUTURE_DEVELOPMENT.md`
**Section:** Copy Number Variation Detection

### Impact of Missing CNV Detection

| Gene | Affected Alleles | Clinical Impact |
|------|------------------|-----------------|
| CYP2D6 | *5 (deletion) | **CRITICAL** - Poor Metabolizer missed, called as Normal |
| CYP2D6 | *1xN, *2xN (duplication) | **CRITICAL** - Ultrarapid Metabolizer missed |
| CYP2D6 | *13, *36 (hybrid) | HIGH - Incorrect allele assignment |
| CYP2A6 | *4 (deletion) | HIGH - No function missed |
| GSTT1/GSTM1 | null | MODERATE - Null genotype missed |

### Required Implementation

```python
# Pseudocode for CNV detection
def detect_cnv_by_depth(bam_path, gene_region, reference_regions):
    """
    Detect CNV by comparing gene depth to stable reference regions.
    """
    gene_depth = get_median_depth(bam_path, gene_region)
    ref_depth = get_median_depth(bam_path, reference_regions)

    ratio = gene_depth / ref_depth

    if ratio < 0.7:
        return CNVResult(type="deletion", copy_number=1, confidence=0.8)
    elif ratio < 0.3:
        return CNVResult(type="deletion", copy_number=0, confidence=0.9)
    elif ratio > 1.3:
        estimated_cn = round(ratio * 2)
        return CNVResult(type="duplication", copy_number=estimated_cn, confidence=0.7)
    else:
        return CNVResult(type="normal", copy_number=2, confidence=0.95)


def flag_cnv_uncertainty(cnv_result):
    """
    Flag ambiguous CNV calls for manual review.
    """
    if cnv_result.confidence < 0.7:
        return ConfidenceFlag.CNV_UNCERTAIN
    if cnv_result.type != "normal":
        return ConfidenceFlag.CNV_DETECTED
    return None
```

### Workaround for Users

Until CNV detection is implemented:
1. Use external CNV caller (Stargazer, Aldy) for CYP2D6
2. Manual review of CYP2D6 calls recommended
3. Document limitation in clinical report

---

## 4. Confidence Downgrading

### Requirement
Ambiguity should automatically reduce confidence score.

### Current Implementation Status: ✅ IMPLEMENTED

### What IS Implemented

#### Ambiguity Detection
**File:** `elpgx/allele_caller.py`
**Lines:** 350-352

```python
# Determine call quality and confidence
call_quality = "definitive"
confidence_score = 0.95
if len(allele_copies) > 2:
    call_quality = "ambiguous"
    confidence_score = 0.70
```

#### Ambiguous Call Flag
**File:** `elpgx/confidence.py`
**Lines:** 146-150

```python
# Check call quality
if call_quality == "ambiguous":
    flags.append(ConfidenceFlag.AMBIGUOUS_CALL)
elif call_quality == "no_call":
    flags.append(ConfidenceFlag.NO_CALL)
```

#### Low Depth Confidence Reduction
**File:** `elpgx/allele_caller.py`
**Lines:** 357-362

```python
if all_variants:
    avg_depth = sum(v.depth or 0 for v in all_variants) / len(all_variants)
    if avg_depth < 10:
        confidence_score = min(confidence_score, 0.60)
    elif avg_depth < 20:
        confidence_score = min(confidence_score, 0.80)
```

#### Reference Call Lower Confidence
**File:** `elpgx/allele_caller.py`
**Lines:** 281-282

```python
# No variants found -> homozygous reference
confidence=ConfidenceMetrics(score=0.90, flags=[]),  # Lower confidence for ref/ref
```

#### Low Quality Flag
**File:** `elpgx/confidence.py`
**Lines:** 136-140

```python
# Check quality
for v in variants:
    if v.quality < config.min_quality_threshold:  # default: 20
        flags.append(ConfidenceFlag.LOW_QUALITY)
        break
```

#### Phenotype Blocking Based on Confidence
**File:** `modules/report.nf`
**Lines:** 109-130 (newly added)

```python
def should_block_phenotype(gene, confidence_score, flags):
    strategy = gene_strategies.get(gene)
    if not strategy:
        return False, None, None

    for rule in strategy.get("phenotype_blocking_rules", []):
        condition = rule.get("condition", "")

        if "confidence_score <" in condition:
            threshold = float(condition.split("<")[1].strip())
            if confidence_score < threshold:
                return True, rule.get("phenotype_value"), rule.get("message")
```

### Confidence Score Ranges

| Scenario | Confidence Score | Flags Set |
|----------|------------------|-----------|
| Definitive call, good depth | 0.95 | None |
| Reference-only call | 0.90 | None |
| No variants, uncertain | 0.85 | None |
| Ambiguous (>2 allele matches) | 0.70 | AMBIGUOUS_CALL |
| Low depth (10-20x) | ≤0.80 | LOW_DEPTH |
| Very low depth (<10x) | ≤0.60 | LOW_DEPTH |
| No call possible | N/A | NO_CALL |

---

## 5. Explicit Limitations in Output

### Requirement
Report should include warnings for complex gene uncertainty.

### Current Implementation Status: ⚠️ PARTIAL

### What IS Implemented

#### Warnings Array in Report
**File:** `elpgx/models.py`
**Lines:** 185-188

```python
class PgxReport(BaseModel):
    warnings: list[str] = Field(
        default_factory=list,
        description="Pipeline-level warnings"
    )
```

#### Warning Generation in Report Builder
**File:** `elpgx/report.py`
**Lines:** 135-140

```python
warnings: list[str] = []

for call in diplotype_calls:
    gene = call.gene
    phenotype = phenotype_by_gene.get(gene)

    if phenotype is None:
        warnings.append(f"No phenotype translation available for {gene}")
```

#### CYP2D6 Strategy Document
**File:** `data/gene_strategies/CYP2D6_STRATEGY.json`

```json
{
  "known_limitations": [
    "Cannot detect *5 (gene deletion) - will incorrectly call as *1/*1 or *1/*X",
    "Cannot detect gene duplications - ultrarapid metabolizers may be missed",
    "Hybrid alleles (*13, *36) will be partially called based on SNV component only",
    "Activity score assumes diploid state without CNV",
    "Phase information not available - rare cis arrangements may be missed"
  ]
}
```

#### Phenotype Blocking Warnings
**File:** `modules/report.nf`
**Lines:** 124-126

```python
if block:
    phenotype = blocked_phenotype
    warnings.append(f"{gene}: {block_message}")
```

### What is NOT Implemented

1. **Automatic CYP2D6 CNV warning**
   - Missing: Every CYP2D6 call should warn about CNV limitation
   - Current: No automatic warning added

2. **Gene-specific limitation warnings**
   - Missing: Automatic warnings from gene strategy files
   - Current: Strategy files exist but not automatically included in output

3. **Structured limitation reporting**
   - Missing: Limitations as structured data, not just text warnings
   - Current: Only text warnings in array

### Required Implementation

```python
# Add to report generation
def add_gene_limitation_warnings(gene, strategy, warnings):
    """
    Automatically add limitation warnings for complex genes.
    """
    if strategy is None:
        return

    unsupported = strategy.get("unsupported_features", {})
    for feature, info in unsupported.items():
        if info.get("clinical_impact") in ["critical", "high"]:
            warnings.append(
                f"{gene}: {info['description']} - "
                f"Clinical impact: {info['clinical_impact']}"
            )

# Automatically add CYP2D6 warning
if gene == "CYP2D6":
    warnings.append(
        "CYP2D6: CNV detection not performed. "
        "*5 deletions and gene duplications cannot be detected. "
        "Consider confirmatory testing for clinical decisions."
    )
```

---

## Implementation Roadmap

### Phase 1: Critical (Before Clinical Use)

| Item | Priority | Effort | Impact |
|------|----------|--------|--------|
| Add automatic CYP2D6 CNV warning | CRITICAL | Low | Prevents misinterpretation |
| Implement position coverage check | HIGH | Medium | Catches coverage gaps |
| Add gene-specific limitation warnings | HIGH | Low | User awareness |

### Phase 2: Important (v1.2.0)

| Item | Priority | Effort | Impact |
|------|----------|--------|--------|
| Full gene region coverage analysis | HIGH | Medium | Better confidence |
| Read-depth CNV heuristics | HIGH | High | CYP2D6 accuracy |
| Coverage uniformity metrics | MEDIUM | Medium | Quality metrics |

### Phase 3: Enhancement (v1.3.0)

| Item | Priority | Effort | Impact |
|------|----------|--------|--------|
| Full CNV calling module | HIGH | High | Complex gene accuracy |
| Per-gene coverage thresholds | MEDIUM | Low | Gene-specific QC |
| Structured limitation reporting | LOW | Low | Better integration |

---

## Testing Checklist

Before considering safeguards complete, verify:

- [ ] Low coverage position triggers `LOW_DEPTH` flag
- [ ] Missing defining variant triggers `INCOMPLETE_ALLELE` flag
- [ ] Ambiguous calls trigger `AMBIGUOUS_CALL` flag
- [ ] Confidence < 0.50 blocks phenotype (CYP2D6)
- [ ] CYP2D6 report includes CNV limitation warning
- [ ] Provenance includes all required fields
- [ ] Warnings array populated for complex genes

---

## References

1. CPIC Guidelines: https://cpicpgx.org/guidelines/
2. PharmVar CYP2D6: https://www.pharmvar.org/gene/CYP2D6
3. Pratt VM et al. Recommendations for Clinical CYP2D6 Genotyping. J Mol Diagn. 2021.
4. Gaedigk A et al. The Pharmacogene Variation Consortium. Clin Pharmacol Ther. 2018.
