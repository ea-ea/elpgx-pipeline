"""
Validation Framework for ElPGx Pipeline.

Provides tools for:
- Testing against known PGx reference samples (GeT-RM, etc.)
- Calculating concordance metrics
- Generating validation reports

Reference Datasets:
- GeT-RM: CDC Genetic Testing Reference Materials
- Platinum Genomes: Illumina high-confidence calls
- 1000 Genomes: Population diversity validation
"""

import json
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import Optional


class ValidationStatus(str, Enum):
    """Validation result status."""
    CONCORDANT = "concordant"
    DISCORDANT = "discordant"
    PARTIAL = "partial"  # One allele matches
    NOT_EVALUATED = "not_evaluated"


@dataclass
class ReferenceSample:
    """A reference sample with known genotypes."""
    sample_id: str
    source: str  # "GeT-RM", "Platinum", "1000G", etc.
    known_genotypes: dict[str, str]  # gene -> diplotype
    known_phenotypes: dict[str, str]  # gene -> phenotype (optional)
    coriell_id: Optional[str] = None
    population: Optional[str] = None
    notes: Optional[str] = None


@dataclass
class GeneValidationResult:
    """Validation result for a single gene."""
    gene: str
    expected_diplotype: str
    called_diplotype: str
    expected_phenotype: Optional[str]
    called_phenotype: Optional[str]
    status: ValidationStatus
    allele1_match: bool
    allele2_match: bool
    confidence_score: float
    notes: str = ""


@dataclass
class SampleValidationResult:
    """Validation result for a complete sample."""
    sample_id: str
    source: str
    timestamp: datetime
    gene_results: list[GeneValidationResult] = field(default_factory=list)

    @property
    def concordance_rate(self) -> float:
        """Calculate overall concordance rate."""
        if not self.gene_results:
            return 0.0
        concordant = sum(1 for r in self.gene_results if r.status == ValidationStatus.CONCORDANT)
        return concordant / len(self.gene_results)

    @property
    def genes_evaluated(self) -> int:
        return len(self.gene_results)

    @property
    def genes_concordant(self) -> int:
        return sum(1 for r in self.gene_results if r.status == ValidationStatus.CONCORDANT)

    @property
    def genes_discordant(self) -> int:
        return sum(1 for r in self.gene_results if r.status == ValidationStatus.DISCORDANT)


@dataclass
class ValidationReport:
    """Complete validation report across multiple samples."""
    pipeline_version: str
    validation_date: datetime
    reference_dataset: str
    sample_results: list[SampleValidationResult] = field(default_factory=list)

    @property
    def overall_concordance(self) -> float:
        """Calculate overall concordance across all samples and genes."""
        total_genes = sum(r.genes_evaluated for r in self.sample_results)
        total_concordant = sum(r.genes_concordant for r in self.sample_results)
        return total_concordant / total_genes if total_genes > 0 else 0.0

    @property
    def per_gene_concordance(self) -> dict[str, float]:
        """Calculate concordance rate per gene."""
        gene_counts: dict[str, dict[str, int]] = {}

        for sample in self.sample_results:
            for gene_result in sample.gene_results:
                if gene_result.gene not in gene_counts:
                    gene_counts[gene_result.gene] = {"total": 0, "concordant": 0}
                gene_counts[gene_result.gene]["total"] += 1
                if gene_result.status == ValidationStatus.CONCORDANT:
                    gene_counts[gene_result.gene]["concordant"] += 1

        return {
            gene: counts["concordant"] / counts["total"]
            for gene, counts in gene_counts.items()
            if counts["total"] > 0
        }

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return {
            "pipeline_version": self.pipeline_version,
            "validation_date": self.validation_date.isoformat(),
            "reference_dataset": self.reference_dataset,
            "summary": {
                "total_samples": len(self.sample_results),
                "overall_concordance": f"{self.overall_concordance:.2%}",
                "per_gene_concordance": {
                    gene: f"{rate:.2%}"
                    for gene, rate in self.per_gene_concordance.items()
                },
            },
            "sample_results": [
                {
                    "sample_id": sr.sample_id,
                    "source": sr.source,
                    "concordance_rate": f"{sr.concordance_rate:.2%}",
                    "genes_evaluated": sr.genes_evaluated,
                    "genes_concordant": sr.genes_concordant,
                    "genes_discordant": sr.genes_discordant,
                    "gene_details": [
                        {
                            "gene": gr.gene,
                            "expected": gr.expected_diplotype,
                            "called": gr.called_diplotype,
                            "status": gr.status.value,
                            "confidence": gr.confidence_score,
                            "notes": gr.notes,
                        }
                        for gr in sr.gene_results
                    ],
                }
                for sr in self.sample_results
            ],
        }


def normalize_diplotype(diplotype: str) -> tuple[str, str]:
    """
    Normalize diplotype for comparison.

    Handles:
    - Sorting alleles (*1/*2 == *2/*1)
    - Removing whitespace
    - Case normalization

    Returns:
        Tuple of (allele1, allele2) in sorted order.
    """
    diplotype = diplotype.strip()

    # Handle various separators
    for sep in ["/", "|", " "]:
        if sep in diplotype:
            parts = diplotype.split(sep)
            if len(parts) == 2:
                alleles = sorted([p.strip() for p in parts])
                return alleles[0], alleles[1]

    # Single allele (hemizygous or error)
    return diplotype, diplotype


def compare_diplotypes(expected: str, called: str) -> tuple[ValidationStatus, bool, bool]:
    """
    Compare expected and called diplotypes.

    Returns:
        Tuple of (status, allele1_match, allele2_match).
    """
    exp_a1, exp_a2 = normalize_diplotype(expected)
    call_a1, call_a2 = normalize_diplotype(called)

    # Exact match
    if (exp_a1, exp_a2) == (call_a1, call_a2):
        return ValidationStatus.CONCORDANT, True, True

    # Check partial matches
    exp_set = {exp_a1, exp_a2}
    call_set = {call_a1, call_a2}

    matches = exp_set & call_set

    if len(matches) == 2:
        return ValidationStatus.CONCORDANT, True, True
    elif len(matches) == 1:
        # One allele matches
        a1_match = exp_a1 in call_set or call_a1 in exp_set
        a2_match = exp_a2 in call_set or call_a2 in exp_set
        return ValidationStatus.PARTIAL, a1_match, a2_match
    else:
        return ValidationStatus.DISCORDANT, False, False


def validate_sample(
    reference: ReferenceSample,
    pipeline_results: dict,  # Parsed JSON report
) -> SampleValidationResult:
    """
    Validate pipeline results against a reference sample.

    Args:
        reference: Reference sample with known genotypes.
        pipeline_results: Parsed JSON output from pipeline.

    Returns:
        SampleValidationResult with per-gene comparisons.
    """
    result = SampleValidationResult(
        sample_id=reference.sample_id,
        source=reference.source,
        timestamp=datetime.utcnow(),
    )

    # Extract results from pipeline output
    called_results = {}
    for gene_result in pipeline_results.get("results", []):
        gene = gene_result.get("gene")
        called_results[gene] = {
            "diplotype": gene_result.get("diplotype"),
            "phenotype": gene_result.get("phenotype"),
            "confidence": gene_result.get("confidence", {}).get("score", 0),
        }

    # Compare each known genotype
    for gene, expected_diplotype in reference.known_genotypes.items():
        called = called_results.get(gene, {})
        called_diplotype = called.get("diplotype", "NOT_CALLED")
        called_phenotype = called.get("phenotype")
        confidence = called.get("confidence", 0)

        expected_phenotype = reference.known_phenotypes.get(gene)

        if called_diplotype == "NOT_CALLED":
            status = ValidationStatus.NOT_EVALUATED
            a1_match, a2_match = False, False
            notes = "Gene not called by pipeline"
        else:
            status, a1_match, a2_match = compare_diplotypes(expected_diplotype, called_diplotype)
            notes = ""
            if status == ValidationStatus.PARTIAL:
                notes = "One allele matches"
            elif status == ValidationStatus.DISCORDANT:
                notes = "No alleles match"

        gene_validation = GeneValidationResult(
            gene=gene,
            expected_diplotype=expected_diplotype,
            called_diplotype=called_diplotype,
            expected_phenotype=expected_phenotype,
            called_phenotype=called_phenotype,
            status=status,
            allele1_match=a1_match,
            allele2_match=a2_match,
            confidence_score=confidence,
            notes=notes,
        )
        result.gene_results.append(gene_validation)

    return result


# === Reference Datasets ===

# GeT-RM Reference Samples (CDC Genetic Testing Reference Materials)
# These are consensus genotypes from multiple laboratories
GETRM_SAMPLES = [
    ReferenceSample(
        sample_id="NA12878",
        source="GeT-RM",
        coriell_id="NA12878",
        population="CEU",
        known_genotypes={
            "CYP2C19": "*1/*2",
            "CYP2C9": "*1/*1",
            "CYP2D6": "*1/*4",
            "VKORC1": "rs9923231 C/T",
            "SLCO1B1": "*1/*1",
            "TPMT": "*1/*1",
            "DPYD": "*1/*1",
        },
        known_phenotypes={
            "CYP2C19": "Intermediate Metabolizer",
            "CYP2C9": "Normal Metabolizer",
            "CYP2D6": "Intermediate Metabolizer",
        },
    ),
    ReferenceSample(
        sample_id="NA19685",
        source="GeT-RM",
        coriell_id="NA19685",
        population="MXL",
        known_genotypes={
            "CYP2C19": "*1/*1",
            "CYP2C9": "*1/*3",
            "CYP2D6": "*2/*2",
        },
        known_phenotypes={
            "CYP2C19": "Normal Metabolizer",
            "CYP2C9": "Intermediate Metabolizer",
            "CYP2D6": "Normal Metabolizer",
        },
    ),
    ReferenceSample(
        sample_id="HG00436",
        source="1000G",
        population="CHS",
        known_genotypes={
            "CYP2C19": "*2/*3",
            "CYP2C9": "*1/*1",
        },
        known_phenotypes={
            "CYP2C19": "Poor Metabolizer",
            "CYP2C9": "Normal Metabolizer",
        },
    ),
]


def load_reference_samples(reference_file: Path) -> list[ReferenceSample]:
    """
    Load reference samples from JSON file.

    File format:
    {
        "source": "GeT-RM",
        "samples": [
            {
                "sample_id": "NA12878",
                "coriell_id": "NA12878",
                "population": "CEU",
                "known_genotypes": {"CYP2C19": "*1/*2", ...},
                "known_phenotypes": {"CYP2C19": "Intermediate Metabolizer", ...}
            },
            ...
        ]
    }
    """
    with open(reference_file) as f:
        data = json.load(f)

    source = data.get("source", "Unknown")
    samples = []

    for sample_data in data.get("samples", []):
        sample = ReferenceSample(
            sample_id=sample_data["sample_id"],
            source=source,
            coriell_id=sample_data.get("coriell_id"),
            population=sample_data.get("population"),
            known_genotypes=sample_data.get("known_genotypes", {}),
            known_phenotypes=sample_data.get("known_phenotypes", {}),
            notes=sample_data.get("notes"),
        )
        samples.append(sample)

    return samples


def run_validation(
    reference_samples: list[ReferenceSample],
    results_dir: Path,
    pipeline_version: str,
) -> ValidationReport:
    """
    Run validation across multiple reference samples.

    Args:
        reference_samples: List of reference samples.
        results_dir: Directory containing pipeline output JSON files.
        pipeline_version: Version of pipeline being validated.

    Returns:
        ValidationReport with all results.
    """
    report = ValidationReport(
        pipeline_version=pipeline_version,
        validation_date=datetime.utcnow(),
        reference_dataset=reference_samples[0].source if reference_samples else "Unknown",
    )

    for ref_sample in reference_samples:
        # Find corresponding pipeline output
        result_file = results_dir / f"{ref_sample.sample_id}_pgx_report.json"

        if not result_file.exists():
            print(f"Warning: No results found for {ref_sample.sample_id}")
            continue

        with open(result_file) as f:
            pipeline_results = json.load(f)

        sample_result = validate_sample(ref_sample, pipeline_results)
        report.sample_results.append(sample_result)

    return report


def write_validation_report(report: ValidationReport, output_path: Path) -> None:
    """Write validation report to JSON file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as f:
        json.dump(report.to_dict(), f, indent=2)


# === Performance Metrics ===

@dataclass
class PerformanceMetrics:
    """Performance metrics for pipeline validation."""
    sensitivity: float  # True positives / (True positives + False negatives)
    specificity: float  # True negatives / (True negatives + False positives)
    accuracy: float     # (TP + TN) / Total
    concordance: float  # Exact match rate
    no_call_rate: float # Proportion of genes with no call


def calculate_metrics(report: ValidationReport) -> PerformanceMetrics:
    """
    Calculate performance metrics from validation report.

    For PGx calling, we define:
    - True Positive: Correct non-reference call
    - True Negative: Correct reference (*1/*1) call
    - False Positive: Incorrect non-reference call
    - False Negative: Missing variant (called *1/*1 when should be variant)
    """
    tp = tn = fp = fn = 0
    total_evaluated = 0
    no_calls = 0

    for sample in report.sample_results:
        for gene in sample.gene_results:
            if gene.status == ValidationStatus.NOT_EVALUATED:
                no_calls += 1
                continue

            total_evaluated += 1

            # Is expected a variant (non-reference)?
            exp_a1, exp_a2 = normalize_diplotype(gene.expected_diplotype)
            is_expected_variant = not (exp_a1 == "*1" and exp_a2 == "*1")

            # Is called a variant?
            call_a1, call_a2 = normalize_diplotype(gene.called_diplotype)
            is_called_variant = not (call_a1 == "*1" and call_a2 == "*1")

            if gene.status == ValidationStatus.CONCORDANT:
                if is_expected_variant:
                    tp += 1
                else:
                    tn += 1
            else:
                if is_expected_variant and not is_called_variant:
                    fn += 1  # Missed a variant
                elif not is_expected_variant and is_called_variant:
                    fp += 1  # Called variant when reference
                else:
                    # Both are variants but discordant
                    fp += 1  # Wrong variant call

    total = tp + tn + fp + fn
    total_with_nocalls = total + no_calls

    return PerformanceMetrics(
        sensitivity=tp / (tp + fn) if (tp + fn) > 0 else 0.0,
        specificity=tn / (tn + fp) if (tn + fp) > 0 else 0.0,
        accuracy=(tp + tn) / total if total > 0 else 0.0,
        concordance=report.overall_concordance,
        no_call_rate=no_calls / total_with_nocalls if total_with_nocalls > 0 else 0.0,
    )
