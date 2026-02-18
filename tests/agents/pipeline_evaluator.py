"""
Pipeline Evaluator Agent

Runs end-to-end pipeline validations using known-good datasets/fixtures.
Checks correctness against expected outputs.
Reports pass/fail with reproducible steps.

Deliverables:
- Evaluation plan (what "correct" means)
- Fixture inputs + expected outputs
- End-to-end test runner entry point
"""

import json
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


@dataclass
class EvaluationResult:
    """Result of a single evaluation test."""

    test_name: str
    passed: bool
    expected: Any = None
    actual: Any = None
    error_message: str | None = None
    details: dict = field(default_factory=dict)


@dataclass
class EvaluationReport:
    """Complete evaluation report."""

    total_tests: int
    passed: int
    failed: int
    results: list[EvaluationResult]
    summary: str


class PipelineEvaluator:
    """
    End-to-end pipeline evaluator.

    Validates pipeline outputs against known-good expected results.
    """

    def __init__(
        self,
        fixtures_dir: Path,
        output_dir: Path,
    ):
        self.fixtures_dir = fixtures_dir
        self.output_dir = output_dir
        self.results: list[EvaluationResult] = []

    def load_expected_output(self, gene: str) -> dict | None:
        """Load expected output for a gene from fixtures."""
        expected_file = self.fixtures_dir / f"expected_{gene}.json"
        if not expected_file.exists():
            return None
        with open(expected_file) as f:
            return json.load(f)

    def load_actual_output(self, sample_id: str) -> dict | None:
        """Load actual pipeline output."""
        output_file = self.output_dir / f"{sample_id}_pgx_report.json"
        if not output_file.exists():
            return None
        with open(output_file) as f:
            return json.load(f)

    def evaluate_diplotype(
        self,
        gene: str,
        expected: str,
        actual: str,
    ) -> EvaluationResult:
        """Evaluate diplotype call correctness."""
        # Normalize diplotype order
        def normalize_diplotype(d: str) -> str:
            parts = sorted(d.split("/"))
            return "/".join(parts)

        exp_norm = normalize_diplotype(expected)
        act_norm = normalize_diplotype(actual)
        passed = exp_norm == act_norm

        return EvaluationResult(
            test_name=f"diplotype_{gene}",
            passed=passed,
            expected=expected,
            actual=actual,
            error_message=None if passed else f"Diplotype mismatch for {gene}",
        )

    def evaluate_phenotype(
        self,
        gene: str,
        expected: str,
        actual: str,
    ) -> EvaluationResult:
        """Evaluate phenotype prediction correctness."""
        passed = expected.lower() == actual.lower()

        return EvaluationResult(
            test_name=f"phenotype_{gene}",
            passed=passed,
            expected=expected,
            actual=actual,
            error_message=None if passed else f"Phenotype mismatch for {gene}",
        )

    def evaluate_confidence_range(
        self,
        gene: str,
        score: float,
        min_expected: float = 0.5,
    ) -> EvaluationResult:
        """Evaluate confidence score is within acceptable range."""
        passed = score >= min_expected

        return EvaluationResult(
            test_name=f"confidence_{gene}",
            passed=passed,
            expected=f">= {min_expected}",
            actual=score,
            error_message=None if passed else f"Low confidence for {gene}: {score}",
        )

    def evaluate_report_structure(self, report: dict) -> EvaluationResult:
        """Evaluate report has required structure."""
        required_keys = ["meta", "results"]
        missing = [k for k in required_keys if k not in report]

        passed = len(missing) == 0

        return EvaluationResult(
            test_name="report_structure",
            passed=passed,
            expected=required_keys,
            actual=list(report.keys()),
            error_message=f"Missing keys: {missing}" if missing else None,
        )

    def evaluate_sample(
        self,
        sample_id: str,
        expected_results: dict,
    ) -> list[EvaluationResult]:
        """Evaluate all aspects of a sample's pipeline output."""
        results = []

        # Load actual output
        actual = self.load_actual_output(sample_id)
        if actual is None:
            results.append(
                EvaluationResult(
                    test_name="output_exists",
                    passed=False,
                    error_message=f"Output file not found for {sample_id}",
                )
            )
            return results

        # Check structure
        results.append(self.evaluate_report_structure(actual))

        # Check each gene result
        actual_by_gene = {r["gene"]: r for r in actual.get("results", [])}

        for gene, expected in expected_results.items():
            if gene not in actual_by_gene:
                results.append(
                    EvaluationResult(
                        test_name=f"gene_present_{gene}",
                        passed=False,
                        error_message=f"Gene {gene} not in output",
                    )
                )
                continue

            actual_gene = actual_by_gene[gene]

            # Check diplotype
            if "diplotype" in expected:
                results.append(
                    self.evaluate_diplotype(
                        gene,
                        expected["diplotype"],
                        actual_gene.get("diplotype", ""),
                    )
                )

            # Check phenotype
            if "phenotype" in expected:
                results.append(
                    self.evaluate_phenotype(
                        gene,
                        expected["phenotype"],
                        actual_gene.get("phenotype", ""),
                    )
                )

            # Check confidence
            if "confidence" in actual_gene:
                score = actual_gene["confidence"].get("score", 0)
                results.append(self.evaluate_confidence_range(gene, score))

        return results

    def run_evaluation(
        self,
        test_cases: list[dict],
    ) -> EvaluationReport:
        """
        Run full evaluation across all test cases.

        Args:
            test_cases: List of dicts with sample_id and expected results.

        Returns:
            EvaluationReport with all results.
        """
        all_results = []

        for case in test_cases:
            sample_id = case.get("sample_id")
            expected = case.get("expected", {})
            results = self.evaluate_sample(sample_id, expected)
            all_results.extend(results)

        passed = sum(1 for r in all_results if r.passed)
        failed = len(all_results) - passed

        summary = f"Evaluation complete: {passed}/{len(all_results)} tests passed"
        if failed > 0:
            summary += f" ({failed} failures)"

        return EvaluationReport(
            total_tests=len(all_results),
            passed=passed,
            failed=failed,
            results=all_results,
            summary=summary,
        )


def create_golden_fixtures(fixtures_dir: Path) -> None:
    """Create example golden fixtures for testing."""
    fixtures_dir.mkdir(parents=True, exist_ok=True)

    # Example expected output for CYP2C19 homozygous reference
    cyp2c19_homref = {
        "gene": "CYP2C19",
        "diplotype": "*1/*1",
        "phenotype": "Normal Metabolizer",
        "activity_score": 2.0,
    }

    with open(fixtures_dir / "expected_CYP2C19_homref.json", "w") as f:
        json.dump(cyp2c19_homref, f, indent=2)

    # Example expected output for CYP2C19 *1/*17
    cyp2c19_rapid = {
        "gene": "CYP2C19",
        "diplotype": "*1/*17",
        "phenotype": "Rapid Metabolizer",
        "activity_score": 2.5,
    }

    with open(fixtures_dir / "expected_CYP2C19_rapid.json", "w") as f:
        json.dump(cyp2c19_rapid, f, indent=2)

    print(f"Created golden fixtures in {fixtures_dir}")


def main():
    """Entry point for pipeline evaluator."""
    import argparse

    parser = argparse.ArgumentParser(description="PGx Pipeline Evaluator")
    parser.add_argument(
        "--fixtures",
        type=Path,
        default=Path("tests/fixtures"),
        help="Path to fixtures directory",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results"),
        help="Path to pipeline output directory",
    )
    parser.add_argument(
        "--create-fixtures",
        action="store_true",
        help="Create example golden fixtures",
    )
    args = parser.parse_args()

    if args.create_fixtures:
        create_golden_fixtures(args.fixtures)
        return

    evaluator = PipelineEvaluator(
        fixtures_dir=args.fixtures,
        output_dir=args.output,
    )

    # Define test cases
    test_cases = [
        {
            "sample_id": "test_sample",
            "expected": {
                "CYP2C19": {"diplotype": "*1/*1", "phenotype": "Normal Metabolizer"},
            },
        },
    ]

    report = evaluator.run_evaluation(test_cases)
    print(report.summary)

    for result in report.results:
        status = "PASS" if result.passed else "FAIL"
        print(f"  [{status}] {result.test_name}")
        if result.error_message:
            print(f"         {result.error_message}")

    sys.exit(0 if report.failed == 0 else 1)


if __name__ == "__main__":
    main()
