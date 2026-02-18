"""
Test Results Analyzer Agent

Aggregates results from all test agents.
Produces human-readable summary + machine-readable artifact (JSON).
Flags failures with likely root-cause clusters.

Deliverables:
- Unified test report format
- Summary output (console + artifact)
- Trend tracking hooks
"""

import json
import sys
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any


@dataclass
class TestCategory:
    """Category of test results."""

    name: str
    total: int
    passed: int
    failed: int
    skipped: int = 0
    duration_seconds: float = 0.0


@dataclass
class FailureCluster:
    """Cluster of related failures."""

    root_cause: str
    failures: list[str]
    severity: str  # "critical", "high", "medium", "low"
    suggested_action: str


@dataclass
class UnifiedTestReport:
    """Unified test report aggregating all test agents."""

    timestamp: str
    pipeline_version: str
    categories: list[TestCategory]
    total_tests: int
    total_passed: int
    total_failed: int
    total_skipped: int
    failure_clusters: list[FailureCluster]
    summary: str
    raw_results: dict[str, Any]


class TestResultsAnalyzer:
    """
    Aggregates and analyzes results from all test agents.

    Produces unified reports and identifies failure patterns.
    """

    def __init__(self, results_dir: Path):
        self.results_dir = results_dir
        self.categories: list[TestCategory] = []
        self.failures: list[dict] = []
        self.raw_results: dict[str, Any] = {}

    def load_pipeline_evaluator_results(self, results_file: Path) -> TestCategory:
        """Load results from pipeline evaluator."""
        if not results_file.exists():
            return TestCategory(
                name="pipeline_evaluator",
                total=0,
                passed=0,
                failed=0,
                skipped=1,
            )

        with open(results_file) as f:
            data = json.load(f)

        self.raw_results["pipeline_evaluator"] = data

        total = data.get("total_tests", 0)
        passed = data.get("passed", 0)
        failed = data.get("failed", 0)

        # Collect failures
        for result in data.get("results", []):
            if not result.get("passed", True):
                self.failures.append(
                    {
                        "agent": "pipeline_evaluator",
                        "test": result.get("test_name", "unknown"),
                        "error": result.get("error_message", ""),
                        "expected": result.get("expected"),
                        "actual": result.get("actual"),
                    }
                )

        return TestCategory(
            name="pipeline_evaluator",
            total=total,
            passed=passed,
            failed=failed,
        )

    def load_api_tester_results(self, results_file: Path) -> TestCategory:
        """Load results from API tester."""
        if not results_file.exists():
            return TestCategory(
                name="api_tester",
                total=0,
                passed=0,
                failed=0,
                skipped=1,
            )

        with open(results_file) as f:
            data = json.load(f)

        self.raw_results["api_tester"] = data

        total = len(data.get("results", []))
        passed = sum(1 for r in data.get("results", []) if r.get("passed"))
        failed = total - passed

        for result in data.get("results", []):
            if not result.get("passed"):
                for error in result.get("errors", []):
                    self.failures.append(
                        {
                            "agent": "api_tester",
                            "test": result.get("contract_name", "unknown"),
                            "error": error,
                        }
                    )

        return TestCategory(
            name="api_tester",
            total=total,
            passed=passed,
            failed=failed,
        )

    def load_benchmark_results(self, results_file: Path) -> TestCategory:
        """Load results from performance benchmarker."""
        if not results_file.exists():
            return TestCategory(
                name="performance_benchmarker",
                total=0,
                passed=0,
                failed=0,
                skipped=1,
            )

        with open(results_file) as f:
            data = json.load(f)

        self.raw_results["performance_benchmarker"] = data

        results = data.get("results", [])
        total = len(results)
        passed = sum(1 for r in results if r.get("passed_threshold", True))
        failed = total - passed

        for result in results:
            if not result.get("passed_threshold", True):
                self.failures.append(
                    {
                        "agent": "performance_benchmarker",
                        "test": result.get("name", "unknown"),
                        "error": f"Exceeded threshold: {result.get('mean_time_ms')}ms > {result.get('threshold_ms')}ms",
                    }
                )

        return TestCategory(
            name="performance_benchmarker",
            total=total,
            passed=passed,
            failed=failed,
            duration_seconds=data.get("total_time_seconds", 0),
        )

    def cluster_failures(self) -> list[FailureCluster]:
        """Identify patterns in failures and cluster by root cause."""
        clusters = []

        # Cluster by agent
        by_agent: dict[str, list[dict]] = {}
        for failure in self.failures:
            agent = failure.get("agent", "unknown")
            if agent not in by_agent:
                by_agent[agent] = []
            by_agent[agent].append(failure)

        # Analyze patterns
        schema_failures = [
            f for f in self.failures if "schema" in f.get("error", "").lower()
        ]
        if schema_failures:
            clusters.append(
                FailureCluster(
                    root_cause="Schema validation failures",
                    failures=[f.get("test", "") for f in schema_failures],
                    severity="high",
                    suggested_action="Review JSON schema compliance in report generation",
                )
            )

        diplotype_failures = [
            f for f in self.failures if "diplotype" in f.get("test", "").lower()
        ]
        if diplotype_failures:
            clusters.append(
                FailureCluster(
                    root_cause="Diplotype calling errors",
                    failures=[f.get("test", "") for f in diplotype_failures],
                    severity="critical",
                    suggested_action="Review allele calling logic and variant matching",
                )
            )

        phenotype_failures = [
            f for f in self.failures if "phenotype" in f.get("test", "").lower()
        ]
        if phenotype_failures:
            clusters.append(
                FailureCluster(
                    root_cause="Phenotype translation errors",
                    failures=[f.get("test", "") for f in phenotype_failures],
                    severity="high",
                    suggested_action="Review phenotype mapping tables",
                )
            )

        perf_failures = [
            f for f in self.failures if f.get("agent") == "performance_benchmarker"
        ]
        if perf_failures:
            clusters.append(
                FailureCluster(
                    root_cause="Performance regression",
                    failures=[f.get("test", "") for f in perf_failures],
                    severity="medium",
                    suggested_action="Profile slow operations and optimize",
                )
            )

        return clusters

    def generate_summary(self) -> str:
        """Generate human-readable summary."""
        total = sum(c.total for c in self.categories)
        passed = sum(c.passed for c in self.categories)
        failed = sum(c.failed for c in self.categories)
        skipped = sum(c.skipped for c in self.categories)

        lines = [
            "=" * 60,
            "PGx Pipeline Test Results Summary",
            "=" * 60,
            f"Timestamp: {datetime.utcnow().isoformat()}",
            "",
            f"Total Tests:  {total}",
            f"Passed:       {passed} ({100*passed/total:.1f}%)" if total > 0 else "Passed: 0",
            f"Failed:       {failed}",
            f"Skipped:      {skipped}",
            "",
            "By Category:",
        ]

        for cat in self.categories:
            status = "PASS" if cat.failed == 0 else "FAIL"
            lines.append(
                f"  [{status}] {cat.name}: {cat.passed}/{cat.total} passed"
            )

        if self.failures:
            lines.append("")
            lines.append("Failures:")
            for failure in self.failures[:10]:  # Show first 10
                lines.append(f"  - [{failure.get('agent')}] {failure.get('test')}")
                if failure.get("error"):
                    lines.append(f"    Error: {failure.get('error')[:80]}")

            if len(self.failures) > 10:
                lines.append(f"  ... and {len(self.failures) - 10} more")

        lines.append("")
        lines.append("=" * 60)

        return "\n".join(lines)

    def analyze(self, pipeline_version: str = "1.0.0") -> UnifiedTestReport:
        """Run complete analysis and generate unified report."""
        # Load results from each agent
        self.categories.append(
            self.load_pipeline_evaluator_results(
                self.results_dir / "pipeline_evaluator_results.json"
            )
        )
        self.categories.append(
            self.load_api_tester_results(
                self.results_dir / "api_tester_results.json"
            )
        )
        self.categories.append(
            self.load_benchmark_results(
                self.results_dir / "benchmark_results.json"
            )
        )

        # Cluster failures
        clusters = self.cluster_failures()

        # Calculate totals
        total = sum(c.total for c in self.categories)
        passed = sum(c.passed for c in self.categories)
        failed = sum(c.failed for c in self.categories)
        skipped = sum(c.skipped for c in self.categories)

        return UnifiedTestReport(
            timestamp=datetime.utcnow().isoformat(),
            pipeline_version=pipeline_version,
            categories=self.categories,
            total_tests=total,
            total_passed=passed,
            total_failed=failed,
            total_skipped=skipped,
            failure_clusters=clusters,
            summary=self.generate_summary(),
            raw_results=self.raw_results,
        )

    def save_report(self, report: UnifiedTestReport, output_file: Path) -> None:
        """Save report to JSON file."""
        output_file.parent.mkdir(parents=True, exist_ok=True)

        report_dict = {
            "timestamp": report.timestamp,
            "pipeline_version": report.pipeline_version,
            "total_tests": report.total_tests,
            "total_passed": report.total_passed,
            "total_failed": report.total_failed,
            "total_skipped": report.total_skipped,
            "categories": [
                {
                    "name": c.name,
                    "total": c.total,
                    "passed": c.passed,
                    "failed": c.failed,
                    "skipped": c.skipped,
                }
                for c in report.categories
            ],
            "failure_clusters": [
                {
                    "root_cause": fc.root_cause,
                    "failures": fc.failures,
                    "severity": fc.severity,
                    "suggested_action": fc.suggested_action,
                }
                for fc in report.failure_clusters
            ],
        }

        with open(output_file, "w") as f:
            json.dump(report_dict, f, indent=2)


def main():
    """Entry point for test results analyzer."""
    import argparse

    parser = argparse.ArgumentParser(description="PGx Pipeline Test Results Analyzer")
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=Path("tests/results"),
        help="Directory containing test agent results",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("tests/unified_test_report.json"),
        help="Output file for unified report",
    )
    parser.add_argument(
        "--version",
        type=str,
        default="1.0.0",
        help="Pipeline version",
    )
    args = parser.parse_args()

    analyzer = TestResultsAnalyzer(results_dir=args.results_dir)
    report = analyzer.analyze(pipeline_version=args.version)

    # Print summary
    print(report.summary)

    # Print failure clusters
    if report.failure_clusters:
        print("\nFailure Analysis:")
        for cluster in report.failure_clusters:
            print(f"\n  [{cluster.severity.upper()}] {cluster.root_cause}")
            print(f"    Affected: {len(cluster.failures)} tests")
            print(f"    Action: {cluster.suggested_action}")

    # Save report
    analyzer.save_report(report, args.output)
    print(f"\nUnified report saved to {args.output}")

    sys.exit(0 if report.total_failed == 0 else 1)


if __name__ == "__main__":
    main()
