"""
Performance Benchmarker Agent

Measures latency, throughput, memory/CPU, and cost proxies.
Produces reproducible benchmarks.

Deliverables:
- Benchmark suite + command to run
- Baseline metrics
- Regression detection thresholds
"""

import json
import statistics
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable


@dataclass
class BenchmarkResult:
    """Result of a single benchmark run."""

    name: str
    iterations: int
    mean_time_ms: float
    std_dev_ms: float
    min_time_ms: float
    max_time_ms: float
    throughput_per_sec: float | None = None
    memory_peak_mb: float | None = None
    passed_threshold: bool = True
    threshold_ms: float | None = None


@dataclass
class BenchmarkSuite:
    """Collection of benchmark results."""

    suite_name: str
    results: list[BenchmarkResult]
    total_time_seconds: float
    passed: int
    failed: int


class PerformanceBenchmarker:
    """
    Performance benchmarking framework for PGx Pipeline.

    Measures execution time of individual components and detects regressions.
    """

    def __init__(
        self,
        baseline_file: Path | None = None,
        threshold_multiplier: float = 1.5,
    ):
        self.baseline_file = baseline_file
        self.threshold_multiplier = threshold_multiplier
        self.baseline: dict[str, float] = {}
        self.results: list[BenchmarkResult] = []

        if baseline_file and baseline_file.exists():
            self._load_baseline()

    def _load_baseline(self) -> None:
        """Load baseline metrics from file."""
        if self.baseline_file and self.baseline_file.exists():
            with open(self.baseline_file) as f:
                data = json.load(f)
                self.baseline = data.get("benchmarks", {})

    def save_baseline(self, output_file: Path) -> None:
        """Save current results as new baseline."""
        baseline_data = {
            "benchmarks": {r.name: r.mean_time_ms for r in self.results},
            "threshold_multiplier": self.threshold_multiplier,
        }
        output_file.parent.mkdir(parents=True, exist_ok=True)
        with open(output_file, "w") as f:
            json.dump(baseline_data, f, indent=2)

    def benchmark(
        self,
        name: str,
        func: Callable[[], Any],
        iterations: int = 10,
        warmup: int = 2,
    ) -> BenchmarkResult:
        """
        Benchmark a function.

        Args:
            name: Benchmark name.
            func: Function to benchmark (takes no arguments).
            iterations: Number of timed iterations.
            warmup: Number of warmup iterations (not timed).

        Returns:
            BenchmarkResult with timing statistics.
        """
        # Warmup
        for _ in range(warmup):
            func()

        # Timed iterations
        times_ms = []
        for _ in range(iterations):
            start = time.perf_counter()
            func()
            end = time.perf_counter()
            times_ms.append((end - start) * 1000)

        mean_time = statistics.mean(times_ms)
        std_dev = statistics.stdev(times_ms) if len(times_ms) > 1 else 0.0
        min_time = min(times_ms)
        max_time = max(times_ms)
        throughput = 1000 / mean_time if mean_time > 0 else None

        # Check threshold
        threshold = self.baseline.get(name)
        if threshold:
            threshold *= self.threshold_multiplier
            passed = mean_time <= threshold
        else:
            threshold = None
            passed = True

        result = BenchmarkResult(
            name=name,
            iterations=iterations,
            mean_time_ms=mean_time,
            std_dev_ms=std_dev,
            min_time_ms=min_time,
            max_time_ms=max_time,
            throughput_per_sec=throughput,
            passed_threshold=passed,
            threshold_ms=threshold,
        )

        self.results.append(result)
        return result

    def benchmark_json_parsing(self, json_file: Path) -> BenchmarkResult:
        """Benchmark JSON file parsing."""

        def parse_json():
            with open(json_file) as f:
                json.load(f)

        return self.benchmark(f"json_parse_{json_file.name}", parse_json)

    def benchmark_vcf_parsing(
        self,
        vcf_file: Path,
        iterations: int = 5,
    ) -> BenchmarkResult:
        """Benchmark VCF file parsing."""
        sys.path.insert(0, str(Path(__file__).parent.parent.parent))

        try:
            from elpgx.vcf_parser import parse_vcf

            def parse():
                parse_vcf(vcf_file)

            return self.benchmark(f"vcf_parse_{vcf_file.name}", parse, iterations)
        except ImportError:
            return BenchmarkResult(
                name=f"vcf_parse_{vcf_file.name}",
                iterations=0,
                mean_time_ms=0,
                std_dev_ms=0,
                min_time_ms=0,
                max_time_ms=0,
                passed_threshold=False,
            )

    def benchmark_allele_calling(
        self,
        allele_defs_dir: Path,
        iterations: int = 10,
    ) -> BenchmarkResult:
        """Benchmark allele definition loading."""
        sys.path.insert(0, str(Path(__file__).parent.parent.parent))

        try:
            from elpgx.allele_caller import load_gene_definitions

            def load():
                load_gene_definitions(allele_defs_dir)

            return self.benchmark("allele_definitions_load", load, iterations)
        except ImportError:
            return BenchmarkResult(
                name="allele_definitions_load",
                iterations=0,
                mean_time_ms=0,
                std_dev_ms=0,
                min_time_ms=0,
                max_time_ms=0,
                passed_threshold=False,
            )

    def benchmark_model_serialization(self, iterations: int = 100) -> BenchmarkResult:
        """Benchmark Pydantic model serialization."""
        sys.path.insert(0, str(Path(__file__).parent.parent.parent))

        try:
            from elpgx.models import (
                ConfidenceMetrics,
                GeneResult,
                GeneTier,
                BiomarkerType,
                PgxReport,
                PipelineMeta,
                Variant,
            )

            # Create sample data
            variant = Variant(
                rsid="rs12248560",
                chrom="chr10",
                pos=94762706,
                ref="C",
                alt="T",
                genotype="0/1",
                quality=38.5,
                depth=42,
            )

            confidence = ConfidenceMetrics(score=0.95, flags=[])

            result = GeneResult(
                gene="CYP2C19",
                gene_tier=GeneTier.GROUP1,
                biomarker_type=BiomarkerType.METABOLIC_ENZYME,
                diplotype="*1/*17",
                functional_status="Normal Function / Increased Function",
                phenotype="Rapid Metabolizer",
                prescribing_implication="Consider alternative.",
                confidence=confidence,
                variants=[variant],
            )

            meta = PipelineMeta(
                sample_id="TEST",
                pipeline_version="1.0.0",
                pgx_db_version="CPIC v1.22",
            )

            report = PgxReport(meta=meta, results=[result])

            def serialize():
                report.to_json()

            return self.benchmark("model_serialization", serialize, iterations)
        except ImportError:
            return BenchmarkResult(
                name="model_serialization",
                iterations=0,
                mean_time_ms=0,
                std_dev_ms=0,
                min_time_ms=0,
                max_time_ms=0,
                passed_threshold=False,
            )

    def run_suite(self, data_dir: Path) -> BenchmarkSuite:
        """Run all benchmarks."""
        start_time = time.time()

        # Run benchmarks
        self.benchmark_model_serialization()
        self.benchmark_allele_calling(data_dir / "allele_definitions")

        # Benchmark JSON files
        for json_file in (data_dir / "allele_definitions").glob("*.json"):
            self.benchmark_json_parsing(json_file)
            break  # Just one file for benchmark

        total_time = time.time() - start_time
        passed = sum(1 for r in self.results if r.passed_threshold)
        failed = len(self.results) - passed

        return BenchmarkSuite(
            suite_name="pgx_pipeline_benchmarks",
            results=self.results,
            total_time_seconds=total_time,
            passed=passed,
            failed=failed,
        )


def main():
    """Entry point for performance benchmarker."""
    import argparse

    parser = argparse.ArgumentParser(description="PGx Pipeline Performance Benchmarker")
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path("data"),
        help="Path to data directory",
    )
    parser.add_argument(
        "--baseline",
        type=Path,
        default=Path("tests/benchmark_baseline.json"),
        help="Path to baseline file",
    )
    parser.add_argument(
        "--save-baseline",
        action="store_true",
        help="Save results as new baseline",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=1.5,
        help="Regression threshold multiplier (default: 1.5x)",
    )
    args = parser.parse_args()

    benchmarker = PerformanceBenchmarker(
        baseline_file=args.baseline if args.baseline.exists() else None,
        threshold_multiplier=args.threshold,
    )

    suite = benchmarker.run_suite(args.data_dir)

    print(f"Benchmark Suite: {suite.suite_name}")
    print(f"Total Time: {suite.total_time_seconds:.2f}s")
    print(f"Results: {suite.passed}/{len(suite.results)} passed")
    print()

    for result in suite.results:
        status = "PASS" if result.passed_threshold else "FAIL"
        threshold_info = f" (threshold: {result.threshold_ms:.2f}ms)" if result.threshold_ms else ""
        print(f"  [{status}] {result.name}")
        print(f"         Mean: {result.mean_time_ms:.2f}ms (std: {result.std_dev_ms:.2f}ms){threshold_info}")
        print(f"         Min/Max: {result.min_time_ms:.2f}ms / {result.max_time_ms:.2f}ms")

    if args.save_baseline:
        benchmarker.save_baseline(args.baseline)
        print(f"\nBaseline saved to {args.baseline}")

    sys.exit(0 if suite.failed == 0 else 1)


if __name__ == "__main__":
    main()
