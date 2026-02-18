#!/usr/bin/env python3
"""
ElPGx Validation Script

Runs validation against reference samples and generates concordance report.

Usage:
    python scripts/run_validation.py --results-dir results/ --output validation_report.json

Requirements:
    - Pipeline results in results/ directory
    - Reference samples in tests/validation/reference_samples.json
"""

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from tests.validation.validation_framework import (
    load_reference_samples,
    run_validation,
    calculate_metrics,
    write_validation_report,
    GETRM_SAMPLES,
)


def main():
    parser = argparse.ArgumentParser(
        description="Validate ElPGx pipeline against reference samples"
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        required=True,
        help="Directory containing pipeline output JSON files",
    )
    parser.add_argument(
        "--reference-file",
        type=Path,
        default=Path(__file__).parent.parent / "tests/validation/reference_samples.json",
        help="JSON file with reference samples (default: tests/validation/reference_samples.json)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("validation_report.json"),
        help="Output validation report path",
    )
    parser.add_argument(
        "--pipeline-version",
        type=str,
        default="1.1.0",
        help="Pipeline version being validated",
    )
    parser.add_argument(
        "--use-builtin",
        action="store_true",
        help="Use built-in GeT-RM samples instead of file",
    )

    args = parser.parse_args()

    print("=" * 60)
    print("ElPGx Pipeline Validation")
    print("=" * 60)

    # Load reference samples
    if args.use_builtin:
        print(f"\nUsing built-in GeT-RM reference samples ({len(GETRM_SAMPLES)} samples)")
        reference_samples = GETRM_SAMPLES
    else:
        if not args.reference_file.exists():
            print(f"Error: Reference file not found: {args.reference_file}")
            sys.exit(1)
        reference_samples = load_reference_samples(args.reference_file)
        print(f"\nLoaded {len(reference_samples)} reference samples from {args.reference_file}")

    # Check results directory
    if not args.results_dir.exists():
        print(f"Error: Results directory not found: {args.results_dir}")
        sys.exit(1)

    # Find available result files
    result_files = list(args.results_dir.glob("*_pgx_report.json"))
    print(f"Found {len(result_files)} result files in {args.results_dir}")

    # Run validation
    print("\nRunning validation...")
    report = run_validation(
        reference_samples=reference_samples,
        results_dir=args.results_dir,
        pipeline_version=args.pipeline_version,
    )

    # Calculate metrics
    metrics = calculate_metrics(report)

    # Print summary
    print("\n" + "=" * 60)
    print("VALIDATION SUMMARY")
    print("=" * 60)
    print(f"Pipeline Version:     {report.pipeline_version}")
    print(f"Validation Date:      {report.validation_date.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Reference Dataset:    {report.reference_dataset}")
    print(f"Samples Evaluated:    {len(report.sample_results)}")
    print()
    print("OVERALL METRICS:")
    print(f"  Concordance:        {report.overall_concordance:.1%}")
    print(f"  Sensitivity:        {metrics.sensitivity:.1%}")
    print(f"  Specificity:        {metrics.specificity:.1%}")
    print(f"  Accuracy:           {metrics.accuracy:.1%}")
    print(f"  No-Call Rate:       {metrics.no_call_rate:.1%}")
    print()
    print("PER-GENE CONCORDANCE:")
    for gene, rate in sorted(report.per_gene_concordance.items()):
        status = "PASS" if rate >= 0.95 else "WARN" if rate >= 0.80 else "FAIL"
        print(f"  {gene:12s} {rate:6.1%}  [{status}]")

    print()

    # Print discordant calls
    discordant_count = sum(r.genes_discordant for r in report.sample_results)
    if discordant_count > 0:
        print("DISCORDANT CALLS:")
        for sample in report.sample_results:
            for gene in sample.gene_results:
                if gene.status.value == "discordant":
                    print(f"  {sample.sample_id} / {gene.gene}:")
                    print(f"    Expected: {gene.expected_diplotype}")
                    print(f"    Called:   {gene.called_diplotype}")
                    print(f"    Confidence: {gene.confidence_score:.2f}")
        print()

    # Write report
    write_validation_report(report, args.output)
    print(f"Validation report written to: {args.output}")

    # Exit code based on concordance
    if report.overall_concordance >= 0.95:
        print("\n[PASS] Validation passed (concordance >= 95%)")
        sys.exit(0)
    elif report.overall_concordance >= 0.80:
        print("\n[WARN] Validation passed with warnings (concordance 80-95%)")
        sys.exit(0)
    else:
        print("\n[FAIL] Validation failed (concordance < 80%)")
        sys.exit(1)


if __name__ == "__main__":
    main()
