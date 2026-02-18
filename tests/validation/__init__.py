"""
Validation module for ElPGx Pipeline.

Provides tools for testing against reference samples and
calculating concordance metrics.
"""

from .validation_framework import (
    ReferenceSample,
    ValidationReport,
    SampleValidationResult,
    GeneValidationResult,
    ValidationStatus,
    PerformanceMetrics,
    run_validation,
    load_reference_samples,
    calculate_metrics,
    write_validation_report,
)

__all__ = [
    "ReferenceSample",
    "ValidationReport",
    "SampleValidationResult",
    "GeneValidationResult",
    "ValidationStatus",
    "PerformanceMetrics",
    "run_validation",
    "load_reference_samples",
    "calculate_metrics",
    "write_validation_report",
]
