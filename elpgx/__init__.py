"""
ElPGx - Clinical Pharmacogenetics Analysis Pipeline

A comprehensive pharmacogenetic analysis toolkit for clinical NGS platforms.
Processes WGS/WES data to predict drug response phenotypes based on genetic variants.

This package provides:
- Input validation for BAM/CRAM files (hg38/GRCh38)
- VCF parsing with dual backend support (cyvcf2/pysam)
- Star-allele diplotype calling
- CPIC-based phenotype translation
- Confidence scoring with quality metrics
- Clinical report generation with JSON schema validation

Author: Elif Arslan
License: MIT
"""

__version__ = "1.1.0"
__author__ = "Elif Arslan"
__package_name__ = "elpgx"

from elpgx.models import (
    BiomarkerType,
    ConfidenceFlag,
    ConfidenceMetrics,
    GeneResult,
    GeneTier,
    PgxReport,
    PipelineMeta,
    Variant,
)

__all__ = [
    "BiomarkerType",
    "ConfidenceFlag",
    "ConfidenceMetrics",
    "GeneResult",
    "GeneTier",
    "PgxReport",
    "PipelineMeta",
    "Variant",
    "__version__",
    "__author__",
]
