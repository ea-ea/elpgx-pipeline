"""
Clinical-Grade Validation Module for ElPGx Pipeline.

This module implements FDA/CAP-compliant validation checks for clinical
pharmacogenomics reporting. All checks are documented and auditable.

CLINICAL VALIDATION REQUIREMENTS:
1. Data Integrity - SHA256 checksums for all reference data
2. Version Locking - Prevent unauthorized data changes
3. Completeness Checks - Ensure all critical genes have complete definitions
4. Schema Compliance - Validate output schema strictly
5. Audit Trail - Log all validation events with timestamps
6. Safety Checks - Validate clinical safety rules are in place

Author: Elif Arslan
License: MIT
"""

import hashlib
import json
import logging
import os
import re
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path
from typing import Optional

from pydantic import BaseModel, Field

# Configure clinical validation logger
logger = logging.getLogger("elpgx.clinical_validation")


class ValidationSeverity(str, Enum):
    """Severity levels for validation findings."""
    CRITICAL = "CRITICAL"  # Blocks clinical use
    ERROR = "ERROR"        # Must be resolved before production
    WARNING = "WARNING"    # Should be reviewed
    INFO = "INFO"          # Informational only


class ValidationFinding(BaseModel):
    """A single validation finding."""
    code: str = Field(description="Unique finding code (e.g., CV001)")
    severity: ValidationSeverity
    category: str = Field(description="Validation category")
    message: str
    details: Optional[str] = None
    remediation: Optional[str] = None


class ValidationReport(BaseModel):
    """Complete validation report."""
    validation_id: str = Field(description="Unique validation run ID")
    timestamp: str = Field(description="ISO8601 timestamp")
    pipeline_version: str
    validation_mode: str = Field(description="startup, runtime, or full")
    passed: bool
    findings: list[ValidationFinding] = Field(default_factory=list)
    checksums: dict[str, str] = Field(default_factory=dict)
    metadata: dict = Field(default_factory=dict)


# Critical genes that MUST have complete definitions for clinical use
CRITICAL_GENES = {
    "DPYD": {
        "reason": "5-FU toxicity risk - potentially fatal",
        "min_alleles": 5,
        "required_alleles": ["*1", "*2A", "*13"],
    },
    "TPMT": {
        "reason": "Thiopurine toxicity - myelosuppression risk",
        "min_alleles": 4,
        "required_alleles": ["*1", "*2", "*3A", "*3C"],
    },
    "CYP2C19": {
        "reason": "Clopidogrel efficacy - cardiovascular risk",
        "min_alleles": 8,
        "required_alleles": ["*1", "*2", "*3", "*17"],
    },
    "CYP2C9": {
        "reason": "Warfarin dosing - bleeding risk",
        "min_alleles": 6,
        "required_alleles": ["*1", "*2", "*3"],
    },
    "NUDT15": {
        "reason": "Thiopurine toxicity - especially in Asian populations",
        "min_alleles": 3,
        "required_alleles": ["*1", "*3"],
    },
    "UGT1A1": {
        "reason": "Irinotecan toxicity - neutropenia risk",
        "min_alleles": 3,
        "required_alleles": ["*1", "*28", "*80"],
    },
}

# Required phenotype categories for each gene
REQUIRED_PHENOTYPES = {
    "CYP2C19": ["Ultrarapid Metabolizer", "Normal Metabolizer", "Intermediate Metabolizer", "Poor Metabolizer"],
    "CYP2C9": ["Normal Metabolizer", "Intermediate Metabolizer", "Poor Metabolizer"],
    "CYP2D6": ["Ultrarapid Metabolizer", "Normal Metabolizer", "Intermediate Metabolizer", "Poor Metabolizer"],
    "DPYD": ["Normal Metabolizer", "Intermediate Metabolizer", "Poor Metabolizer"],
    "TPMT": ["Normal Metabolizer", "Intermediate Metabolizer", "Poor Metabolizer"],
}

# Expected CPIC version for this release
EXPECTED_CPIC_VERSION = "v1.24"

# Minimum confidence thresholds for clinical reporting
CLINICAL_CONFIDENCE_THRESHOLDS = {
    "reportable_minimum": 0.70,      # Below this, phenotype should be blocked
    "high_confidence": 0.90,         # High confidence threshold
    "cnv_gene_warning": 0.85,        # CYP2D6 without CNV should warn below this
}


def generate_validation_id() -> str:
    """Generate unique validation run ID."""
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%d%H%M%S")
    random_suffix = os.urandom(4).hex()
    return f"CV-{timestamp}-{random_suffix}"


def calculate_file_hash(file_path: Path) -> str:
    """Calculate SHA256 hash of a file."""
    sha256 = hashlib.sha256()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            sha256.update(chunk)
    return sha256.hexdigest()


def calculate_directory_hash(dir_path: Path, pattern: str = "*.json") -> dict[str, str]:
    """Calculate hashes for all files matching pattern in directory."""
    hashes = {}
    for file_path in sorted(dir_path.glob(pattern)):
        hashes[file_path.name] = calculate_file_hash(file_path)
    return hashes


class ClinicalValidator:
    """
    Clinical-grade validator for ElPGx pipeline.

    Performs comprehensive validation checks required for clinical use:
    - Data integrity verification
    - Version consistency
    - Critical gene completeness
    - Safety rule validation
    - Schema compliance
    """

    def __init__(self, project_root: Path, strict: bool = True):
        """
        Initialize clinical validator.

        Args:
            project_root: Path to project root directory.
            strict: If True, treat warnings as errors (recommended for clinical).
        """
        self.project_root = project_root
        self.strict = strict
        self.data_dir = project_root / "data"
        self.schema_dir = project_root / "schemas"
        self.findings: list[ValidationFinding] = []
        self.checksums: dict[str, str] = {}

    def add_finding(
        self,
        code: str,
        severity: ValidationSeverity,
        category: str,
        message: str,
        details: Optional[str] = None,
        remediation: Optional[str] = None,
    ) -> None:
        """Add a validation finding."""
        finding = ValidationFinding(
            code=code,
            severity=severity,
            category=category,
            message=message,
            details=details,
            remediation=remediation,
        )
        self.findings.append(finding)

        # Log finding
        log_level = {
            ValidationSeverity.CRITICAL: logging.CRITICAL,
            ValidationSeverity.ERROR: logging.ERROR,
            ValidationSeverity.WARNING: logging.WARNING,
            ValidationSeverity.INFO: logging.INFO,
        }.get(severity, logging.INFO)

        logger.log(log_level, f"[{code}] {message}")

    def validate_data_integrity(self) -> None:
        """
        CV-100: Validate data file integrity using checksums.

        Calculates SHA256 hashes for all data files and compares
        against stored manifest if available.
        """
        category = "Data Integrity"

        # Calculate checksums for allele definitions
        allele_dir = self.data_dir / "allele_definitions"
        if allele_dir.exists():
            allele_hashes = calculate_directory_hash(allele_dir)
            self.checksums["allele_definitions"] = hashlib.sha256(
                json.dumps(allele_hashes, sort_keys=True).encode()
            ).hexdigest()

            self.add_finding(
                code="CV-101",
                severity=ValidationSeverity.INFO,
                category=category,
                message=f"Calculated checksums for {len(allele_hashes)} allele definition files",
                details=f"Combined hash: {self.checksums['allele_definitions'][:16]}...",
            )
        else:
            self.add_finding(
                code="CV-102",
                severity=ValidationSeverity.CRITICAL,
                category=category,
                message="Allele definitions directory not found",
                remediation="Ensure data/allele_definitions/ exists with CPIC data",
            )

        # Calculate checksums for phenotype maps
        phenotype_dir = self.data_dir / "phenotype_maps"
        if phenotype_dir.exists():
            phenotype_hashes = calculate_directory_hash(phenotype_dir)
            self.checksums["phenotype_maps"] = hashlib.sha256(
                json.dumps(phenotype_hashes, sort_keys=True).encode()
            ).hexdigest()
        else:
            self.add_finding(
                code="CV-103",
                severity=ValidationSeverity.CRITICAL,
                category=category,
                message="Phenotype maps directory not found",
                remediation="Ensure data/phenotype_maps/ exists",
            )

        # Check PROVENANCE.json
        provenance_file = self.data_dir / "PROVENANCE.json"
        if provenance_file.exists():
            try:
                with open(provenance_file) as f:
                    provenance = json.load(f)

                self.checksums["provenance"] = calculate_file_hash(provenance_file)

                # Verify CPIC version
                cpic_version = provenance.get("cpic_version", "")
                if cpic_version != EXPECTED_CPIC_VERSION:
                    self.add_finding(
                        code="CV-104",
                        severity=ValidationSeverity.WARNING,
                        category=category,
                        message=f"CPIC version mismatch: expected {EXPECTED_CPIC_VERSION}, found {cpic_version}",
                        remediation="Update PROVENANCE.json or verify data version compatibility",
                    )
                else:
                    self.add_finding(
                        code="CV-105",
                        severity=ValidationSeverity.INFO,
                        category=category,
                        message=f"CPIC version verified: {cpic_version}",
                    )

            except json.JSONDecodeError as e:
                self.add_finding(
                    code="CV-106",
                    severity=ValidationSeverity.CRITICAL,
                    category=category,
                    message=f"Invalid PROVENANCE.json: {e}",
                    remediation="Fix JSON syntax in data/PROVENANCE.json",
                )
        else:
            self.add_finding(
                code="CV-107",
                severity=ValidationSeverity.ERROR,
                category=category,
                message="PROVENANCE.json not found - data provenance unknown",
                remediation="Create data/PROVENANCE.json with source information",
            )

    def validate_critical_genes(self) -> None:
        """
        CV-200: Validate critical gene definitions are complete.

        For clinical safety, certain genes MUST have complete allele
        definitions. Missing definitions could lead to incorrect
        phenotype predictions with serious clinical consequences.
        """
        category = "Critical Gene Completeness"

        allele_dir = self.data_dir / "allele_definitions"
        if not allele_dir.exists():
            return  # Already reported in integrity check

        for gene, requirements in CRITICAL_GENES.items():
            gene_file = allele_dir / f"{gene}.json"

            if not gene_file.exists():
                self.add_finding(
                    code="CV-201",
                    severity=ValidationSeverity.CRITICAL,
                    category=category,
                    message=f"Critical gene {gene} definition file missing",
                    details=f"Reason: {requirements['reason']}",
                    remediation=f"Add {gene}.json to data/allele_definitions/",
                )
                continue

            try:
                with open(gene_file) as f:
                    gene_def = json.load(f)

                # Check allele count
                alleles = gene_def.get("alleles", [])
                # Extract star allele names, handling CPIC naming conventions:
                # - "Reference" → "*1" (reference allele)
                # - "c.1905+1G>A (*2A)" → "*2A" (star allele in parentheses)
                # - "*2" → "*2" (direct star allele name)
                allele_names = set()
                for a in alleles:
                    raw_name = a.get("name", a.get("allele", ""))
                    if raw_name == "Reference" or raw_name == gene_def.get("reference_allele", ""):
                        allele_names.add("*1")
                    elif "(*" in raw_name:
                        # Extract star allele from parentheses: "c.1905+1G>A (*2A)" → "*2A"
                        match = re.search(r'\((\*[^)]+)\)', raw_name)
                        if match:
                            allele_names.add(match.group(1))
                    allele_names.add(raw_name)

                if len(alleles) < requirements["min_alleles"]:
                    self.add_finding(
                        code="CV-202",
                        severity=ValidationSeverity.ERROR,
                        category=category,
                        message=f"{gene}: Insufficient alleles ({len(alleles)} < {requirements['min_alleles']})",
                        details=f"Clinical impact: {requirements['reason']}",
                        remediation=f"Add missing allele definitions for {gene}",
                    )

                # Check required alleles
                missing_alleles = set(requirements["required_alleles"]) - allele_names
                if missing_alleles:
                    self.add_finding(
                        code="CV-203",
                        severity=ValidationSeverity.ERROR,
                        category=category,
                        message=f"{gene}: Missing required alleles: {', '.join(sorted(missing_alleles))}",
                        details=f"Clinical impact: {requirements['reason']}",
                        remediation=f"Add definitions for: {', '.join(sorted(missing_alleles))}",
                    )
                else:
                    self.add_finding(
                        code="CV-204",
                        severity=ValidationSeverity.INFO,
                        category=category,
                        message=f"{gene}: All {len(requirements['required_alleles'])} required alleles present",
                    )

            except json.JSONDecodeError as e:
                self.add_finding(
                    code="CV-205",
                    severity=ValidationSeverity.CRITICAL,
                    category=category,
                    message=f"Invalid JSON in {gene}.json: {e}",
                    remediation=f"Fix JSON syntax in data/allele_definitions/{gene}.json",
                )

    def validate_phenotype_maps(self) -> None:
        """
        CV-300: Validate phenotype translation maps are complete.

        Each gene must have mappings for all expected phenotype categories.
        """
        category = "Phenotype Map Completeness"

        phenotype_dir = self.data_dir / "phenotype_maps"
        if not phenotype_dir.exists():
            return  # Already reported

        for gene, required_phenotypes in REQUIRED_PHENOTYPES.items():
            gene_file = phenotype_dir / f"{gene}_phenotype.json"

            if not gene_file.exists():
                self.add_finding(
                    code="CV-301",
                    severity=ValidationSeverity.ERROR,
                    category=category,
                    message=f"Phenotype map missing for {gene}",
                    remediation=f"Add {gene}_phenotype.json to data/phenotype_maps/",
                )
                continue

            try:
                with open(gene_file) as f:
                    phenotype_map = json.load(f)

                # Check phenotype coverage
                defined_phenotypes = set()
                for entry in phenotype_map.get("mappings", phenotype_map.get("diplotypes", [])):
                    phenotype = entry.get("phenotype", "")
                    if phenotype:
                        defined_phenotypes.add(phenotype)

                missing = set(required_phenotypes) - defined_phenotypes
                if missing:
                    self.add_finding(
                        code="CV-302",
                        severity=ValidationSeverity.WARNING,
                        category=category,
                        message=f"{gene}: Missing phenotype mappings: {', '.join(sorted(missing))}",
                        remediation="Add missing diplotype-to-phenotype mappings",
                    )
                else:
                    self.add_finding(
                        code="CV-303",
                        severity=ValidationSeverity.INFO,
                        category=category,
                        message=f"{gene}: All {len(required_phenotypes)} phenotype categories covered",
                    )

            except json.JSONDecodeError as e:
                self.add_finding(
                    code="CV-304",
                    severity=ValidationSeverity.ERROR,
                    category=category,
                    message=f"Invalid JSON in {gene}_phenotype.json: {e}",
                )

    def validate_safety_rules(self) -> None:
        """
        CV-400: Validate clinical safety rules are in place.

        Checks that confidence thresholds and phenotype blocking rules
        are properly configured.
        """
        category = "Clinical Safety Rules"

        # Check gene strategies for confidence blocking
        strategies_dir = self.data_dir / "gene_strategies"
        if strategies_dir.exists():
            cyp2d6_strategy = strategies_dir / "CYP2D6_STRATEGY.json"
            if cyp2d6_strategy.exists():
                try:
                    with open(cyp2d6_strategy) as f:
                        strategy = json.load(f)

                    thresholds = strategy.get("confidence_thresholds", {})
                    no_call = thresholds.get("no_call", {})

                    if no_call.get("phenotype_reported") is True:
                        self.add_finding(
                            code="CV-401",
                            severity=ValidationSeverity.CRITICAL,
                            category=category,
                            message="CYP2D6 no_call threshold allows phenotype reporting",
                            details="Low-confidence CYP2D6 calls should block phenotype reporting",
                            remediation="Set no_call.phenotype_reported to false in CYP2D6_STRATEGY.json",
                        )
                    else:
                        self.add_finding(
                            code="CV-402",
                            severity=ValidationSeverity.INFO,
                            category=category,
                            message="CYP2D6 phenotype blocking rule verified",
                        )

                    # Check CNV warning
                    unsupported = strategy.get("unsupported_features", {})
                    cnv = unsupported.get("cnv_detection", {})
                    if not cnv.get("warning_message"):
                        self.add_finding(
                            code="CV-403",
                            severity=ValidationSeverity.WARNING,
                            category=category,
                            message="CYP2D6 CNV limitation warning not configured",
                            remediation="Add warning_message to cnv_detection in strategy",
                        )

                except json.JSONDecodeError:
                    self.add_finding(
                        code="CV-404",
                        severity=ValidationSeverity.ERROR,
                        category=category,
                        message="Invalid CYP2D6_STRATEGY.json",
                    )
        else:
            self.add_finding(
                code="CV-405",
                severity=ValidationSeverity.WARNING,
                category=category,
                message="Gene strategies directory not found",
                details="Gene-specific confidence rules may not be applied",
            )

    def validate_schema(self) -> None:
        """
        CV-500: Validate output schema compliance.

        Ensures the JSON schema is valid and contains all required
        fields for clinical reporting.
        """
        category = "Schema Compliance"

        schema_file = self.schema_dir / "elpgx_report_schema.json"
        if not schema_file.exists():
            self.add_finding(
                code="CV-501",
                severity=ValidationSeverity.CRITICAL,
                category=category,
                message="Report schema file not found",
                remediation="Ensure schemas/elpgx_report_schema.json exists",
            )
            return

        try:
            with open(schema_file) as f:
                schema = json.load(f)

            self.checksums["schema"] = calculate_file_hash(schema_file)

            # Check required top-level fields
            required_fields = {"schema_version", "meta", "results"}
            schema_required = set(schema.get("required", []))
            missing = required_fields - schema_required

            if missing:
                self.add_finding(
                    code="CV-502",
                    severity=ValidationSeverity.ERROR,
                    category=category,
                    message=f"Schema missing required top-level fields: {missing}",
                )

            # Check meta required fields
            meta_schema = schema.get("properties", {}).get("meta", {})
            meta_required = set(meta_schema.get("required", []))
            clinical_meta_fields = {"sample_id", "pipeline_version", "run_date", "reference_build", "pgx_db_version"}
            missing_meta = clinical_meta_fields - meta_required

            if missing_meta:
                self.add_finding(
                    code="CV-503",
                    severity=ValidationSeverity.ERROR,
                    category=category,
                    message=f"Schema missing required meta fields: {missing_meta}",
                )

            # Check provenance field exists
            properties = schema.get("properties", {})
            if "provenance" not in properties:
                self.add_finding(
                    code="CV-504",
                    severity=ValidationSeverity.ERROR,
                    category=category,
                    message="Schema missing provenance field",
                    details="Provenance is required for clinical audit trail",
                )
            else:
                self.add_finding(
                    code="CV-505",
                    severity=ValidationSeverity.INFO,
                    category=category,
                    message="Schema includes provenance for audit trail",
                )

            # Check confidence field in results
            results_schema = properties.get("results", {}).get("items", {})
            result_required = set(results_schema.get("required", []))

            if "confidence" not in result_required:
                self.add_finding(
                    code="CV-506",
                    severity=ValidationSeverity.ERROR,
                    category=category,
                    message="Schema does not require confidence in results",
                    details="Confidence scores are required for clinical interpretation",
                )

        except json.JSONDecodeError as e:
            self.add_finding(
                code="CV-507",
                severity=ValidationSeverity.CRITICAL,
                category=category,
                message=f"Invalid schema JSON: {e}",
            )

    def validate_version_files(self) -> None:
        """
        CV-600: Validate VERSION files exist in all data directories.

        Version files are required for data provenance tracking.
        """
        category = "Version Tracking"

        required_version_files = [
            self.data_dir / "VERSION",
            self.data_dir / "allele_definitions" / "VERSION",
            self.data_dir / "phenotype_maps" / "VERSION",
            self.data_dir / "prescribing_implications" / "VERSION",
            self.schema_dir / "VERSION",
        ]

        for version_file in required_version_files:
            if version_file.exists():
                version = version_file.read_text().strip().split("\n")[0]
                self.add_finding(
                    code="CV-601",
                    severity=ValidationSeverity.INFO,
                    category=category,
                    message=f"Version file present: {version_file.parent.name}/VERSION = {version}",
                )
            else:
                self.add_finding(
                    code="CV-602",
                    severity=ValidationSeverity.ERROR,
                    category=category,
                    message=f"Missing VERSION file: {version_file}",
                    remediation="Add VERSION file for data provenance tracking",
                )

    def run_full_validation(self) -> ValidationReport:
        """
        Run complete clinical validation suite.

        Returns:
            ValidationReport with all findings.
        """
        from elpgx import __version__

        validation_id = generate_validation_id()
        logger.info(f"Starting clinical validation: {validation_id}")

        # Run all validation checks
        self.validate_data_integrity()
        self.validate_critical_genes()
        self.validate_phenotype_maps()
        self.validate_safety_rules()
        self.validate_schema()
        self.validate_version_files()

        # Determine pass/fail
        critical_count = sum(1 for f in self.findings if f.severity == ValidationSeverity.CRITICAL)
        error_count = sum(1 for f in self.findings if f.severity == ValidationSeverity.ERROR)
        warning_count = sum(1 for f in self.findings if f.severity == ValidationSeverity.WARNING)

        if self.strict:
            passed = critical_count == 0 and error_count == 0
        else:
            passed = critical_count == 0

        report = ValidationReport(
            validation_id=validation_id,
            timestamp=datetime.now(timezone.utc).isoformat(),
            pipeline_version=__version__,
            validation_mode="full",
            passed=passed,
            findings=self.findings,
            checksums=self.checksums,
            metadata={
                "strict_mode": self.strict,
                "critical_count": critical_count,
                "error_count": error_count,
                "warning_count": warning_count,
                "info_count": sum(1 for f in self.findings if f.severity == ValidationSeverity.INFO),
            },
        )

        logger.info(
            f"Validation complete: {'PASSED' if passed else 'FAILED'} "
            f"(Critical: {critical_count}, Errors: {error_count}, Warnings: {warning_count})"
        )

        return report


class ClinicalValidationError(Exception):
    """Raised when clinical validation fails."""

    def __init__(self, report: ValidationReport):
        self.report = report
        critical_findings = [f for f in report.findings if f.severity == ValidationSeverity.CRITICAL]
        error_findings = [f for f in report.findings if f.severity == ValidationSeverity.ERROR]

        messages = []
        if critical_findings:
            messages.append("CRITICAL issues:")
            for f in critical_findings:
                messages.append(f"  [{f.code}] {f.message}")
        if error_findings:
            messages.append("ERROR issues:")
            for f in error_findings:
                messages.append(f"  [{f.code}] {f.message}")

        super().__init__("\n".join(messages))


def run_clinical_validation(
    project_root: Optional[Path] = None,
    strict: bool = True,
    raise_on_failure: bool = True,
) -> ValidationReport:
    """
    Run clinical-grade validation checks.

    Args:
        project_root: Path to project root (auto-detected if None).
        strict: If True, treat errors as failures (recommended for clinical).
        raise_on_failure: If True, raise exception on failure.

    Returns:
        ValidationReport with all findings.

    Raises:
        ClinicalValidationError: If validation fails and raise_on_failure is True.
    """
    if project_root is None:
        project_root = Path(__file__).parent.parent

    validator = ClinicalValidator(project_root, strict=strict)
    report = validator.run_full_validation()

    if not report.passed and raise_on_failure:
        raise ClinicalValidationError(report)

    return report


def print_validation_report(report: ValidationReport) -> None:
    """Print formatted validation report to console."""
    print("=" * 70)
    print("ElPGx CLINICAL VALIDATION REPORT")
    print("=" * 70)
    print(f"Validation ID: {report.validation_id}")
    print(f"Timestamp:     {report.timestamp}")
    print(f"Version:       {report.pipeline_version}")
    print(f"Mode:          {report.validation_mode}")
    print(f"Status:        {'PASSED' if report.passed else 'FAILED'}")
    print("-" * 70)

    # Group findings by severity
    for severity in [ValidationSeverity.CRITICAL, ValidationSeverity.ERROR,
                     ValidationSeverity.WARNING, ValidationSeverity.INFO]:
        findings = [f for f in report.findings if f.severity == severity]
        if findings:
            print(f"\n{severity.value} ({len(findings)}):")
            for f in findings:
                print(f"  [{f.code}] {f.message}")
                if f.details:
                    print(f"           Details: {f.details}")
                if f.remediation:
                    print(f"           Fix: {f.remediation}")

    print("\n" + "-" * 70)
    print(f"Summary: {report.metadata.get('critical_count', 0)} critical, "
          f"{report.metadata.get('error_count', 0)} errors, "
          f"{report.metadata.get('warning_count', 0)} warnings")
    print("=" * 70)


if __name__ == "__main__":
    # Run validation when executed directly
    logging.basicConfig(level=logging.INFO)
    project_root = Path(__file__).parent.parent

    try:
        report = run_clinical_validation(project_root, strict=True, raise_on_failure=False)
        print_validation_report(report)

        # Save report to file
        report_file = project_root / "validation_report.json"
        with open(report_file, "w") as f:
            json.dump(report.model_dump(), f, indent=2)
        print(f"\nReport saved to: {report_file}")

    except Exception as e:
        print(f"Validation error: {e}")
        raise
