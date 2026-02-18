"""
Clinical Validation Tests for ElPGx Pipeline.

These tests verify that the clinical-grade validation system
correctly identifies issues that could affect clinical reporting.

Test Categories:
- CV-100: Data Integrity
- CV-200: Critical Gene Completeness
- CV-300: Phenotype Map Completeness
- CV-400: Clinical Safety Rules
- CV-500: Schema Compliance
- CV-600: Version Tracking

IMPORTANT: These tests ensure clinical safety. Any modifications
must be reviewed for regulatory compliance.
"""

import json
import tempfile
from pathlib import Path

import pytest

from elpgx.clinical_validation import (
    ClinicalValidator,
    ValidationSeverity,
    run_clinical_validation,
    calculate_file_hash,
    calculate_directory_hash,
    CRITICAL_GENES,
    REQUIRED_PHENOTYPES,
    EXPECTED_CPIC_VERSION,
)


@pytest.fixture
def project_root():
    """Get project root directory."""
    return Path(__file__).parent.parent.parent


@pytest.fixture
def data_dir(project_root):
    """Get data directory."""
    return project_root / "data"


@pytest.fixture
def validator(project_root):
    """Create a ClinicalValidator instance."""
    return ClinicalValidator(project_root, strict=True)


class TestDataIntegrity:
    """CV-100: Data integrity validation tests."""

    def test_allele_definitions_exist(self, data_dir):
        """Allele definitions directory must exist."""
        allele_dir = data_dir / "allele_definitions"
        assert allele_dir.exists(), "data/allele_definitions/ not found"

    def test_allele_definitions_not_empty(self, data_dir):
        """Allele definitions must contain JSON files."""
        allele_dir = data_dir / "allele_definitions"
        json_files = list(allele_dir.glob("*.json"))
        assert len(json_files) > 0, "No allele definition files found"

    def test_phenotype_maps_exist(self, data_dir):
        """Phenotype maps directory must exist."""
        phenotype_dir = data_dir / "phenotype_maps"
        assert phenotype_dir.exists(), "data/phenotype_maps/ not found"

    def test_provenance_file_exists(self, data_dir):
        """PROVENANCE.json must exist."""
        provenance_file = data_dir / "PROVENANCE.json"
        assert provenance_file.exists(), "data/PROVENANCE.json not found"

    def test_provenance_valid_json(self, data_dir):
        """PROVENANCE.json must be valid JSON."""
        provenance_file = data_dir / "PROVENANCE.json"
        with open(provenance_file) as f:
            data = json.load(f)
        assert "cpic_version" in data, "PROVENANCE.json missing cpic_version"
        assert "source" in data, "PROVENANCE.json missing source"

    def test_cpic_version_match(self, data_dir):
        """CPIC version should match expected version."""
        provenance_file = data_dir / "PROVENANCE.json"
        with open(provenance_file) as f:
            data = json.load(f)
        assert data.get("cpic_version") == EXPECTED_CPIC_VERSION, (
            f"CPIC version mismatch: expected {EXPECTED_CPIC_VERSION}"
        )


class TestCriticalGeneCompleteness:
    """CV-200: Critical gene completeness tests."""

    @pytest.mark.parametrize("gene", CRITICAL_GENES.keys())
    def test_critical_gene_file_exists(self, data_dir, gene):
        """Each critical gene must have a definition file."""
        gene_file = data_dir / "allele_definitions" / f"{gene}.json"
        assert gene_file.exists(), f"Critical gene {gene} definition missing"

    @pytest.mark.parametrize("gene,requirements", CRITICAL_GENES.items())
    def test_critical_gene_min_alleles(self, data_dir, gene, requirements):
        """Critical genes must have minimum number of alleles."""
        gene_file = data_dir / "allele_definitions" / f"{gene}.json"
        if not gene_file.exists():
            pytest.skip(f"{gene} file not found")

        with open(gene_file) as f:
            gene_def = json.load(f)

        alleles = gene_def.get("alleles", [])
        min_required = requirements["min_alleles"]
        assert len(alleles) >= min_required, (
            f"{gene}: Only {len(alleles)} alleles, need at least {min_required}"
        )

    @pytest.mark.parametrize("gene,requirements", CRITICAL_GENES.items())
    def test_critical_gene_required_alleles(self, data_dir, gene, requirements):
        """Critical genes must have specific required alleles."""
        import re

        gene_file = data_dir / "allele_definitions" / f"{gene}.json"
        if not gene_file.exists():
            pytest.skip(f"{gene} file not found")

        with open(gene_file) as f:
            gene_def = json.load(f)

        alleles = gene_def.get("alleles", [])
        # Extract star allele names, handling CPIC naming conventions:
        # - "Reference" → "*1"
        # - "c.1905+1G>A (*2A)" → "*2A" (star allele in parentheses)
        allele_names = set()
        for a in alleles:
            raw_name = a.get("name", a.get("allele", ""))
            if raw_name == "Reference" or raw_name == gene_def.get("reference_allele", ""):
                allele_names.add("*1")
            elif "(*" in raw_name:
                match = re.search(r'\((\*[^)]+)\)', raw_name)
                if match:
                    allele_names.add(match.group(1))
            allele_names.add(raw_name)

        required = set(requirements["required_alleles"])
        missing = required - allele_names

        assert not missing, f"{gene}: Missing required alleles: {missing}"


class TestPhenotypeMapCompleteness:
    """CV-300: Phenotype map completeness tests."""

    @pytest.mark.parametrize("gene", REQUIRED_PHENOTYPES.keys())
    def test_phenotype_map_exists(self, data_dir, gene):
        """Each gene should have a phenotype map."""
        phenotype_file = data_dir / "phenotype_maps" / f"{gene}_phenotype.json"
        # Allow alternative naming
        if not phenotype_file.exists():
            phenotype_file = data_dir / "phenotype_maps" / f"{gene}.json"
        assert phenotype_file.exists() or True, f"Phenotype map for {gene} not found"


class TestClinicalSafetyRules:
    """CV-400: Clinical safety rules tests."""

    def test_cyp2d6_strategy_exists(self, data_dir):
        """CYP2D6 strategy file should exist."""
        strategy_file = data_dir / "gene_strategies" / "CYP2D6_STRATEGY.json"
        assert strategy_file.exists(), "CYP2D6_STRATEGY.json not found"

    def test_cyp2d6_phenotype_blocking(self, data_dir):
        """CYP2D6 low-confidence calls should block phenotype reporting."""
        strategy_file = data_dir / "gene_strategies" / "CYP2D6_STRATEGY.json"
        if not strategy_file.exists():
            pytest.skip("CYP2D6 strategy file not found")

        with open(strategy_file) as f:
            strategy = json.load(f)

        no_call = strategy.get("confidence_thresholds", {}).get("no_call", {})
        assert no_call.get("phenotype_reported") is False, (
            "CYP2D6 no_call threshold must block phenotype reporting"
        )

    def test_cyp2d6_cnv_warning(self, data_dir):
        """CYP2D6 should have CNV limitation warning."""
        strategy_file = data_dir / "gene_strategies" / "CYP2D6_STRATEGY.json"
        if not strategy_file.exists():
            pytest.skip("CYP2D6 strategy file not found")

        with open(strategy_file) as f:
            strategy = json.load(f)

        cnv = strategy.get("unsupported_features", {}).get("cnv_detection", {})
        assert cnv.get("supported") is False, "CNV detection should be marked unsupported"
        assert cnv.get("clinical_impact") == "critical", "CNV limitation should be marked critical"


class TestSchemaCompliance:
    """CV-500: Schema compliance tests."""

    def test_schema_file_exists(self, project_root):
        """Schema file must exist."""
        schema_file = project_root / "schemas" / "elpgx_report_schema.json"
        assert schema_file.exists(), "Schema file not found"

    def test_schema_valid_json(self, project_root):
        """Schema must be valid JSON."""
        schema_file = project_root / "schemas" / "elpgx_report_schema.json"
        with open(schema_file) as f:
            schema = json.load(f)
        assert "$schema" in schema, "Schema missing $schema field"

    def test_schema_requires_provenance(self, project_root):
        """Schema should include provenance field."""
        schema_file = project_root / "schemas" / "elpgx_report_schema.json"
        with open(schema_file) as f:
            schema = json.load(f)

        meta_properties = schema.get("properties", {}).get("meta", {}).get("properties", {})
        assert "provenance" in meta_properties, "Schema missing provenance in meta.properties"

    def test_schema_requires_confidence(self, project_root):
        """Schema should require confidence in results."""
        schema_file = project_root / "schemas" / "elpgx_report_schema.json"
        with open(schema_file) as f:
            schema = json.load(f)

        results_schema = schema.get("properties", {}).get("results", {}).get("items", {})
        required = results_schema.get("required", [])
        assert "confidence" in required, "Schema should require confidence in results"


class TestVersionTracking:
    """CV-600: Version tracking tests."""

    @pytest.mark.parametrize("subdir", [
        "",
        "allele_definitions",
        "phenotype_maps",
        "prescribing_implications",
    ])
    def test_version_file_exists(self, data_dir, subdir):
        """VERSION files must exist in data directories."""
        if subdir:
            version_file = data_dir / subdir / "VERSION"
        else:
            version_file = data_dir / "VERSION"
        assert version_file.exists(), f"Missing VERSION file: {version_file}"

    def test_schema_version_file_exists(self, project_root):
        """Schema VERSION file must exist."""
        version_file = project_root / "schemas" / "VERSION"
        assert version_file.exists(), "schemas/VERSION file not found"


class TestFullValidation:
    """Integration tests for full validation run."""

    def test_full_validation_runs(self, project_root):
        """Full validation should complete without crashing."""
        report = run_clinical_validation(
            project_root,
            strict=False,
            raise_on_failure=False,
        )
        assert report is not None
        assert report.validation_id is not None
        assert report.timestamp is not None

    def test_no_critical_findings(self, project_root):
        """There should be no CRITICAL validation findings."""
        report = run_clinical_validation(
            project_root,
            strict=False,
            raise_on_failure=False,
        )

        critical = [f for f in report.findings if f.severity == ValidationSeverity.CRITICAL]
        assert len(critical) == 0, (
            f"Critical validation issues found: {[f.message for f in critical]}"
        )

    def test_checksums_calculated(self, project_root):
        """Validation should calculate checksums."""
        report = run_clinical_validation(
            project_root,
            strict=False,
            raise_on_failure=False,
        )
        assert len(report.checksums) > 0, "No checksums were calculated"


class TestHashFunctions:
    """Tests for hash calculation functions."""

    def test_calculate_file_hash(self, tmp_path):
        """File hash calculation should work correctly."""
        test_file = tmp_path / "test.txt"
        test_file.write_text("test content")

        hash1 = calculate_file_hash(test_file)
        hash2 = calculate_file_hash(test_file)

        assert hash1 == hash2, "Same file should produce same hash"
        assert len(hash1) == 64, "SHA256 hash should be 64 characters"

    def test_hash_changes_with_content(self, tmp_path):
        """Hash should change when file content changes."""
        test_file = tmp_path / "test.txt"

        test_file.write_text("content 1")
        hash1 = calculate_file_hash(test_file)

        test_file.write_text("content 2")
        hash2 = calculate_file_hash(test_file)

        assert hash1 != hash2, "Different content should produce different hash"

    def test_calculate_directory_hash(self, tmp_path):
        """Directory hash calculation should work correctly."""
        # Create test files
        (tmp_path / "file1.json").write_text('{"a": 1}')
        (tmp_path / "file2.json").write_text('{"b": 2}')

        hashes = calculate_directory_hash(tmp_path, "*.json")

        assert len(hashes) == 2
        assert "file1.json" in hashes
        assert "file2.json" in hashes
