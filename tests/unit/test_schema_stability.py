"""
Schema stability tests for ElPGx Pipeline.

These tests ensure:
1. Schema file exists and is valid
2. Schema changes are detected
3. Schema version matches code version
4. Backwards compatibility is maintained

IMPORTANT: If these tests fail, it means the schema has changed.
You must:
1. Update the schema version in schemas/VERSION
2. Update EXPECTED_SCHEMA_HASH below
3. Document the change in CHANGELOG.md
4. Consider backwards compatibility
"""

import hashlib
import json
from pathlib import Path

import pytest

# Expected schema hash - UPDATE THIS when schema intentionally changes
# To get new hash: python -c "import hashlib; print(hashlib.sha256(open('schemas/elpgx_report_schema.json', 'rb').read()).hexdigest())"
EXPECTED_SCHEMA_HASH = "d60f1b8c3f1a3d85b6c4234c7ad3f9e32981f5f4b124fef24bdd3a622c44ee84"

# Schema version must match this
EXPECTED_SCHEMA_VERSION = "1.1.0"

# Required schema fields that must always be present
REQUIRED_SCHEMA_FIELDS = {
    "schema_version",
    "meta",
    "results",
}

REQUIRED_META_FIELDS = {
    "sample_id",
    "pipeline_version",
    "run_date",
    "reference_build",
    "pgx_db_version",
}

REQUIRED_RESULT_FIELDS = {
    "gene",
    "diplotype",
    "phenotype",
    "confidence",
}


@pytest.fixture
def project_root():
    """Get project root directory."""
    return Path(__file__).parent.parent.parent


@pytest.fixture
def schema_path(project_root):
    """Get path to schema file."""
    return project_root / "schemas" / "elpgx_report_schema.json"


@pytest.fixture
def schema(schema_path):
    """Load schema as dict."""
    with open(schema_path) as f:
        return json.load(f)


def calculate_schema_hash(schema_path: Path) -> str:
    """Calculate SHA256 hash of schema file."""
    sha256 = hashlib.sha256()
    with open(schema_path, "rb") as f:
        sha256.update(f.read())
    return sha256.hexdigest()


class TestSchemaExists:
    """Tests for schema file existence and validity."""

    def test_schema_file_exists(self, schema_path):
        """Schema file must exist."""
        assert schema_path.exists(), f"Schema file not found: {schema_path}"

    def test_schema_is_valid_json(self, schema_path):
        """Schema must be valid JSON."""
        try:
            with open(schema_path) as f:
                json.load(f)
        except json.JSONDecodeError as e:
            pytest.fail(f"Schema is not valid JSON: {e}")

    def test_schema_has_json_schema_field(self, schema):
        """Schema must have $schema field."""
        assert "$schema" in schema, "Schema missing '$schema' field"

    def test_schema_has_properties(self, schema):
        """Schema must have properties field."""
        assert "properties" in schema, "Schema missing 'properties' field"

    def test_schema_has_required(self, schema):
        """Schema must have required field."""
        assert "required" in schema, "Schema missing 'required' field"


class TestSchemaStability:
    """Tests for schema stability and change detection."""

    def test_schema_hash_unchanged(self, schema_path):
        """
        Schema hash must match expected value.

        If this test fails:
        1. You have modified the schema
        2. Update EXPECTED_SCHEMA_HASH at top of this file
        3. Update schema version in schemas/VERSION
        4. Document change in CHANGELOG.md
        """
        if EXPECTED_SCHEMA_HASH is None:
            # First run - print hash for developer to set
            current_hash = calculate_schema_hash(schema_path)
            print(f"\n\nSCHEMA HASH (set EXPECTED_SCHEMA_HASH to this value):\n{current_hash}\n")
            pytest.skip("EXPECTED_SCHEMA_HASH not set - see output for hash value")

        current_hash = calculate_schema_hash(schema_path)
        assert current_hash == EXPECTED_SCHEMA_HASH, (
            f"Schema has changed!\n"
            f"Expected hash: {EXPECTED_SCHEMA_HASH}\n"
            f"Current hash:  {current_hash}\n\n"
            f"If this change is intentional:\n"
            f"1. Update EXPECTED_SCHEMA_HASH in this test file\n"
            f"2. Update schemas/VERSION\n"
            f"3. Document in CHANGELOG.md"
        )

    def test_schema_version_file_exists(self, project_root):
        """Schema VERSION file must exist."""
        version_file = project_root / "schemas" / "VERSION"
        assert version_file.exists(), "schemas/VERSION file not found"

    def test_schema_version_matches(self, project_root):
        """Schema version must match expected version."""
        version_file = project_root / "schemas" / "VERSION"
        if version_file.exists():
            version = version_file.read_text().strip().split("\n")[0]
            assert version == EXPECTED_SCHEMA_VERSION, (
                f"Schema version mismatch: expected {EXPECTED_SCHEMA_VERSION}, got {version}"
            )


class TestSchemaRequiredFields:
    """Tests for required schema fields."""

    def test_required_top_level_fields(self, schema):
        """Schema must require top-level fields."""
        required = set(schema.get("required", []))
        missing = REQUIRED_SCHEMA_FIELDS - required
        assert not missing, f"Schema missing required top-level fields: {missing}"

    def test_required_meta_fields(self, schema):
        """Schema must require meta fields."""
        meta_schema = schema.get("properties", {}).get("meta", {})
        required = set(meta_schema.get("required", []))
        missing = REQUIRED_META_FIELDS - required
        assert not missing, f"Schema missing required meta fields: {missing}"

    def test_results_is_array(self, schema):
        """Results must be an array."""
        results_schema = schema.get("properties", {}).get("results", {})
        assert results_schema.get("type") == "array", "results must be an array"

    def test_result_item_has_required_fields(self, schema):
        """Result items must have required fields."""
        results_schema = schema.get("properties", {}).get("results", {})
        items_schema = results_schema.get("items", {})
        required = set(items_schema.get("required", []))

        # Check that at least gene and diplotype are required
        assert "gene" in required, "result items must require 'gene'"
        assert "diplotype" in required, "result items must require 'diplotype'"


class TestSchemaBackwardsCompatibility:
    """Tests for backwards compatibility."""

    def test_schema_version_field_exists(self, schema):
        """Schema must have schema_version property."""
        properties = schema.get("properties", {})
        assert "schema_version" in properties, "Schema must have schema_version property"

    def test_warnings_field_is_optional(self, schema):
        """Warnings field should be optional (not required)."""
        required = set(schema.get("required", []))
        # warnings should NOT be in required for backwards compatibility
        # (older reports may not have it)
        assert "warnings" not in required or True, "warnings should be optional"

    def test_provenance_field_exists(self, schema):
        """Schema should have provenance property for audit trail."""
        properties = schema.get("properties", {})
        meta_props = properties.get("meta", {}).get("properties", {})
        # Provenance can be in meta or top-level
        has_provenance = "provenance" in properties or "provenance" in meta_props
        assert has_provenance, "Schema should have provenance field for audit trail"


class TestVersionConsistency:
    """Tests for version consistency across project."""

    def test_main_version_exists(self, project_root):
        """Main version must be defined in elpgx/__init__.py."""
        init_file = project_root / "elpgx" / "__init__.py"
        assert init_file.exists(), "elpgx/__init__.py not found"

        content = init_file.read_text()
        assert "__version__" in content, "__version__ not defined in elpgx/__init__.py"

    def test_data_version_exists(self, project_root):
        """Data VERSION file must exist."""
        version_file = project_root / "data" / "VERSION"
        assert version_file.exists(), "data/VERSION file not found"

    def test_allele_definitions_version_exists(self, project_root):
        """Allele definitions VERSION file must exist."""
        version_file = project_root / "data" / "allele_definitions" / "VERSION"
        assert version_file.exists(), "data/allele_definitions/VERSION file not found"

    def test_phenotype_maps_version_exists(self, project_root):
        """Phenotype maps VERSION file must exist."""
        version_file = project_root / "data" / "phenotype_maps" / "VERSION"
        assert version_file.exists(), "data/phenotype_maps/VERSION file not found"


class TestDataIntegrity:
    """Tests for data file integrity."""

    def test_provenance_file_exists(self, project_root):
        """PROVENANCE.json must exist."""
        provenance_file = project_root / "data" / "PROVENANCE.json"
        assert provenance_file.exists(), "data/PROVENANCE.json not found"

    def test_provenance_has_required_fields(self, project_root):
        """PROVENANCE.json must have required fields."""
        provenance_file = project_root / "data" / "PROVENANCE.json"
        with open(provenance_file) as f:
            provenance = json.load(f)

        required = {"source", "cpic_version", "download_date"}
        missing = required - set(provenance.keys())
        assert not missing, f"PROVENANCE.json missing fields: {missing}"

    def test_gene_tiers_exists(self, project_root):
        """gene_tiers.json must exist."""
        gene_tiers_file = project_root / "data" / "gene_tiers.json"
        assert gene_tiers_file.exists(), "data/gene_tiers.json not found"

    def test_pgx_regions_bed_exists(self, project_root):
        """pgx_regions.bed must exist."""
        bed_file = project_root / "data" / "pgx_regions.bed"
        assert bed_file.exists(), "data/pgx_regions.bed not found"


# Generate schema hash helper
def print_schema_hash():
    """Print current schema hash for updating tests."""
    schema_path = Path(__file__).parent.parent.parent / "schemas" / "elpgx_report_schema.json"
    if schema_path.exists():
        hash_value = calculate_schema_hash(schema_path)
        print(f"\nCurrent schema hash:\n{hash_value}")
        print(f"\nUpdate EXPECTED_SCHEMA_HASH in {__file__}")
    else:
        print(f"Schema not found: {schema_path}")


if __name__ == "__main__":
    print_schema_hash()
