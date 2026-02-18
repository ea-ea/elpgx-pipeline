"""
Startup validation module for ElPGx Pipeline.

Performs critical validation checks on pipeline initialization:
- Schema validation
- Data directory integrity
- Version consistency
- Required file presence

This module should be called at CLI/GUI startup and pipeline initialization.
"""

import hashlib
import json
import sys
from pathlib import Path
from typing import Optional

# Store schema hash for change detection
SCHEMA_HASH_FILE = ".schema_hash"


class StartupValidationError(Exception):
    """Raised when startup validation fails."""
    pass


def calculate_file_hash(file_path: Path) -> str:
    """Calculate SHA256 hash of a file."""
    sha256 = hashlib.sha256()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            sha256.update(chunk)
    return sha256.hexdigest()


def calculate_schema_hash(schema_path: Path) -> str:
    """Calculate hash of schema file for change detection."""
    return calculate_file_hash(schema_path)


def get_stored_schema_hash(project_root: Path) -> Optional[str]:
    """Get previously stored schema hash."""
    hash_file = project_root / SCHEMA_HASH_FILE
    if hash_file.exists():
        return hash_file.read_text().strip()
    return None


def store_schema_hash(project_root: Path, hash_value: str) -> None:
    """Store schema hash for future comparison."""
    hash_file = project_root / SCHEMA_HASH_FILE
    hash_file.write_text(hash_value)


def validate_schema_exists(schema_path: Path) -> list[str]:
    """Validate that schema file exists and is valid JSON."""
    errors = []

    if not schema_path.exists():
        errors.append(f"Schema file not found: {schema_path}")
        return errors

    try:
        with open(schema_path) as f:
            schema = json.load(f)

        # Check required schema fields
        if "$schema" not in schema:
            errors.append("Schema missing '$schema' field")
        if "properties" not in schema:
            errors.append("Schema missing 'properties' field")
        if "required" not in schema:
            errors.append("Schema missing 'required' field")

    except json.JSONDecodeError as e:
        errors.append(f"Invalid JSON in schema: {e}")

    return errors


def validate_data_versions(data_dir: Path) -> list[str]:
    """Validate VERSION files exist in data directories."""
    errors = []

    required_version_files = [
        data_dir / "VERSION",
        data_dir / "allele_definitions" / "VERSION",
        data_dir / "phenotype_maps" / "VERSION",
        data_dir / "prescribing_implications" / "VERSION",
    ]

    for version_file in required_version_files:
        if not version_file.exists():
            errors.append(f"Missing VERSION file: {version_file}")

    return errors


def validate_required_files(data_dir: Path) -> list[str]:
    """Validate required data files exist."""
    errors = []

    required_files = [
        data_dir / "gene_tiers.json",
        data_dir / "pgx_regions.bed",
        data_dir / "PROVENANCE.json",
    ]

    for required_file in required_files:
        if not required_file.exists():
            errors.append(f"Missing required file: {required_file}")

    # Check for at least some allele definitions
    allele_defs = list((data_dir / "allele_definitions").glob("*.json"))
    if len(allele_defs) == 0:
        errors.append("No allele definition files found in data/allele_definitions/")

    # Check for at least some phenotype maps
    phenotype_maps = list((data_dir / "phenotype_maps").glob("*.json"))
    if len(phenotype_maps) == 0:
        errors.append("No phenotype map files found in data/phenotype_maps/")

    return errors


def validate_version_consistency(project_root: Path) -> list[str]:
    """Validate version consistency across project."""
    errors = []
    warnings = []

    # Get main version from elpgx/__init__.py
    init_file = project_root / "elpgx" / "__init__.py"
    main_version = None

    if init_file.exists():
        content = init_file.read_text()
        for line in content.split("\n"):
            if "__version__" in line and "=" in line:
                main_version = line.split("=")[1].strip().strip('"\'')
                break

    if not main_version:
        errors.append("Could not determine main version from elpgx/__init__.py")
        return errors

    # Check data VERSION file
    data_version_file = project_root / "data" / "VERSION"
    if data_version_file.exists():
        data_version = data_version_file.read_text().strip().split("\n")[0]
        if data_version != main_version:
            warnings.append(
                f"Data version ({data_version}) differs from main version ({main_version})"
            )

    # Check schema VERSION file
    schema_version_file = project_root / "schemas" / "VERSION"
    if schema_version_file.exists():
        schema_version = schema_version_file.read_text().strip().split("\n")[0]
        if schema_version != main_version:
            warnings.append(
                f"Schema version ({schema_version}) differs from main version ({main_version})"
            )

    return errors  # Warnings are not failures


def validate_schema_unchanged(project_root: Path, schema_path: Path) -> tuple[bool, str]:
    """
    Check if schema has changed since last validation.

    Returns:
        Tuple of (changed: bool, message: str)
    """
    current_hash = calculate_schema_hash(schema_path)
    stored_hash = get_stored_schema_hash(project_root)

    if stored_hash is None:
        # First run - store the hash
        store_schema_hash(project_root, current_hash)
        return False, "Schema hash stored for first time"

    if current_hash != stored_hash:
        return True, f"Schema changed! Old: {stored_hash[:16]}... New: {current_hash[:16]}..."

    return False, "Schema unchanged"


def run_startup_validation(
    project_root: Path,
    strict: bool = False,
    update_hash: bool = False,
) -> tuple[bool, list[str], list[str]]:
    """
    Run all startup validation checks.

    Args:
        project_root: Path to project root directory.
        strict: If True, treat warnings as errors.
        update_hash: If True, update stored schema hash.

    Returns:
        Tuple of (success: bool, errors: list, warnings: list)
    """
    errors = []
    warnings = []

    data_dir = project_root / "data"
    schema_path = project_root / "schemas" / "elpgx_report_schema.json"

    # 1. Validate schema exists and is valid
    schema_errors = validate_schema_exists(schema_path)
    errors.extend(schema_errors)

    # 2. Validate data VERSION files
    version_errors = validate_data_versions(data_dir)
    errors.extend(version_errors)

    # 3. Validate required files
    file_errors = validate_required_files(data_dir)
    errors.extend(file_errors)

    # 4. Validate version consistency
    consistency_errors = validate_version_consistency(project_root)
    errors.extend(consistency_errors)

    # 5. Check schema changes
    if schema_path.exists():
        changed, message = validate_schema_unchanged(project_root, schema_path)
        if changed:
            if update_hash:
                store_schema_hash(project_root, calculate_schema_hash(schema_path))
                warnings.append(f"Schema changed - hash updated: {message}")
            else:
                errors.append(f"Schema change detected without version update: {message}")

    success = len(errors) == 0
    if strict:
        success = success and len(warnings) == 0

    return success, errors, warnings


def validate_on_startup(project_root: Optional[Path] = None) -> None:
    """
    Validate pipeline on startup. Raises exception on failure.

    This should be called at the start of CLI or GUI.
    """
    if project_root is None:
        # Try to find project root
        project_root = Path(__file__).parent.parent

    success, errors, warnings = run_startup_validation(project_root)

    # Print warnings
    for warning in warnings:
        print(f"[WARNING] {warning}", file=sys.stderr)

    if not success:
        error_msg = "Startup validation failed:\n" + "\n".join(f"  - {e}" for e in errors)
        raise StartupValidationError(error_msg)


def print_validation_status(project_root: Path) -> None:
    """Print validation status for debugging."""
    success, errors, warnings = run_startup_validation(project_root)

    print("=" * 60)
    print("ElPGx Startup Validation")
    print("=" * 60)

    if success:
        print("[PASS] All validation checks passed")
    else:
        print("[FAIL] Validation failed")

    if warnings:
        print("\nWarnings:")
        for w in warnings:
            print(f"  - {w}")

    if errors:
        print("\nErrors:")
        for e in errors:
            print(f"  - {e}")

    print("=" * 60)


if __name__ == "__main__":
    # Run validation when executed directly
    project_root = Path(__file__).parent.parent
    print_validation_status(project_root)
