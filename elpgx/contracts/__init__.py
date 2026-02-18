"""
Contract validation module for PGx Pipeline.

Provides validation functions for intermediate data between pipeline stages.
Ensures deterministic behavior and clinical defensibility.
"""

import hashlib
import json
from pathlib import Path
from typing import Any

import jsonschema


def load_schema(schema_path: Path) -> dict:
    """Load JSON schema from file."""
    with open(schema_path) as f:
        return json.load(f)


def validate_json_schema(data: dict, schema: dict) -> list[str]:
    """
    Validate data against JSON schema.

    Args:
        data: Dictionary to validate.
        schema: JSON schema dict.

    Returns:
        List of validation error messages (empty if valid).
    """
    validator = jsonschema.Draft202012Validator(schema)
    errors = []
    for error in validator.iter_errors(data):
        path = ".".join(str(p) for p in error.path) or "root"
        errors.append(f"{path}: {error.message}")
    return errors


def validate_final_report(report_data: dict, schema_path: Path) -> list[str]:
    """
    Validate final PGx report against schema.

    Args:
        report_data: Report dictionary.
        schema_path: Path to pgx_report_schema.json.

    Returns:
        List of validation errors (empty if valid).
    """
    schema = load_schema(schema_path)
    return validate_json_schema(report_data, schema)


def calculate_file_md5(file_path: Path) -> str:
    """Calculate MD5 hash of a file."""
    md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            md5.update(chunk)
    return md5.hexdigest()


def calculate_dir_checksum(dir_path: Path) -> str:
    """Calculate combined checksum of all files in directory."""
    md5 = hashlib.md5()
    for file_path in sorted(dir_path.glob("**/*")):
        if file_path.is_file():
            md5.update(file_path.name.encode())
            md5.update(calculate_file_md5(file_path).encode())
    return md5.hexdigest()


# --- Contract Schemas ---

VARIANTS_META_SCHEMA = {
    "type": "object",
    "required": ["caller", "caller_version", "sample_id", "reference_build"],
    "properties": {
        "caller": {"type": "string"},
        "caller_version": {"type": "string"},
        "reference_build": {"type": "string"},
        "reference_md5": {"type": "string"},
        "regions_bed_hash": {"type": "string"},
        "sample_id": {"type": "string"},
    },
}

ALLELES_JSON_SCHEMA = {
    "type": "object",
    "required": ["sample_id", "reference", "calls"],
    "properties": {
        "sample_id": {"type": "string"},
        "reference": {"type": "string"},
        "calls": {
            "type": "array",
            "items": {
                "type": "object",
                "required": ["gene", "diplotype", "call_quality"],
                "properties": {
                    "gene": {"type": "string"},
                    "diplotype": {"type": "string"},
                    "haplotypes": {"type": "array", "items": {"type": "string"}},
                    "call_method": {"type": "string"},
                    "call_quality": {
                        "type": "object",
                        "required": ["score"],
                        "properties": {
                            "score": {"type": "number", "minimum": 0, "maximum": 1},
                            "warnings": {"type": "array", "items": {"type": "string"}},
                        },
                    },
                    "evidence": {
                        "type": "object",
                        "properties": {
                            "variants_used": {"type": "array"},
                        },
                    },
                },
            },
        },
    },
}

PHENOTYPES_JSON_SCHEMA = {
    "type": "object",
    "required": ["sample_id", "phenotypes"],
    "properties": {
        "sample_id": {"type": "string"},
        "phenotypes": {
            "type": "array",
            "items": {
                "type": "object",
                "required": ["gene", "diplotype", "phenotype"],
                "properties": {
                    "gene": {"type": "string"},
                    "diplotype": {"type": "string"},
                    "phenotype": {"type": "string"},
                    "activity_score": {"type": ["number", "null"]},
                    "functional_status": {"type": "string"},
                    "prescribing_implication": {"type": "string"},
                },
            },
        },
    },
}


def validate_variants_meta(data: dict) -> list[str]:
    """Validate variants.meta.json contract."""
    return validate_json_schema(data, VARIANTS_META_SCHEMA)


def validate_alleles_json(data: dict) -> list[str]:
    """Validate alleles.json contract."""
    return validate_json_schema(data, ALLELES_JSON_SCHEMA)


def validate_phenotypes_json(data: dict) -> list[str]:
    """Validate phenotypes.json contract."""
    return validate_json_schema(data, PHENOTYPES_JSON_SCHEMA)


def create_provenance(
    pipeline_version: str,
    git_commit: str | None = None,
    nextflow_version: str | None = None,
    docker_images: list[str] | None = None,
    reference_path: Path | None = None,
    pgx_db_dir: Path | None = None,
) -> dict[str, Any]:
    """
    Create provenance metadata for output.

    Args:
        pipeline_version: Pipeline version string.
        git_commit: Git commit hash.
        nextflow_version: Nextflow version.
        docker_images: List of Docker image names.
        reference_path: Path to reference FASTA.
        pgx_db_dir: Path to PGx database directory.

    Returns:
        Provenance dictionary.
    """
    provenance: dict[str, Any] = {
        "pipeline_name": "pgx-pipeline",
        "pipeline_version": pipeline_version,
    }

    if git_commit:
        provenance["git_commit"] = git_commit

    if nextflow_version:
        provenance["nextflow_version"] = nextflow_version

    if docker_images:
        provenance["docker_images"] = docker_images

    if reference_path and reference_path.exists():
        provenance["reference_md5"] = calculate_file_md5(reference_path)

    if pgx_db_dir and pgx_db_dir.exists():
        provenance["pgx_db_checksum"] = calculate_dir_checksum(pgx_db_dir)

    return provenance
