"""
Utility functions for PGx Pipeline.

Shared helpers used across multiple modules.
"""

import json
import logging
from pathlib import Path
from typing import Any

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)

logger = logging.getLogger("pgx")


def setup_logging(level: str = "INFO") -> None:
    """
    Configure logging level.

    Args:
        level: Logging level (DEBUG, INFO, WARNING, ERROR).
    """
    numeric_level = getattr(logging, level.upper(), logging.INFO)
    logger.setLevel(numeric_level)


def load_json(path: Path) -> dict[str, Any]:
    """
    Load JSON file with error handling.

    Args:
        path: Path to JSON file.

    Returns:
        Parsed JSON as dict.

    Raises:
        FileNotFoundError: If file doesn't exist.
        json.JSONDecodeError: If JSON is invalid.
    """
    if not path.exists():
        raise FileNotFoundError(f"JSON file not found: {path}")

    with open(path) as f:
        return json.load(f)


def save_json(data: dict[str, Any], path: Path, indent: int = 2) -> None:
    """
    Save dict to JSON file.

    Args:
        data: Dictionary to save.
        path: Output path.
        indent: JSON indentation level.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(data, f, indent=indent, default=str)


def normalize_chromosome(chrom: str) -> str:
    """
    Normalize chromosome name to hg38 format (with 'chr' prefix).

    Args:
        chrom: Chromosome name (e.g., "10", "chr10", "Chr10").

    Returns:
        Normalized chromosome name (e.g., "chr10").
    """
    chrom = chrom.strip()

    # Remove existing prefix and re-add consistently
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]

    return f"chr{chrom}"


def parse_diplotype(diplotype: str) -> tuple[str, str]:
    """
    Parse diplotype string into two alleles.

    Args:
        diplotype: Diplotype string (e.g., "*1/*17", "*2/*3").

    Returns:
        Tuple of (allele1, allele2).
    """
    if "/" not in diplotype:
        return diplotype, diplotype

    parts = diplotype.split("/")
    if len(parts) != 2:
        return parts[0], parts[0]

    return parts[0].strip(), parts[1].strip()


def format_diplotype(allele1: str, allele2: str) -> str:
    """
    Format two alleles into canonical diplotype string.

    Alleles are sorted alphabetically for consistency.

    Args:
        allele1: First allele.
        allele2: Second allele.

    Returns:
        Formatted diplotype (e.g., "*1/*17").
    """
    # Sort for consistent representation
    sorted_alleles = sorted([allele1, allele2])
    return f"{sorted_alleles[0]}/{sorted_alleles[1]}"


def get_genes_for_tier(gene_tiers: dict, tier: str) -> list[str]:
    """
    Get list of genes belonging to a specific tier.

    Args:
        gene_tiers: Gene tiers dict from gene_tiers.json.
        tier: Tier name (e.g., "group1").

    Returns:
        List of gene names in that tier.
    """
    genes = []
    for gene, info in gene_tiers.get("genes", {}).items():
        if info.get("tier") == tier:
            genes.append(gene)
    return genes


def validate_gene_tiers_config(gene_tiers_file: Path) -> list[str]:
    """
    Validate gene_tiers.json configuration.

    Args:
        gene_tiers_file: Path to gene_tiers.json.

    Returns:
        List of validation error messages (empty if valid).
    """
    errors = []

    if not gene_tiers_file.exists():
        errors.append(f"Gene tiers file not found: {gene_tiers_file}")
        return errors

    try:
        data = load_json(gene_tiers_file)
    except json.JSONDecodeError as e:
        errors.append(f"Invalid JSON in gene tiers file: {e}")
        return errors

    if "genes" not in data:
        errors.append("Missing 'genes' key in gene tiers file")
        return errors

    valid_tiers = {"group1", "group2", "group3", "group4"}
    valid_biomarkers = {"metabolic_enzyme", "transporter", "drug_target", "hla", "other"}

    for gene, info in data["genes"].items():
        tier = info.get("tier")
        if tier not in valid_tiers:
            errors.append(f"Invalid tier '{tier}' for gene {gene}")

        biomarker = info.get("biomarker_type")
        if biomarker not in valid_biomarkers:
            errors.append(f"Invalid biomarker_type '{biomarker}' for gene {gene}")

    return errors
