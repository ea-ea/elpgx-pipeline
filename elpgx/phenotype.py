"""
Phenotype translation module for PGx Pipeline.

Maps diplotypes to phenotypes and prescribing implications:
- Diplotype → Activity Score (where applicable)
- Activity Score → Phenotype
- Phenotype → Prescribing Implication
"""

import json
from dataclasses import dataclass
from pathlib import Path

from elpgx.allele_caller import DiplotypeCall


@dataclass
class PhenotypeResult:
    """Result of phenotype translation."""

    gene: str
    diplotype: str
    activity_score: float | None
    functional_status: str
    phenotype: str
    prescribing_implication: str


def load_phenotype_map(phenotype_maps_dir: Path, gene: str) -> dict | None:
    """
    Load phenotype mapping for a gene.

    Args:
        phenotype_maps_dir: Path to phenotype maps directory.
        gene: Gene name.

    Returns:
        Phenotype map dict or None if not found.
    """
    map_file = phenotype_maps_dir / f"{gene}_phenotype.json"
    if not map_file.exists():
        return None

    with open(map_file) as f:
        return json.load(f)


def load_prescribing_implications(
    implications_dir: Path, gene: str
) -> dict[str, str]:
    """
    Load prescribing implications for a gene.

    Args:
        implications_dir: Path to prescribing implications directory.
        gene: Gene name.

    Returns:
        Dict mapping phenotype to implication text.
    """
    impl_file = implications_dir / f"{gene}_implications.json"
    if not impl_file.exists():
        return {}

    with open(impl_file) as f:
        data = json.load(f)

    implications = {}
    for item in data.get("implications", []):
        phenotype = item.get("phenotype")
        implication = item.get("implication")
        if phenotype and implication:
            implications[phenotype] = implication

    return implications


def calculate_activity_score(
    diplotype_call: DiplotypeCall,
) -> float | None:
    """
    Calculate activity score from diplotype call.

    Sum of individual allele activity values.

    Args:
        diplotype_call: DiplotypeCall with allele activities.

    Returns:
        Activity score or None if not applicable.
    """
    act1 = diplotype_call.allele1_activity
    act2 = diplotype_call.allele2_activity

    if act1 is None or act2 is None:
        return None

    return act1 + act2


def get_functional_status(diplotype_call: DiplotypeCall) -> str:
    """
    Get combined functional status from diplotype call.

    Args:
        diplotype_call: DiplotypeCall with allele functions.

    Returns:
        Combined functional status string.
    """
    func1 = diplotype_call.allele1_function
    func2 = diplotype_call.allele2_function

    if func1 == func2:
        return func1

    return f"{func1} / {func2}"


def lookup_phenotype_by_diplotype(
    diplotype: str,
    phenotype_map: dict,
) -> tuple[str | None, float | None, str | None]:
    """
    Look up phenotype by exact diplotype match.

    Args:
        diplotype: Diplotype string (e.g., "*1/*17").
        phenotype_map: Phenotype map dict.

    Returns:
        Tuple of (phenotype, activity_score, functional_status) or (None, None, None).
    """
    mappings = phenotype_map.get("mappings", [])

    # Try exact match
    for mapping in mappings:
        if mapping.get("diplotype") == diplotype:
            return (
                mapping.get("phenotype"),
                mapping.get("activity_score"),
                mapping.get("functional_status"),
            )

    # Try reversed diplotype (e.g., "*17/*1" -> "*1/*17")
    parts = diplotype.split("/")
    if len(parts) == 2:
        reversed_diplotype = f"{parts[1]}/{parts[0]}"
        for mapping in mappings:
            if mapping.get("diplotype") == reversed_diplotype:
                return (
                    mapping.get("phenotype"),
                    mapping.get("activity_score"),
                    mapping.get("functional_status"),
                )

    return None, None, None


def lookup_phenotype_by_activity_score(
    activity_score: float,
    phenotype_map: dict,
) -> str | None:
    """
    Look up phenotype by activity score.

    Args:
        activity_score: Calculated activity score.
        phenotype_map: Phenotype map dict.

    Returns:
        Phenotype string or None.
    """
    score_map = phenotype_map.get("phenotype_by_activity_score", {})

    # Try exact match
    score_str = str(activity_score)
    if score_str in score_map:
        return score_map[score_str]

    # Try integer version
    if activity_score == int(activity_score):
        int_str = str(int(activity_score))
        if int_str in score_map:
            return score_map[int_str]

    # Find closest score
    available_scores = []
    for key in score_map.keys():
        try:
            available_scores.append(float(key))
        except ValueError:
            continue

    if available_scores:
        closest = min(available_scores, key=lambda x: abs(x - activity_score))
        return score_map.get(str(closest)) or score_map.get(str(int(closest)))

    return None


def translate_diplotype_to_phenotype(
    diplotype_call: DiplotypeCall,
    phenotype_maps_dir: Path,
    implications_dir: Path,
) -> PhenotypeResult:
    """
    Translate a diplotype call to phenotype and prescribing implication.

    Args:
        diplotype_call: DiplotypeCall from allele caller.
        phenotype_maps_dir: Path to phenotype maps directory.
        implications_dir: Path to prescribing implications directory.

    Returns:
        PhenotypeResult with phenotype and implication.
    """
    gene = diplotype_call.gene
    diplotype = diplotype_call.diplotype

    # Load phenotype map
    phenotype_map = load_phenotype_map(phenotype_maps_dir, gene)

    # Calculate activity score
    activity_score = calculate_activity_score(diplotype_call)

    # Get functional status
    functional_status = get_functional_status(diplotype_call)

    # Determine phenotype
    phenotype: str | None = None

    if phenotype_map:
        # Try diplotype lookup first
        looked_up_phenotype, looked_up_score, looked_up_status = (
            lookup_phenotype_by_diplotype(diplotype, phenotype_map)
        )

        if looked_up_phenotype:
            phenotype = looked_up_phenotype
            if looked_up_score is not None:
                activity_score = looked_up_score
            if looked_up_status:
                functional_status = looked_up_status
        elif activity_score is not None and phenotype_map.get("activity_score_based"):
            # Try activity score lookup
            phenotype = lookup_phenotype_by_activity_score(activity_score, phenotype_map)

    # Fallback phenotype based on functional status
    if phenotype is None:
        phenotype = infer_phenotype_from_function(functional_status, gene)

    # Load prescribing implication
    implications = load_prescribing_implications(implications_dir, gene)
    prescribing_implication = implications.get(phenotype, "No specific guidance available.")

    return PhenotypeResult(
        gene=gene,
        diplotype=diplotype,
        activity_score=activity_score,
        functional_status=functional_status,
        phenotype=phenotype,
        prescribing_implication=prescribing_implication,
    )


def infer_phenotype_from_function(functional_status: str, gene: str) -> str:
    """
    Infer phenotype from functional status when no direct mapping exists.

    Args:
        functional_status: Combined functional status string.
        gene: Gene name (for biomarker-type specific labels).

    Returns:
        Inferred phenotype string.
    """
    status_lower = functional_status.lower()

    # Check for specific patterns
    if "no function" in status_lower:
        if status_lower.count("no function") == 2:
            return "Poor Metabolizer"
        return "Intermediate Metabolizer"

    if "decreased function" in status_lower:
        if status_lower.count("decreased function") == 2:
            return "Intermediate Metabolizer"
        return "Intermediate Metabolizer"

    if "increased function" in status_lower:
        if status_lower.count("increased function") == 2:
            return "Ultrarapid Metabolizer"
        if "normal function" in status_lower:
            return "Rapid Metabolizer"
        return "Rapid Metabolizer"

    if "normal function" in status_lower:
        return "Normal Metabolizer"

    # For transporters/drug targets
    if "sensitivity" in status_lower:
        if "high" in status_lower:
            return "High Sensitivity"
        if "increased" in status_lower:
            return "Increased Sensitivity"
        return "Normal Sensitivity"

    # Default
    return "Indeterminate"


def translate_all_calls(
    diplotype_calls: list[DiplotypeCall],
    phenotype_maps_dir: Path,
    implications_dir: Path,
) -> list[PhenotypeResult]:
    """
    Translate all diplotype calls to phenotypes.

    Args:
        diplotype_calls: List of DiplotypeCall from allele caller.
        phenotype_maps_dir: Path to phenotype maps directory.
        implications_dir: Path to prescribing implications directory.

    Returns:
        List of PhenotypeResult for each gene.
    """
    results = []
    for call in diplotype_calls:
        result = translate_diplotype_to_phenotype(
            call,
            phenotype_maps_dir,
            implications_dir,
        )
        results.append(result)
    return results
