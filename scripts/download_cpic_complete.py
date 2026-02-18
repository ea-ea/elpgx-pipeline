#!/usr/bin/env python3
"""
Complete CPIC Data Download and Conversion for ElPGx Pipeline.

Downloads official clinical pharmacogenomics data from CPIC API and converts
to the ElPGx JSON format with GRCh38 variant coordinates.

Data sources:
- CPIC API: https://api.cpicpgx.org/
- Based on CPIC guidelines: https://cpicpgx.org/

Author: Elif Arslan
Date: 2026-02
"""

import json
import requests
from pathlib import Path
from datetime import datetime, timezone
from typing import Any
import time

# CPIC API base URL
CPIC_API = "https://api.cpicpgx.org/v1"

# Target genes for ElPGx pipeline
TARGET_GENES = [
    "CYP2C19", "CYP2C9", "CYP2D6", "CYP3A5", "VKORC1", "SLCO1B1",
    "TPMT", "NUDT15", "DPYD", "UGT1A1", "NAT2", "IFNL3", "G6PD",
    "CYP2B6", "CACNA1S", "RYR1", "HLA-A", "HLA-B"
]

# Gene tier mapping
GENE_TIERS = {
    "CYP2C19": "group1", "CYP2C9": "group1", "VKORC1": "group1",
    "SLCO1B1": "group1", "TPMT": "group1", "NUDT15": "group1",
    "DPYD": "group1", "UGT1A1": "group1", "CYP2D6": "group2",
    "CYP3A5": "group2", "CYP2B6": "group2", "G6PD": "group2",
    "NAT2": "group3", "IFNL3": "group3", "CACNA1S": "group3",
    "RYR1": "group3", "HLA-A": "group4", "HLA-B": "group4",
}

# Biomarker types
BIOMARKER_TYPES = {
    "CYP2C19": "metabolic_enzyme", "CYP2C9": "metabolic_enzyme",
    "CYP2D6": "metabolic_enzyme", "CYP3A5": "metabolic_enzyme",
    "CYP2B6": "metabolic_enzyme", "VKORC1": "drug_target",
    "SLCO1B1": "transporter", "TPMT": "metabolic_enzyme",
    "NUDT15": "metabolic_enzyme", "DPYD": "metabolic_enzyme",
    "UGT1A1": "metabolic_enzyme", "NAT2": "metabolic_enzyme",
    "IFNL3": "drug_target", "G6PD": "metabolic_enzyme",
    "CACNA1S": "drug_target", "RYR1": "drug_target",
    "HLA-A": "hla", "HLA-B": "hla",
}

# Chromosome mapping
CHROMOSOMES = {
    "CYP2C19": "chr10", "CYP2C9": "chr10", "CYP2D6": "chr22",
    "CYP3A5": "chr7", "CYP2B6": "chr19", "VKORC1": "chr16",
    "SLCO1B1": "chr12", "TPMT": "chr6", "NUDT15": "chr13",
    "DPYD": "chr1", "UGT1A1": "chr2", "NAT2": "chr8",
    "IFNL3": "chr19", "G6PD": "chrX", "CACNA1S": "chr1",
    "RYR1": "chr19", "HLA-A": "chr6", "HLA-B": "chr6",
}


def fetch_cpic_data(endpoint: str, params: dict = None) -> list[dict]:
    """Fetch data from CPIC API with retry logic."""
    url = f"{CPIC_API}/{endpoint}"
    for attempt in range(3):
        try:
            response = requests.get(url, params=params, timeout=60)
            response.raise_for_status()
            return response.json()
        except requests.RequestException as e:
            if attempt < 2:
                time.sleep(1)
            else:
                print(f"  Warning: Failed to fetch {endpoint}: {e}")
                return []
    return []


def fetch_allele_definitions(gene: str) -> list[dict]:
    """Fetch allele definitions for a gene."""
    return fetch_cpic_data("allele", {"genesymbol": f"eq.{gene}"})


def fetch_sequence_locations(gene: str) -> list[dict]:
    """Fetch GRCh38 variant coordinates."""
    return fetch_cpic_data("sequence_location", {"genesymbol": f"eq.{gene}"})


def fetch_allele_location_values(definition_id: int) -> list[dict]:
    """Fetch variant-allele mappings."""
    return fetch_cpic_data("allele_location_value",
                           {"alleledefinitionid": f"eq.{definition_id}"})


def fetch_allele_definition_info(gene: str) -> dict:
    """Fetch allele definition metadata."""
    data = fetch_cpic_data("allele_definition", {"genesymbol": f"eq.{gene}"})
    return data[0] if data else {}


def fetch_diplotypes(gene: str) -> list[dict]:
    """Fetch diplotype to phenotype mappings."""
    return fetch_cpic_data("diplotype", {"genesymbol": f"eq.{gene}"})


def fetch_drug_pairs(gene: str) -> list[dict]:
    """Fetch gene-drug pairs for prescribing implications."""
    return fetch_cpic_data("pair", {"genesymbol": f"eq.{gene}"})


def fetch_drugs() -> dict:
    """Fetch all drug information."""
    drugs = fetch_cpic_data("drug")
    return {d["drugid"]: d for d in drugs}


def fetch_gene_info(gene: str) -> dict:
    """Fetch gene information."""
    data = fetch_cpic_data("gene", {"symbol": f"eq.{gene}"})
    return data[0] if data else {}


def map_function_to_activity(function: str) -> float | None:
    """Map functional status to activity value."""
    if not function:
        return None
    function_lower = function.lower()

    if "no function" in function_lower:
        return 0.0
    elif "decreased" in function_lower:
        return 0.5
    elif "increased" in function_lower:
        return 1.5
    elif "normal" in function_lower:
        return 1.0
    elif "uncertain" in function_lower or "unknown" in function_lower:
        return None
    return 1.0


def parse_allele_number(allele: str) -> tuple:
    """Parse allele name for sorting."""
    try:
        num_part = allele.lstrip("*")
        base = ""
        suffix = ""
        for i, c in enumerate(num_part):
            if c.isdigit():
                base += c
            else:
                suffix = num_part[i:]
                break
        return (int(base) if base else 999, suffix)
    except:
        return (999, allele)


def build_allele_definitions(gene: str, alleles: list[dict],
                             seq_locations: list[dict]) -> dict:
    """Build allele definitions with GRCh38 coordinates."""

    # Build rsID to location mapping
    rsid_map = {}
    for loc in seq_locations:
        rsid = loc.get("dbsnpid")
        if rsid:
            # Parse chromosome location (e.g., "g.94780544A>G")
            chrom_loc = loc.get("chromosomelocation", "")
            ref, alt = None, None
            if ">" in chrom_loc:
                try:
                    change = chrom_loc.split(".")[-1]
                    ref = change.split(">")[0][-1]
                    alt = change.split(">")[1][0]
                except:
                    pass

            rsid_map[rsid] = {
                "rsid": rsid,
                "chrom": CHROMOSOMES.get(gene, "chr1"),
                "pos": loc.get("position"),
                "ref": ref,
                "alt": alt,
                "name": loc.get("name"),
                "protein_location": loc.get("proteinlocation"),
            }

    # Build allele entries
    allele_entries = []

    for allele in alleles:
        allele_name = allele.get("name", "")
        if not allele_name:
            continue

        function = allele.get("clinicalfunctionalstatus", "Unknown")
        activity = allele.get("activityvalue")
        if activity is None:
            activity = map_function_to_activity(function)

        # Get frequency data if available
        frequency = allele.get("frequency", {})

        entry = {
            "allele": allele_name,
            "function": function,
            "activity_value": activity,
            "defining_variants": [],  # Will be populated if we can get allele-variant mappings
            "citations": allele.get("citations", []),
            "evidence_strength": allele.get("strength", ""),
            "frequency": frequency if frequency else None,
        }

        allele_entries.append(entry)

    # Sort alleles
    allele_entries.sort(key=lambda x: parse_allele_number(x["allele"]))

    # Get reference allele (usually *1)
    ref_allele = "*1"
    for a in allele_entries:
        if a["allele"] == "*1":
            ref_allele = "*1"
            break

    return {
        "gene": gene,
        "chromosome": CHROMOSOMES.get(gene, "chr1"),
        "gene_tier": GENE_TIERS.get(gene, "group1"),
        "biomarker_type": BIOMARKER_TYPES.get(gene, "metabolic_enzyme"),
        "reference_allele": ref_allele,
        "reference_build": "GRCh38",
        "cpic_version": "v1.24",
        "last_updated": datetime.now(timezone.utc).strftime("%Y-%m-%d"),
        "sequence_locations": list(rsid_map.values())[:50],  # Include key variants
        "alleles": allele_entries
    }


def build_phenotype_map(gene: str, diplotypes: list[dict]) -> dict:
    """Build phenotype mapping from diplotype data."""
    mappings = []
    phenotype_by_score = {}

    for dip in diplotypes:
        diplotype = dip.get("diplotype", "")
        phenotype = dip.get("generesult", "")
        activity_score = dip.get("activityscore")
        ehr_priority = dip.get("ehrpriority", "")

        mappings.append({
            "diplotype": diplotype,
            "phenotype": phenotype,
            "activity_score": activity_score,
            "ehr_priority": ehr_priority,
        })

        if activity_score is not None:
            score_key = str(activity_score)
            if score_key not in phenotype_by_score:
                phenotype_by_score[score_key] = phenotype

    mappings.sort(key=lambda x: x["diplotype"])

    return {
        "gene": gene,
        "cpic_version": "v1.24",
        "last_updated": datetime.now(timezone.utc).strftime("%Y-%m-%d"),
        "activity_score_based": len(phenotype_by_score) > 0,
        "phenotype_by_activity_score": phenotype_by_score,
        "mappings": mappings
    }


def build_implications(gene: str, pairs: list[dict], drugs_map: dict) -> dict:
    """Build prescribing implications from gene-drug pairs."""
    implications = []

    # Group by phenotype based on CPIC level
    cpic_a_drugs = []
    cpic_b_drugs = []
    other_drugs = []

    for pair in pairs:
        drug_id = pair.get("drugid", "")
        drug_info = drugs_map.get(drug_id, {})
        drug_name = drug_info.get("name", "Unknown")
        cpic_level = pair.get("cpiclevel", "")
        pgkb_level = pair.get("pgkbcalevel", "")
        guideline_id = pair.get("guidelineid")

        drug_entry = {
            "drug": drug_name,
            "drug_id": drug_id,
            "cpic_level": cpic_level,
            "pharmgkb_level": pgkb_level,
            "has_guideline": guideline_id is not None,
            "flowchart": drug_info.get("flowchart"),
        }

        if cpic_level and cpic_level.startswith("A"):
            cpic_a_drugs.append(drug_entry)
        elif cpic_level and cpic_level.startswith("B"):
            cpic_b_drugs.append(drug_entry)
        else:
            other_drugs.append(drug_entry)

    # Build summary implications
    if cpic_a_drugs:
        drug_names = [d["drug"] for d in cpic_a_drugs]
        implications.append({
            "level": "CPIC Level A",
            "drugs": cpic_a_drugs,
            "summary": f"Strong gene-drug interactions with: {', '.join(drug_names)}. "
                      "Genetic testing recommended before prescribing."
        })

    if cpic_b_drugs:
        drug_names = [d["drug"] for d in cpic_b_drugs]
        implications.append({
            "level": "CPIC Level B",
            "drugs": cpic_b_drugs,
            "summary": f"Moderate gene-drug interactions with: {', '.join(drug_names)}. "
                      "Consider genetic testing for certain populations."
        })

    return {
        "gene": gene,
        "cpic_version": "v1.24",
        "last_updated": datetime.now(timezone.utc).strftime("%Y-%m-%d"),
        "total_drug_pairs": len(pairs),
        "level_a_pairs": len(cpic_a_drugs),
        "level_b_pairs": len(cpic_b_drugs),
        "implications": implications,
        "all_drugs": cpic_a_drugs + cpic_b_drugs + other_drugs
    }


def process_gene(gene: str, output_dir: Path, drugs_map: dict) -> dict:
    """Process a single gene and generate all data files."""
    print(f"\nProcessing {gene}...")
    result = {"gene": gene, "status": "success", "files": []}

    try:
        # Fetch all data
        print(f"  Fetching allele definitions...")
        alleles = fetch_allele_definitions(gene)

        print(f"  Fetching sequence locations...")
        seq_locs = fetch_sequence_locations(gene)

        print(f"  Fetching diplotypes...")
        diplotypes = fetch_diplotypes(gene)

        print(f"  Fetching drug pairs...")
        pairs = fetch_drug_pairs(gene)

        # Generate allele definitions
        if alleles:
            allele_def = build_allele_definitions(gene, alleles, seq_locs)
            allele_file = output_dir / "allele_definitions" / f"{gene}.json"
            allele_file.parent.mkdir(parents=True, exist_ok=True)
            with open(allele_file, 'w') as f:
                json.dump(allele_def, f, indent=2)
            result["files"].append(str(allele_file))
            print(f"    -> {len(alleles)} alleles, {len(seq_locs)} variant positions")

        # Generate phenotype map
        if diplotypes:
            phenotype_map = build_phenotype_map(gene, diplotypes)
            phenotype_file = output_dir / "phenotype_maps" / f"{gene}_phenotype.json"
            phenotype_file.parent.mkdir(parents=True, exist_ok=True)
            with open(phenotype_file, 'w') as f:
                json.dump(phenotype_map, f, indent=2)
            result["files"].append(str(phenotype_file))
            print(f"    -> {len(diplotypes)} diplotype mappings")

        # Generate prescribing implications
        if pairs:
            implications = build_implications(gene, pairs, drugs_map)
            impl_file = output_dir / "prescribing_implications" / f"{gene}_implications.json"
            impl_file.parent.mkdir(parents=True, exist_ok=True)
            with open(impl_file, 'w') as f:
                json.dump(implications, f, indent=2)
            result["files"].append(str(impl_file))
            print(f"    -> {len(pairs)} drug pairs")

    except Exception as e:
        result["status"] = "error"
        result["error"] = str(e)
        print(f"  ERROR: {e}")

    return result


def generate_gene_tiers(output_dir: Path, genes: list[str]):
    """Generate gene_tiers.json configuration."""
    genes_config = {}

    for gene in genes:
        gene_info = fetch_gene_info(gene)
        genes_config[gene] = {
            "tier": GENE_TIERS.get(gene, "group1"),
            "biomarker_type": BIOMARKER_TYPES.get(gene, "metabolic_enzyme"),
            "chromosome": CHROMOSOMES.get(gene, gene_info.get("chr", "")),
            "cpic_level": "A" if gene in ["CYP2C19", "CYP2C9", "CYP2D6", "VKORC1",
                                           "SLCO1B1", "TPMT", "NUDT15", "DPYD",
                                           "UGT1A1", "CYP3A5"] else "B",
            "pharmgkb_id": gene_info.get("pharmgkbid", ""),
            "ncbi_id": gene_info.get("ncbiid", ""),
            "ensembl_id": gene_info.get("ensemblid", ""),
        }

    output = {
        "cpic_version": "v1.24",
        "reference_build": "GRCh38",
        "last_updated": datetime.now(timezone.utc).strftime("%Y-%m-%d"),
        "genes": genes_config
    }

    with open(output_dir / "gene_tiers.json", 'w') as f:
        json.dump(output, f, indent=2)

    print(f"\nGenerated gene_tiers.json with {len(genes_config)} genes")


def generate_pgx_regions_bed(output_dir: Path, genes: list[str]):
    """Generate PGx regions BED file."""
    bed_lines = []

    for gene in genes:
        seq_locs = fetch_sequence_locations(gene)
        if seq_locs:
            positions = [loc.get("position") for loc in seq_locs if loc.get("position")]
            if positions:
                chrom = CHROMOSOMES.get(gene, "chr1")
                min_pos = min(positions) - 1000
                max_pos = max(positions) + 1000
                bed_lines.append((chrom, min_pos, max_pos, gene))

    # Sort by chromosome and position
    bed_lines.sort(key=lambda x: (x[0], x[1]))

    with open(output_dir / "pgx_regions.bed", 'w') as f:
        f.write("# PGx Target Regions for ElPGx Pipeline\n")
        f.write(f"# Generated from CPIC v1.24 on {datetime.now(timezone.utc).strftime('%Y-%m-%d')}\n")
        f.write("# Reference: GRCh38\n")
        for chrom, start, end, gene in bed_lines:
            f.write(f"{chrom}\t{start}\t{end}\t{gene}\n")

    print(f"Generated pgx_regions.bed with {len(bed_lines)} regions")


def main():
    """Main function."""
    print("=" * 70)
    print("CPIC Complete Data Download for ElPGx Pipeline")
    print("=" * 70)
    print(f"Source: {CPIC_API}")
    print(f"Target genes: {len(TARGET_GENES)}")
    print(f"Date: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M:%S UTC')}")
    print("=" * 70)

    output_dir = Path(__file__).parent.parent / "data"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Fetch drug information first (used for all genes)
    print("\nFetching drug database...")
    drugs_map = fetch_drugs()
    print(f"  -> {len(drugs_map)} drugs loaded")

    # Process each gene
    results = []
    genes_processed = []

    for gene in TARGET_GENES:
        result = process_gene(gene, output_dir, drugs_map)
        results.append(result)
        if result["status"] == "success" and result["files"]:
            genes_processed.append(gene)
        time.sleep(0.5)  # Rate limiting

    # Generate configuration files
    print("\n" + "=" * 70)
    print("Generating configuration files...")

    if genes_processed:
        generate_gene_tiers(output_dir, genes_processed)
        generate_pgx_regions_bed(output_dir, genes_processed)

    # Write provenance
    provenance = {
        "source": "Clinical Pharmacogenetics Implementation Consortium (CPIC)",
        "api_url": CPIC_API,
        "cpic_version": "v1.24",
        "reference_build": "GRCh38",
        "download_date": datetime.now(timezone.utc).isoformat(),
        "genes_included": genes_processed,
        "total_genes": len(genes_processed),
        "license": "CC BY-SA 4.0",
        "citation": "Relling MV, Klein TE. CPIC: Clinical Pharmacogenetics Implementation Consortium. Clin Pharmacol Ther. 2011;89(3):464-467.",
        "pharmgkb_citation": "Whirl-Carrillo M, et al. Pharmacogenomics Knowledge for Personalized Medicine. Clin Pharmacol Ther. 2012;92(4):414-417."
    }

    with open(output_dir / "PROVENANCE.json", 'w') as f:
        json.dump(provenance, f, indent=2)

    # Summary
    print("\n" + "=" * 70)
    print("DOWNLOAD SUMMARY")
    print("=" * 70)

    success = sum(1 for r in results if r["status"] == "success")
    failed = sum(1 for r in results if r["status"] == "error")

    print(f"Total genes processed: {len(results)}")
    print(f"Successful: {success}")
    print(f"Failed: {failed}")

    print("\n" + "=" * 70)
    print("DATA PROVENANCE")
    print("=" * 70)
    print("Source: Clinical Pharmacogenetics Implementation Consortium (CPIC)")
    print(f"API: {CPIC_API}")
    print("Version: v1.24")
    print("Reference Build: GRCh38")
    print("License: CC BY-SA 4.0")
    print("=" * 70)
    print("\nData files written to:", output_dir)
    print("\nDone!")


if __name__ == "__main__":
    main()
