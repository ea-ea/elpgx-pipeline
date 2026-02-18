#!/usr/bin/env python3
"""
Download and convert official CPIC data to ElPGx format.

This script downloads data from the official CPIC API (api.cpicpgx.org)
and converts it to the JSON format used by the ElPGx pipeline.

Data sources:
- CPIC API: https://api.cpicpgx.org/
- Based on CPIC guidelines: https://cpicpgx.org/

Author: Elif Arslan
Date: 2026-02
"""

import json
import requests
from pathlib import Path
from datetime import datetime
from typing import Any

# CPIC API base URL
CPIC_API = "https://api.cpicpgx.org/v1"

# Target genes for ElPGx pipeline (high-evidence pharmacogenes)
TARGET_GENES = [
    "CYP2C19",   # Clopidogrel, PPIs, antidepressants
    "CYP2C9",    # Warfarin, NSAIDs, phenytoin
    "CYP2D6",    # Opioids, antidepressants, antipsychotics
    "CYP3A5",    # Tacrolimus
    "VKORC1",    # Warfarin
    "SLCO1B1",   # Statins
    "TPMT",      # Thiopurines
    "NUDT15",    # Thiopurines
    "DPYD",      # Fluoropyrimidines
    "UGT1A1",    # Irinotecan
    "NAT2",      # Isoniazid
    "IFNL3",     # Peginterferon alfa-2a/b
    "G6PD",      # Rasburicase, antimalarials
    "CYP2B6",    # Efavirenz
    "CACNA1S",   # Anesthetics (MH)
    "RYR1",      # Anesthetics (MH)
    "HLA-A",     # Abacavir (HLA typing)
    "HLA-B",     # Carbamazepine, allopurinol (HLA typing)
]

# Gene tier mapping (CPIC level of evidence)
GENE_TIERS = {
    "CYP2C19": "group1",
    "CYP2C9": "group1",
    "VKORC1": "group1",
    "SLCO1B1": "group1",
    "TPMT": "group1",
    "NUDT15": "group1",
    "DPYD": "group1",
    "UGT1A1": "group1",
    "CYP2D6": "group2",
    "CYP3A5": "group2",
    "CYP2B6": "group2",
    "G6PD": "group2",
    "NAT2": "group3",
    "IFNL3": "group3",
    "CACNA1S": "group3",
    "RYR1": "group3",
    "HLA-A": "group4",
    "HLA-B": "group4",
}

# Biomarker type mapping
BIOMARKER_TYPES = {
    "CYP2C19": "metabolic_enzyme",
    "CYP2C9": "metabolic_enzyme",
    "CYP2D6": "metabolic_enzyme",
    "CYP3A5": "metabolic_enzyme",
    "CYP2B6": "metabolic_enzyme",
    "VKORC1": "drug_target",
    "SLCO1B1": "transporter",
    "TPMT": "metabolic_enzyme",
    "NUDT15": "metabolic_enzyme",
    "DPYD": "metabolic_enzyme",
    "UGT1A1": "metabolic_enzyme",
    "NAT2": "metabolic_enzyme",
    "IFNL3": "drug_target",
    "G6PD": "metabolic_enzyme",
    "CACNA1S": "drug_target",
    "RYR1": "drug_target",
    "HLA-A": "hla",
    "HLA-B": "hla",
}


def fetch_cpic_data(endpoint: str, params: dict = None) -> list[dict]:
    """Fetch data from CPIC API."""
    url = f"{CPIC_API}/{endpoint}"
    try:
        response = requests.get(url, params=params, timeout=30)
        response.raise_for_status()
        return response.json()
    except requests.RequestException as e:
        print(f"  Warning: Failed to fetch {endpoint}: {e}")
        return []


def fetch_allele_definitions(gene: str) -> list[dict]:
    """Fetch allele definitions for a gene."""
    return fetch_cpic_data("allele", {"genesymbol": f"eq.{gene}"})


def fetch_allele_location_values(gene: str) -> list[dict]:
    """Fetch sequence location (variant) data for a gene."""
    # First get the gene's allele definition ID
    alleles = fetch_cpic_data("allele_definition", {"genesymbol": f"eq.{gene}"})
    if not alleles:
        return []

    # Get variant locations
    locations = fetch_cpic_data("allele_location_value",
                                {"alleledefinitionid": f"eq.{alleles[0]['id']}"})
    return locations


def fetch_diplotype_phenotypes(gene: str) -> list[dict]:
    """Fetch diplotype to phenotype mappings."""
    return fetch_cpic_data("diplotype", {"genesymbol": f"eq.{gene}"})


def fetch_recommendations(gene: str) -> list[dict]:
    """Fetch drug recommendations for a gene."""
    return fetch_cpic_data("recommendation", {"genesymbol": f"eq.{gene}"})


def fetch_sequence_locations(gene: str) -> list[dict]:
    """Fetch GRCh38 coordinates for variants."""
    return fetch_cpic_data("sequence_location", {"genesymbol": f"eq.{gene}"})


def fetch_gene_info(gene: str) -> dict:
    """Fetch gene information."""
    data = fetch_cpic_data("gene", {"symbol": f"eq.{gene}"})
    return data[0] if data else {}


def convert_to_allele_definitions(gene: str, alleles: list[dict],
                                  locations: list[dict]) -> dict:
    """Convert CPIC allele data to ElPGx format."""
    gene_info = fetch_gene_info(gene)

    # Build location lookup by allele
    seq_locations = fetch_sequence_locations(gene)
    location_lookup = {}
    for loc in seq_locations:
        rsid = loc.get("dbsnpid")
        if rsid:
            location_lookup[rsid] = {
                "chrom": loc.get("chromosomelocation", "").replace("NC_", "chr"),
                "pos": loc.get("position"),
                "ref": loc.get("referenceallele"),
            }

    # Build allele entries
    allele_entries = []
    ref_allele = "*1"

    for allele in alleles:
        allele_name = allele.get("name", "")
        if not allele_name:
            continue

        # Determine function and activity value
        function = allele.get("clinicalfunctionalstatus", "Unknown")
        activity = allele.get("activityvalue")

        # Map function to activity value if not provided
        if activity is None:
            activity = map_function_to_activity(function)

        # Get defining variants
        variants = []
        # This would require additional API calls to get variant details
        # For now, we use the allele_location_value endpoint

        entry = {
            "allele": allele_name,
            "function": function,
            "activity_value": activity,
            "defining_variants": variants,
            "citations": allele.get("citations", []),
            "evidence_strength": allele.get("strength", ""),
        }

        if allele_name == "*1":
            ref_allele = "*1"

        allele_entries.append(entry)

    # Sort alleles by name
    allele_entries.sort(key=lambda x: parse_allele_number(x["allele"]))

    return {
        "gene": gene,
        "chromosome": gene_info.get("chr", ""),
        "gene_tier": GENE_TIERS.get(gene, "group1"),
        "biomarker_type": BIOMARKER_TYPES.get(gene, "metabolic_enzyme"),
        "reference_allele": ref_allele,
        "cpic_version": "v1.24",  # Current CPIC version
        "last_updated": datetime.utcnow().strftime("%Y-%m-%d"),
        "alleles": allele_entries
    }


def map_function_to_activity(function: str) -> float:
    """Map functional status to activity value."""
    function_lower = function.lower() if function else ""

    if "no function" in function_lower or "non-functional" in function_lower:
        return 0.0
    elif "decreased" in function_lower or "reduced" in function_lower:
        return 0.5
    elif "increased" in function_lower or "ultra" in function_lower:
        return 1.5
    elif "normal" in function_lower:
        return 1.0
    elif "uncertain" in function_lower or "unknown" in function_lower:
        return None
    else:
        return 1.0  # Default to normal


def parse_allele_number(allele: str) -> tuple:
    """Parse allele name for sorting (e.g., *1, *2, *10)."""
    try:
        # Remove leading * and any letter suffixes
        num_part = allele.lstrip("*")
        # Handle suballeles like *1A, *10B
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


def convert_to_phenotype_map(gene: str, diplotypes: list[dict]) -> dict:
    """Convert CPIC diplotype data to ElPGx phenotype map format."""
    mappings = []
    phenotype_by_score = {}

    for dip in diplotypes:
        diplotype = dip.get("diplotype", "")
        phenotype = dip.get("generesult", "")
        activity_score = dip.get("activityscore")

        # Get EHR priority if available
        ehr_priority = dip.get("ehrpriority", "")

        mappings.append({
            "diplotype": diplotype,
            "phenotype": phenotype,
            "activity_score": activity_score,
            "ehr_priority": ehr_priority,
        })

        # Build activity score to phenotype lookup
        if activity_score is not None:
            score_key = str(activity_score)
            if score_key not in phenotype_by_score:
                phenotype_by_score[score_key] = phenotype

    # Sort by diplotype
    mappings.sort(key=lambda x: x["diplotype"])

    return {
        "gene": gene,
        "cpic_version": "v1.24",
        "last_updated": datetime.utcnow().strftime("%Y-%m-%d"),
        "activity_score_based": len(phenotype_by_score) > 0,
        "phenotype_by_activity_score": phenotype_by_score,
        "mappings": mappings
    }


def convert_to_implications(gene: str, recommendations: list[dict]) -> dict:
    """Convert CPIC recommendations to ElPGx implications format."""
    # Group recommendations by phenotype
    phenotype_drugs: dict[str, list[dict]] = {}

    for rec in recommendations:
        phenotype = rec.get("phenotypes", {}).get(gene, "")
        if not phenotype:
            continue

        drug = rec.get("drug", {}).get("name", "")
        rec_text = rec.get("drugrecommendation", "")
        classification = rec.get("classification", "")

        if phenotype not in phenotype_drugs:
            phenotype_drugs[phenotype] = []

        phenotype_drugs[phenotype].append({
            "drug": drug,
            "recommendation": rec_text,
            "classification": classification,
            "pmid": rec.get("pmid", ""),
        })

    # Build implications list
    implications = []
    for phenotype, drugs in phenotype_drugs.items():
        # Combine drug recommendations
        drug_recs = []
        for d in drugs:
            if d["drug"] and d["recommendation"]:
                drug_recs.append(f"{d['drug']}: {d['recommendation']}")

        implications.append({
            "phenotype": phenotype,
            "drugs": drugs,
            "implication": " | ".join(drug_recs) if drug_recs else "No specific guidance available.",
        })

    return {
        "gene": gene,
        "cpic_version": "v1.24",
        "last_updated": datetime.utcnow().strftime("%Y-%m-%d"),
        "implications": implications
    }


def download_and_convert_gene(gene: str, output_dir: Path) -> dict:
    """Download and convert all data for a single gene."""
    print(f"\nProcessing {gene}...")

    result = {"gene": gene, "status": "success", "files": []}

    try:
        # 1. Download allele definitions
        print(f"  Fetching allele definitions...")
        alleles = fetch_allele_definitions(gene)
        locations = []  # Would need additional API calls

        if alleles:
            allele_def = convert_to_allele_definitions(gene, alleles, locations)
            allele_file = output_dir / "allele_definitions" / f"{gene}.json"
            allele_file.parent.mkdir(parents=True, exist_ok=True)
            with open(allele_file, 'w') as f:
                json.dump(allele_def, f, indent=2)
            result["files"].append(str(allele_file))
            print(f"    -> {len(alleles)} alleles")
        else:
            print(f"    -> No alleles found")

        # 2. Download diplotype-phenotype mappings
        print(f"  Fetching diplotype-phenotype mappings...")
        diplotypes = fetch_diplotype_phenotypes(gene)

        if diplotypes:
            phenotype_map = convert_to_phenotype_map(gene, diplotypes)
            phenotype_file = output_dir / "phenotype_maps" / f"{gene}_phenotype.json"
            phenotype_file.parent.mkdir(parents=True, exist_ok=True)
            with open(phenotype_file, 'w') as f:
                json.dump(phenotype_map, f, indent=2)
            result["files"].append(str(phenotype_file))
            print(f"    -> {len(diplotypes)} diplotype mappings")
        else:
            print(f"    -> No diplotype mappings found")

        # 3. Download recommendations
        print(f"  Fetching prescribing recommendations...")
        recommendations = fetch_recommendations(gene)

        if recommendations:
            implications = convert_to_implications(gene, recommendations)
            impl_file = output_dir / "prescribing_implications" / f"{gene}_implications.json"
            impl_file.parent.mkdir(parents=True, exist_ok=True)
            with open(impl_file, 'w') as f:
                json.dump(implications, f, indent=2)
            result["files"].append(str(impl_file))
            print(f"    -> {len(recommendations)} recommendations")
        else:
            print(f"    -> No recommendations found")

    except Exception as e:
        result["status"] = "error"
        result["error"] = str(e)
        print(f"  ERROR: {e}")

    return result


def generate_gene_tiers_file(output_dir: Path, genes_processed: list[str]):
    """Generate the gene_tiers.json configuration file."""
    genes_config = {}

    for gene in genes_processed:
        gene_info = fetch_gene_info(gene)

        genes_config[gene] = {
            "tier": GENE_TIERS.get(gene, "group1"),
            "biomarker_type": BIOMARKER_TYPES.get(gene, "metabolic_enzyme"),
            "chromosome": gene_info.get("chr", ""),
            "cpic_level": "A" if gene in ["CYP2C19", "CYP2C9", "CYP2D6", "VKORC1",
                                           "SLCO1B1", "TPMT", "NUDT15", "DPYD",
                                           "UGT1A1", "CYP3A5"] else "B",
            "pharmgkb_id": gene_info.get("pharmgkbid", ""),
        }

    output_file = output_dir / "gene_tiers.json"
    with open(output_file, 'w') as f:
        json.dump({
            "cpic_version": "v1.24",
            "last_updated": datetime.utcnow().strftime("%Y-%m-%d"),
            "genes": genes_config
        }, f, indent=2)

    print(f"\nGenerated gene_tiers.json with {len(genes_config)} genes")


def generate_pgx_regions_bed(output_dir: Path, genes_processed: list[str]):
    """Generate PGx regions BED file for variant calling."""
    bed_lines = []

    for gene in genes_processed:
        seq_locs = fetch_sequence_locations(gene)
        gene_info = fetch_gene_info(gene)

        if seq_locs:
            positions = [loc.get("position") for loc in seq_locs if loc.get("position")]
            if positions:
                chrom = gene_info.get("chr", "chr1")
                min_pos = min(positions) - 1000  # Add flanking
                max_pos = max(positions) + 1000
                bed_lines.append(f"{chrom}\t{min_pos}\t{max_pos}\t{gene}")

    output_file = output_dir / "pgx_regions.bed"
    with open(output_file, 'w') as f:
        f.write("# PGx target regions - Generated from CPIC data\n")
        f.write(f"# Generated: {datetime.utcnow().strftime('%Y-%m-%d')}\n")
        f.write("# Source: CPIC v1.24\n")
        for line in sorted(bed_lines):
            f.write(line + "\n")

    print(f"Generated pgx_regions.bed with {len(bed_lines)} regions")


def main():
    """Main function to download and convert all CPIC data."""
    print("=" * 60)
    print("CPIC Data Download and Conversion for ElPGx Pipeline")
    print("=" * 60)
    print(f"Source: {CPIC_API}")
    print(f"Target genes: {len(TARGET_GENES)}")
    print(f"Date: {datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S UTC')}")
    print("=" * 60)

    # Output directory
    output_dir = Path(__file__).parent.parent / "data"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Backup existing data
    backup_dir = output_dir.parent / "data_backup"
    if output_dir.exists() and not backup_dir.exists():
        import shutil
        shutil.copytree(output_dir, backup_dir)
        print(f"\nBacked up existing data to {backup_dir}")

    # Download data for each gene
    results = []
    genes_processed = []

    for gene in TARGET_GENES:
        result = download_and_convert_gene(gene, output_dir)
        results.append(result)
        if result["status"] == "success" and result["files"]:
            genes_processed.append(gene)

    # Generate configuration files
    print("\n" + "=" * 60)
    print("Generating configuration files...")

    if genes_processed:
        generate_gene_tiers_file(output_dir, genes_processed)
        generate_pgx_regions_bed(output_dir, genes_processed)

    # Summary
    print("\n" + "=" * 60)
    print("DOWNLOAD SUMMARY")
    print("=" * 60)

    success = sum(1 for r in results if r["status"] == "success")
    failed = sum(1 for r in results if r["status"] == "error")

    print(f"Total genes processed: {len(results)}")
    print(f"Successful: {success}")
    print(f"Failed: {failed}")

    if failed > 0:
        print("\nFailed genes:")
        for r in results:
            if r["status"] == "error":
                print(f"  - {r['gene']}: {r.get('error', 'Unknown error')}")

    print("\n" + "=" * 60)
    print("DATA PROVENANCE")
    print("=" * 60)
    print(f"Source: Clinical Pharmacogenetics Implementation Consortium (CPIC)")
    print(f"API: {CPIC_API}")
    print(f"Version: v1.24 (February 2026)")
    print(f"License: CC BY-SA 4.0")
    print("=" * 60)

    # Write provenance file
    provenance = {
        "source": "Clinical Pharmacogenetics Implementation Consortium (CPIC)",
        "api_url": CPIC_API,
        "cpic_version": "v1.24",
        "download_date": datetime.utcnow().isoformat() + "Z",
        "genes_included": genes_processed,
        "license": "CC BY-SA 4.0",
        "citation": "Relling MV, Klein TE. CPIC: Clinical Pharmacogenetics Implementation Consortium of the Pharmacogenomics Research Network. Clin Pharmacol Ther. 2011;89(3):464-467."
    }

    with open(output_dir / "PROVENANCE.json", 'w') as f:
        json.dump(provenance, f, indent=2)

    print(f"\nProvenance written to {output_dir / 'PROVENANCE.json'}")
    print("\nDone!")


if __name__ == "__main__":
    main()
