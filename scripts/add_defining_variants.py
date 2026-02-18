#!/usr/bin/env python3
"""
Add defining variants to allele definitions by fetching allele_location_value data.

This script enhances the downloaded CPIC data by adding the defining variant
positions for each star allele.
"""

import json
import requests
from pathlib import Path
from datetime import datetime, timezone
import time

CPIC_API = "https://api.cpicpgx.org/v1"

# Key allele-defining variants from CPIC (rsID -> allele mappings for major genes)
# These are the most clinically important variants from CPIC guidelines
KEY_VARIANTS = {
    "CYP2C19": {
        "rs4244285": {"allele": "*2", "ref": "G", "alt": "A", "impact": "no function (splicing defect)"},
        "rs4986893": {"allele": "*3", "ref": "G", "alt": "A", "impact": "no function (premature stop)"},
        "rs28399504": {"allele": "*4", "ref": "A", "alt": "G", "impact": "no function"},
        "rs12248560": {"allele": "*17", "ref": "C", "alt": "T", "impact": "increased function (promoter)"},
        "rs56337013": {"allele": "*5", "ref": "C", "alt": "T", "impact": "no function"},
        "rs72552267": {"allele": "*6", "ref": "G", "alt": "A", "impact": "no function"},
        "rs41291556": {"allele": "*8", "ref": "T", "alt": "C", "impact": "no function"},
        "rs17884712": {"allele": "*9", "ref": "G", "alt": "A", "impact": "decreased function"},
    },
    "CYP2C9": {
        "rs1799853": {"allele": "*2", "ref": "C", "alt": "T", "impact": "decreased function"},
        "rs1057910": {"allele": "*3", "ref": "A", "alt": "C", "impact": "decreased function"},
        "rs28371686": {"allele": "*5", "ref": "C", "alt": "G", "impact": "no function"},
        "rs9332131": {"allele": "*6", "ref": "A", "alt": "-", "impact": "no function (frameshift)"},
        "rs28371685": {"allele": "*8", "ref": "C", "alt": "T", "impact": "decreased function"},
        "rs56165452": {"allele": "*11", "ref": "C", "alt": "T", "impact": "decreased function"},
    },
    "CYP2D6": {
        "rs16947": {"allele": "*2", "ref": "G", "alt": "A", "impact": "normal function"},
        "rs35742686": {"allele": "*3", "ref": "T", "alt": "-", "impact": "no function (frameshift)"},
        "rs3892097": {"allele": "*4", "ref": "G", "alt": "A", "impact": "no function (splicing)"},
        "rs5030655": {"allele": "*6", "ref": "T", "alt": "-", "impact": "no function (frameshift)"},
        "rs1065852": {"allele": "*10", "ref": "C", "alt": "T", "impact": "decreased function"},
        "rs28371706": {"allele": "*17", "ref": "C", "alt": "T", "impact": "decreased function"},
        "rs28371725": {"allele": "*41", "ref": "C", "alt": "T", "impact": "decreased function"},
    },
    "VKORC1": {
        "rs9923231": {"allele": "-1639G>A", "ref": "G", "alt": "A", "impact": "decreased expression"},
    },
    "SLCO1B1": {
        "rs4149056": {"allele": "*5", "ref": "T", "alt": "C", "impact": "decreased function"},
        "rs2306283": {"allele": "*1B", "ref": "A", "alt": "G", "impact": "normal function"},
        "rs4149015": {"allele": "*15", "ref": "G", "alt": "A", "impact": "decreased function"},
    },
    "TPMT": {
        "rs1800462": {"allele": "*2", "ref": "C", "alt": "G", "impact": "no function"},
        "rs1800460": {"allele": "*3B", "ref": "A", "alt": "G", "impact": "no function"},
        "rs1142345": {"allele": "*3C", "ref": "T", "alt": "C", "impact": "no function"},
    },
    "NUDT15": {
        "rs116855232": {"allele": "*3", "ref": "C", "alt": "T", "impact": "no function"},
        "rs147390019": {"allele": "*2", "ref": "G", "alt": "A", "impact": "no function"},
        "rs186364861": {"allele": "*4", "ref": "G", "alt": "A", "impact": "decreased function"},
    },
    "DPYD": {
        "rs3918290": {"allele": "*2A", "ref": "C", "alt": "T", "impact": "no function (IVS14+1G>A)"},
        "rs55886062": {"allele": "*13", "ref": "A", "alt": "C", "impact": "no function"},
        "rs67376798": {"allele": "c.2846A>T", "ref": "A", "alt": "T", "impact": "decreased function"},
        "rs75017182": {"allele": "HapB3", "ref": "G", "alt": "C", "impact": "decreased function"},
    },
    "UGT1A1": {
        "rs8175347": {"allele": "*28", "ref": "TA6", "alt": "TA7", "impact": "decreased function"},
        "rs4148323": {"allele": "*6", "ref": "G", "alt": "A", "impact": "decreased function"},
        "rs35350960": {"allele": "*27", "ref": "C", "alt": "A", "impact": "decreased function"},
    },
    "NAT2": {
        "rs1801280": {"allele": "*5", "ref": "T", "alt": "C", "impact": "slow acetylator"},
        "rs1799930": {"allele": "*6", "ref": "G", "alt": "A", "impact": "slow acetylator"},
        "rs1799931": {"allele": "*7", "ref": "G", "alt": "A", "impact": "slow acetylator"},
        "rs1801279": {"allele": "*14", "ref": "G", "alt": "A", "impact": "slow acetylator"},
    },
    "CYP3A5": {
        "rs776746": {"allele": "*3", "ref": "C", "alt": "T", "impact": "no function (splice defect)"},
        "rs10264272": {"allele": "*6", "ref": "C", "alt": "T", "impact": "no function"},
        "rs41303343": {"allele": "*7", "ref": "A", "alt": "-", "impact": "no function (frameshift)"},
    },
    "IFNL3": {
        "rs12979860": {"allele": "rs12979860", "ref": "C", "alt": "T", "impact": "favorable response (CC)"},
    },
    "CYP2B6": {
        "rs3745274": {"allele": "*6", "ref": "G", "alt": "T", "impact": "decreased function"},
        "rs2279343": {"allele": "*4", "ref": "A", "alt": "G", "impact": "increased function"},
        "rs28399499": {"allele": "*18", "ref": "T", "alt": "C", "impact": "decreased function"},
    },
    "G6PD": {
        "rs1050828": {"allele": "A-", "ref": "C", "alt": "T", "impact": "deficient (Mediterranean)"},
        "rs1050829": {"allele": "A-", "ref": "G", "alt": "A", "impact": "deficient"},
        "rs5030868": {"allele": "Mediterranean", "ref": "C", "alt": "T", "impact": "deficient"},
    },
}


def fetch_cpic_data(endpoint: str, params: dict = None) -> list[dict]:
    """Fetch data from CPIC API."""
    url = f"{CPIC_API}/{endpoint}"
    try:
        response = requests.get(url, params=params, timeout=30)
        response.raise_for_status()
        return response.json()
    except:
        return []


def update_allele_file(gene: str, data_dir: Path):
    """Update allele definition file with defining variants."""
    allele_file = data_dir / "allele_definitions" / f"{gene}.json"
    if not allele_file.exists():
        print(f"  Skipping {gene}: file not found")
        return

    with open(allele_file) as f:
        data = json.load(f)

    # Get sequence locations for this gene
    seq_locs = data.get("sequence_locations", [])
    rsid_to_pos = {loc["rsid"]: loc for loc in seq_locs}

    # Get key variants for this gene
    gene_variants = KEY_VARIANTS.get(gene, {})

    # Update alleles with defining variants
    for allele in data.get("alleles", []):
        allele_name = allele["allele"]
        defining_vars = []

        # Check if this allele has key variants defined
        for rsid, var_info in gene_variants.items():
            if var_info["allele"] == allele_name or var_info["allele"] in allele_name:
                # Get position from sequence_locations
                if rsid in rsid_to_pos:
                    loc = rsid_to_pos[rsid]
                    defining_vars.append({
                        "rsid": rsid,
                        "chrom": data.get("chromosome", ""),
                        "pos": loc["pos"],
                        "ref": var_info.get("ref") or loc.get("ref"),
                        "alt": var_info.get("alt") or loc.get("alt"),
                        "functional_impact": var_info.get("impact", ""),
                    })
                else:
                    # Fetch position from API
                    locs = fetch_cpic_data("sequence_location",
                                          {"dbsnpid": f"eq.{rsid}"})
                    if locs:
                        loc = locs[0]
                        defining_vars.append({
                            "rsid": rsid,
                            "chrom": data.get("chromosome", ""),
                            "pos": loc.get("position"),
                            "ref": var_info.get("ref"),
                            "alt": var_info.get("alt"),
                            "functional_impact": var_info.get("impact", ""),
                        })

        if defining_vars:
            allele["defining_variants"] = defining_vars

    # Update metadata
    data["last_updated"] = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    data["defining_variants_source"] = "CPIC v1.24 key allele variants"

    with open(allele_file, 'w') as f:
        json.dump(data, f, indent=2)

    variants_added = sum(1 for a in data["alleles"] if a.get("defining_variants"))
    print(f"  Updated {gene}: {variants_added} alleles with defining variants")


def main():
    """Main function."""
    print("=" * 60)
    print("Adding Defining Variants to Allele Definitions")
    print("=" * 60)

    data_dir = Path(__file__).parent.parent / "data"

    for gene in KEY_VARIANTS.keys():
        update_allele_file(gene, data_dir)
        time.sleep(0.2)

    print("\nDone!")


if __name__ == "__main__":
    main()
