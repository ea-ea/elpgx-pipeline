"""Shared fixtures for unit tests."""

import json
from pathlib import Path

import pytest

from elpgx.allele_caller import DiplotypeCall
from elpgx.models import (
    AlleleDefinition,
    AlleleVariant,
    BiomarkerType,
    ConfidenceMetrics,
    GeneAlleleDefinitions,
    GeneTier,
    Variant,
)
from elpgx.phenotype import PhenotypeResult
from elpgx.vcf_parser import ParsedVariant


@pytest.fixture
def cyp2c19_gene_definition():
    """Create CYP2C19 gene definition for testing."""
    return GeneAlleleDefinitions(
        gene="CYP2C19",
        chromosome="chr10",
        gene_tier=GeneTier.GROUP1,
        biomarker_type=BiomarkerType.METABOLIC_ENZYME,
        reference_allele="*1",
        alleles=[
            AlleleDefinition(
                allele="*1",
                function="Normal Function",
                activity_value=1.0,
                defining_variants=[],
            ),
            AlleleDefinition(
                allele="*2",
                function="No Function",
                activity_value=0.0,
                defining_variants=[
                    AlleleVariant(
                        rsid="rs4244285",
                        chrom="chr10",
                        pos=94781859,
                        ref="G",
                        alt="A",
                        functional_impact="no function",
                    )
                ],
            ),
            AlleleDefinition(
                allele="*17",
                function="Increased Function",
                activity_value=1.5,
                defining_variants=[
                    AlleleVariant(
                        rsid="rs12248560",
                        chrom="chr10",
                        pos=94761900,
                        ref="C",
                        alt="T",
                        functional_impact="increased function",
                    )
                ],
            ),
        ],
    )


@pytest.fixture
def cyp2c19_star17_variant():
    """Create *17 defining variant."""
    return ParsedVariant(
        chrom="chr10",
        pos=94761900,
        ref="C",
        alt="T",
        rsid="rs12248560",
        genotype="0/1",
        quality=38.5,
        depth=42,
        genotype_quality=99,
    )


@pytest.fixture
def cyp2c19_star2_variant_hom():
    """Create homozygous *2 defining variant."""
    return ParsedVariant(
        chrom="chr10",
        pos=94781859,
        ref="G",
        alt="A",
        rsid="rs4244285",
        genotype="1/1",
        quality=50.0,
        depth=55,
        genotype_quality=99,
    )


@pytest.fixture
def sample_diplotype_call():
    """Create sample diplotype call."""
    return DiplotypeCall(
        gene="CYP2C19",
        allele1="*1",
        allele2="*17",
        diplotype="*1/*17",
        variants=[
            Variant(
                rsid="rs12248560",
                chrom="chr10",
                pos=94761900,
                ref="C",
                alt="T",
                genotype="0/1",
                quality=38.5,
                depth=42,
                allele_defining=True,
                star_allele="*17",
            )
        ],
        allele1_function="Normal Function",
        allele2_function="Increased Function",
        allele1_activity=1.0,
        allele2_activity=1.5,
        call_quality="definitive",
    )


@pytest.fixture
def sample_phenotype_result():
    """Create sample phenotype result."""
    return PhenotypeResult(
        gene="CYP2C19",
        diplotype="*1/*17",
        activity_score=2.5,
        functional_status="Normal Function / Increased Function",
        phenotype="Rapid Metabolizer",
        prescribing_implication="Consider alternative to clopidogrel.",
    )


@pytest.fixture
def high_quality_variants():
    """Create high-quality variant list."""
    return [
        Variant(
            chrom="chr10",
            pos=94762706,
            ref="C",
            alt="T",
            genotype="0/1",
            quality=50.0,
            depth=45,
        ),
        Variant(
            chrom="chr10",
            pos=94781859,
            ref="G",
            alt="A",
            genotype="0/1",
            quality=45.0,
            depth=40,
        ),
    ]


@pytest.fixture
def low_quality_variants():
    """Create low-quality variant list."""
    return [
        Variant(
            chrom="chr10",
            pos=94762706,
            ref="C",
            alt="T",
            genotype="0/1",
            quality=15.0,
            depth=5,
        ),
    ]


@pytest.fixture
def sample_confidence():
    """Create sample confidence metrics."""
    return ConfidenceMetrics(score=0.95, flags=[])


@pytest.fixture
def test_data_dir(tmp_path):
    """Create temporary data directory with all required files."""
    # Create gene_tiers.json
    gene_tiers = {
        "genes": {
            "CYP2C19": {
                "tier": "group1",
                "biomarker_type": "metabolic_enzyme",
                "chromosome": "chr10",
                "cpic_level": "A",
            },
            "CYP2C9": {
                "tier": "group1",
                "biomarker_type": "metabolic_enzyme",
                "chromosome": "chr10",
                "cpic_level": "A",
            },
            "SLCO1B1": {
                "tier": "group1",
                "biomarker_type": "transporter",
                "chromosome": "chr12",
                "cpic_level": "A",
            },
        }
    }
    tiers_file = tmp_path / "gene_tiers.json"
    with open(tiers_file, "w") as f:
        json.dump(gene_tiers, f)

    # Create phenotype map directory
    phenotype_maps_dir = tmp_path / "phenotype_maps"
    phenotype_maps_dir.mkdir()

    cyp2c19_phenotype = {
        "gene": "CYP2C19",
        "activity_score_based": True,
        "phenotype_by_activity_score": {
            "0": "Poor Metabolizer",
            "1": "Intermediate Metabolizer",
            "2": "Normal Metabolizer",
            "2.5": "Rapid Metabolizer",
            "3": "Ultrarapid Metabolizer",
        },
        "mappings": [
            {
                "diplotype": "*1/*1",
                "phenotype": "Normal Metabolizer",
                "activity_score": 2.0,
            },
            {
                "diplotype": "*1/*17",
                "phenotype": "Rapid Metabolizer",
                "activity_score": 2.5,
            },
        ],
    }
    with open(phenotype_maps_dir / "CYP2C19_phenotype.json", "w") as f:
        json.dump(cyp2c19_phenotype, f)

    # Create prescribing implications directory
    implications_dir = tmp_path / "prescribing_implications"
    implications_dir.mkdir()

    cyp2c19_implications = {
        "gene": "CYP2C19",
        "implications": [
            {
                "phenotype": "Normal Metabolizer",
                "implication": "Standard dosing appropriate.",
            },
            {
                "phenotype": "Rapid Metabolizer",
                "implication": "Consider alternative to clopidogrel.",
            },
        ],
    }
    with open(implications_dir / "CYP2C19_implications.json", "w") as f:
        json.dump(cyp2c19_implications, f)

    # Create allele definitions directory
    allele_defs_dir = tmp_path / "allele_definitions"
    allele_defs_dir.mkdir()

    cyp2c19_alleles = {
        "gene": "CYP2C19",
        "chromosome": "chr10",
        "gene_tier": "group1",
        "biomarker_type": "metabolic_enzyme",
        "reference_allele": "*1",
        "alleles": [
            {
                "allele": "*1",
                "function": "Normal Function",
                "activity_value": 1.0,
                "defining_variants": [],
            },
            {
                "allele": "*17",
                "function": "Increased Function",
                "activity_value": 1.5,
                "defining_variants": [
                    {
                        "rsid": "rs12248560",
                        "chrom": "chr10",
                        "pos": 94761900,
                        "ref": "C",
                        "alt": "T",
                        "functional_impact": "increased function",
                    }
                ],
            },
        ],
    }
    with open(allele_defs_dir / "CYP2C19.json", "w") as f:
        json.dump(cyp2c19_alleles, f)

    return tmp_path


@pytest.fixture
def sample_bed_file(tmp_path):
    """Create sample BED file with PGx regions."""
    bed_content = """# PGx target regions for Group 1 genes
chr10\t94760600\t94853000\tCYP2C19
chr10\t94938600\t94989500\tCYP2C9
chr12\t21130000\t21240000\tSLCO1B1
chr16\t88922500\t89004500\tVKORC1
chr6\t18128000\t18155000\tTPMT
chr1\t97078000\t97098000\tDPYD
chr2\t233757000\t233773000\tUGT1A1
chr13\t48037000\t48052000\tNUDT15
"""
    bed_path = tmp_path / "pgx_regions.bed"
    bed_path.write_text(bed_content)
    return bed_path
