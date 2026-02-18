"""Unit tests for pgx.allele_caller module."""

import json
from pathlib import Path

import pytest

from elpgx.allele_caller import (
    AlleleMatch,
    DiplotypeCall,
    count_allele_copies,
    determine_diplotype,
    find_matching_alleles,
    load_gene_definitions,
    match_allele,
)
from elpgx.models import (
    AlleleDefinition,
    AlleleVariant,
    BiomarkerType,
    GeneAlleleDefinitions,
    GeneTier,
)
from elpgx.vcf_parser import ParsedVariant


@pytest.fixture
def sample_gene_definition():
    """Create a sample gene definition for testing."""
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
def sample_variants():
    """Create sample parsed variants for testing."""
    return [
        ParsedVariant(
            chrom="chr10",
            pos=94761900,
            ref="C",
            alt="T",
            rsid="rs12248560",
            genotype="0/1",
            quality=38.5,
            depth=42,
            genotype_quality=99,
        ),
    ]


@pytest.fixture
def homozygous_variant():
    """Create a homozygous variant."""
    return [
        ParsedVariant(
            chrom="chr10",
            pos=94781859,
            ref="G",
            alt="A",
            rsid="rs4244285",
            genotype="1/1",
            quality=50.0,
            depth=55,
            genotype_quality=99,
        ),
    ]


class TestMatchAllele:
    """Tests for match_allele function."""

    def test_match_allele_with_variant(self, sample_gene_definition, sample_variants):
        """Test matching when variant is present."""
        star17_def = sample_gene_definition.alleles[2]  # *17
        match = match_allele(star17_def, sample_variants)

        assert match.allele == "*17"
        assert len(match.matched_variants) == 1
        assert match.is_complete_match is True

    def test_match_allele_no_variant(self, sample_gene_definition):
        """Test matching when no variants present."""
        star2_def = sample_gene_definition.alleles[1]  # *2
        match = match_allele(star2_def, [])

        assert match.allele == "*2"
        assert len(match.matched_variants) == 0
        assert match.is_complete_match is False

    def test_match_reference_allele(self, sample_gene_definition):
        """Test matching reference allele (no defining variants)."""
        star1_def = sample_gene_definition.alleles[0]  # *1
        match = match_allele(star1_def, [])

        assert match.allele == "*1"
        assert len(match.matched_variants) == 0
        # Reference allele has no defining variants, so not complete
        assert match.is_complete_match is False


class TestFindMatchingAlleles:
    """Tests for find_matching_alleles function."""

    def test_find_single_match(self, sample_gene_definition, sample_variants):
        """Test finding a single matching allele."""
        matches = find_matching_alleles(sample_gene_definition, sample_variants)

        assert len(matches) == 1
        assert matches[0].allele == "*17"

    def test_find_no_matches(self, sample_gene_definition):
        """Test when no alleles match."""
        matches = find_matching_alleles(sample_gene_definition, [])
        assert len(matches) == 0


class TestCountAlleleCopies:
    """Tests for count_allele_copies function."""

    def test_heterozygous_single_copy(self, sample_variants):
        """Test counting copies for heterozygous genotype."""
        match = AlleleMatch(
            allele="*17",
            function="Increased Function",
            activity_value=1.5,
            matched_variants=sample_variants,
            total_defining_variants=1,
            is_complete_match=True,
        )
        copies = count_allele_copies(match, sample_variants)
        assert copies == 1

    def test_homozygous_two_copies(self, homozygous_variant):
        """Test counting copies for homozygous genotype."""
        match = AlleleMatch(
            allele="*2",
            function="No Function",
            activity_value=0.0,
            matched_variants=homozygous_variant,
            total_defining_variants=1,
            is_complete_match=True,
        )
        copies = count_allele_copies(match, homozygous_variant)
        assert copies == 2

    def test_no_match_zero_copies(self):
        """Test counting copies when no match."""
        match = AlleleMatch(
            allele="*2",
            function="No Function",
            activity_value=0.0,
            matched_variants=[],
            total_defining_variants=1,
            is_complete_match=False,
        )
        copies = count_allele_copies(match, [])
        assert copies == 0


class TestDetermineDiplotype:
    """Tests for determine_diplotype function."""

    def test_diplotype_heterozygous(self, sample_gene_definition, sample_variants):
        """Test diplotype determination for heterozygous *1/*17."""
        call = determine_diplotype(sample_gene_definition, sample_variants)

        assert call.gene == "CYP2C19"
        assert "*1" in call.diplotype
        assert "*17" in call.diplotype
        assert call.call_quality == "definitive"

    def test_diplotype_homozygous_reference(self, sample_gene_definition):
        """Test diplotype determination for homozygous reference *1/*1."""
        call = determine_diplotype(sample_gene_definition, [])

        assert call.gene == "CYP2C19"
        assert call.diplotype == "*1/*1"
        assert call.allele1 == "*1"
        assert call.allele2 == "*1"

    def test_diplotype_homozygous_variant(self, sample_gene_definition, homozygous_variant):
        """Test diplotype determination for homozygous *2/*2."""
        call = determine_diplotype(sample_gene_definition, homozygous_variant)

        assert call.gene == "CYP2C19"
        assert call.diplotype == "*2/*2"
        assert call.allele1 == "*2"
        assert call.allele2 == "*2"


class TestLoadGeneDefinitions:
    """Tests for load_gene_definitions function."""

    def test_load_definitions(self, tmp_path):
        """Test loading gene definitions from JSON files."""
        # Create a test JSON file
        gene_data = {
            "gene": "TEST_GENE",
            "chromosome": "chr1",
            "gene_tier": "group1",
            "biomarker_type": "metabolic_enzyme",
            "reference_allele": "*1",
            "alleles": [
                {
                    "allele": "*1",
                    "function": "Normal Function",
                    "activity_value": 1.0,
                    "defining_variants": [],
                }
            ],
        }

        json_path = tmp_path / "TEST_GENE.json"
        with open(json_path, "w") as f:
            json.dump(gene_data, f)

        definitions = load_gene_definitions(tmp_path)

        assert "TEST_GENE" in definitions
        assert definitions["TEST_GENE"].gene == "TEST_GENE"
        assert definitions["TEST_GENE"].reference_allele == "*1"

    def test_load_definitions_empty_dir(self, tmp_path):
        """Test loading from empty directory."""
        definitions = load_gene_definitions(tmp_path)
        assert len(definitions) == 0


class TestDiplotypeCall:
    """Tests for DiplotypeCall dataclass."""

    def test_diplotype_call_creation(self):
        """Test creating a DiplotypeCall."""
        call = DiplotypeCall(
            gene="CYP2C19",
            allele1="*1",
            allele2="*17",
            diplotype="*1/*17",
            variants=[],
            allele1_function="Normal Function",
            allele2_function="Increased Function",
            allele1_activity=1.0,
            allele2_activity=1.5,
            call_quality="definitive",
        )

        assert call.gene == "CYP2C19"
        assert call.diplotype == "*1/*17"
        assert call.allele1_activity + call.allele2_activity == 2.5
