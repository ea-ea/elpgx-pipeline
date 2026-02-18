"""Unit tests for pgx.vcf_parser module."""

import json
from pathlib import Path

import pytest

from elpgx.vcf_parser import (
    ParsedVariant,
    filter_by_positions,
    has_variant,
    load_target_positions,
    to_variant_model,
    variants_at_position,
)


@pytest.fixture
def sample_parsed_variants():
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
        ParsedVariant(
            chrom="chr12",
            pos=21178615,
            ref="T",
            alt="C",
            rsid="rs4149056",
            genotype="0/1",
            quality=45.0,
            depth=40,
            genotype_quality=95,
        ),
    ]


@pytest.fixture
def target_positions():
    """Create sample target positions."""
    return {
        "chr10": {94761900, 94781859, 94762706},
        "chr12": {21178615},
    }


@pytest.fixture
def allele_defs_dir(tmp_path):
    """Create temporary allele definitions directory."""
    gene_data = {
        "gene": "CYP2C19",
        "chromosome": "chr10",
        "alleles": [
            {
                "allele": "*17",
                "defining_variants": [
                    {"chrom": "chr10", "pos": 94761900, "ref": "C", "alt": "T"},
                ],
            },
            {
                "allele": "*2",
                "defining_variants": [
                    {"chrom": "chr10", "pos": 94781859, "ref": "G", "alt": "A"},
                ],
            },
        ],
    }

    json_path = tmp_path / "CYP2C19.json"
    with open(json_path, "w") as f:
        json.dump(gene_data, f)

    return tmp_path


class TestParsedVariant:
    """Tests for ParsedVariant dataclass."""

    def test_create_parsed_variant(self):
        """Test creating a ParsedVariant."""
        pv = ParsedVariant(
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

        assert pv.chrom == "chr10"
        assert pv.pos == 94761900
        assert pv.ref == "C"
        assert pv.alt == "T"
        assert pv.rsid == "rs12248560"
        assert pv.genotype == "0/1"
        assert pv.quality == 38.5
        assert pv.depth == 42
        assert pv.genotype_quality == 99

    def test_create_without_rsid(self):
        """Test creating a ParsedVariant without rsid."""
        pv = ParsedVariant(
            chrom="chr10",
            pos=94761900,
            ref="C",
            alt="T",
            rsid=None,
            genotype="0/1",
            quality=38.5,
            depth=42,
            genotype_quality=99,
        )

        assert pv.rsid is None


class TestLoadTargetPositions:
    """Tests for load_target_positions function."""

    def test_load_from_allele_defs(self, allele_defs_dir):
        """Test loading target positions from allele definition files."""
        positions = load_target_positions(allele_defs_dir)

        assert "chr10" in positions
        assert 94761900 in positions["chr10"]
        assert 94781859 in positions["chr10"]

    def test_load_empty_dir(self, tmp_path):
        """Test loading from empty directory."""
        positions = load_target_positions(tmp_path)
        assert len(positions) == 0


class TestFilterByPositions:
    """Tests for filter_by_positions function."""

    def test_filter_matching(self, sample_parsed_variants, target_positions):
        """Test filtering variants to target positions."""
        filtered = filter_by_positions(sample_parsed_variants, target_positions)

        assert len(filtered) == 3  # All variants are at target positions

    def test_filter_partial(self, sample_parsed_variants):
        """Test filtering when only some variants match."""
        limited_targets = {"chr10": {94761900}}
        filtered = filter_by_positions(sample_parsed_variants, limited_targets)

        assert len(filtered) == 1
        assert filtered[0].pos == 94761900

    def test_filter_no_match(self, sample_parsed_variants):
        """Test filtering when no variants match."""
        no_targets = {"chr99": {1234}}
        filtered = filter_by_positions(sample_parsed_variants, no_targets)

        assert len(filtered) == 0

    def test_filter_empty_variants(self, target_positions):
        """Test filtering empty variant list."""
        filtered = filter_by_positions([], target_positions)
        assert len(filtered) == 0


class TestToVariantModel:
    """Tests for to_variant_model function."""

    def test_basic_conversion(self, sample_parsed_variants):
        """Test basic conversion to Variant model."""
        parsed = sample_parsed_variants[0]
        variant = to_variant_model(parsed)

        assert variant.chrom == parsed.chrom
        assert variant.pos == parsed.pos
        assert variant.ref == parsed.ref
        assert variant.alt == parsed.alt
        assert variant.rsid == parsed.rsid
        assert variant.genotype == parsed.genotype
        assert variant.quality == parsed.quality
        assert variant.depth == parsed.depth
        assert variant.allele_defining is False

    def test_conversion_with_star_allele(self, sample_parsed_variants):
        """Test conversion with allele-defining metadata."""
        parsed = sample_parsed_variants[0]
        variant = to_variant_model(
            parsed,
            allele_defining=True,
            star_allele="*17",
            functional_impact="increased function",
        )

        assert variant.allele_defining is True
        assert variant.star_allele == "*17"
        assert variant.functional_impact == "increased function"


class TestVariantsAtPosition:
    """Tests for variants_at_position function."""

    def test_find_single_variant(self, sample_parsed_variants):
        """Test finding variant at specific position."""
        results = variants_at_position(sample_parsed_variants, "chr10", 94761900)

        assert len(results) == 1
        assert results[0].pos == 94761900

    def test_find_no_variant(self, sample_parsed_variants):
        """Test when no variant at position."""
        results = variants_at_position(sample_parsed_variants, "chr10", 99999999)

        assert len(results) == 0

    def test_find_wrong_chromosome(self, sample_parsed_variants):
        """Test finding variant on wrong chromosome."""
        results = variants_at_position(sample_parsed_variants, "chr99", 94761900)

        assert len(results) == 0


class TestHasVariant:
    """Tests for has_variant function."""

    def test_has_heterozygous_variant(self, sample_parsed_variants):
        """Test finding heterozygous variant (0/1)."""
        result = has_variant(
            sample_parsed_variants,
            "chr10",
            94761900,
            "C",
            "T",
        )

        assert result is not None
        assert result.genotype == "0/1"

    def test_has_homozygous_variant(self, sample_parsed_variants):
        """Test finding homozygous variant (1/1)."""
        result = has_variant(
            sample_parsed_variants,
            "chr10",
            94781859,
            "G",
            "A",
        )

        assert result is not None
        assert result.genotype == "1/1"

    def test_no_variant_at_position(self, sample_parsed_variants):
        """Test when no variant at position."""
        result = has_variant(
            sample_parsed_variants,
            "chr10",
            99999999,
            "A",
            "T",
        )

        assert result is None

    def test_wrong_alt_allele(self, sample_parsed_variants):
        """Test when position matches but alt allele differs."""
        result = has_variant(
            sample_parsed_variants,
            "chr10",
            94761900,
            "C",
            "G",  # Wrong alt
        )

        assert result is None

    def test_reference_genotype_not_found(self):
        """Test that reference genotype (0/0) is not found."""
        ref_variant = [
            ParsedVariant(
                chrom="chr10",
                pos=100,
                ref="A",
                alt="T",
                rsid=None,
                genotype="0/0",  # Reference genotype
                quality=30.0,
                depth=30,
                genotype_quality=99,
            )
        ]

        result = has_variant(ref_variant, "chr10", 100, "A", "T")
        assert result is None  # Should not find 0/0 genotype


class TestGenotypeFormats:
    """Tests for various genotype format handling."""

    def test_phased_genotype(self):
        """Test variant with phased genotype (|)."""
        pv = ParsedVariant(
            chrom="chr10",
            pos=100,
            ref="A",
            alt="T",
            rsid=None,
            genotype="0|1",  # Phased
            quality=30.0,
            depth=30,
            genotype_quality=99,
        )

        result = has_variant([pv], "chr10", 100, "A", "T")
        assert result is not None

    def test_missing_genotype(self):
        """Test variant with missing genotype (./.)."""
        pv = ParsedVariant(
            chrom="chr10",
            pos=100,
            ref="A",
            alt="T",
            rsid=None,
            genotype="./.",
            quality=0.0,
            depth=0,
            genotype_quality=0,
        )

        result = has_variant([pv], "chr10", 100, "A", "T")
        assert result is None  # Missing genotype should not match


class TestMultiAllelicHandling:
    """Tests for multi-allelic site handling."""

    def test_multi_allelic_variants(self):
        """Test handling of multi-allelic sites split into multiple records."""
        variants = [
            ParsedVariant(
                chrom="chr10",
                pos=100,
                ref="A",
                alt="T",
                rsid="rs123",
                genotype="0/1",
                quality=30.0,
                depth=30,
                genotype_quality=99,
            ),
            ParsedVariant(
                chrom="chr10",
                pos=100,
                ref="A",
                alt="C",
                rsid="rs123",
                genotype="0/2",
                quality=30.0,
                depth=30,
                genotype_quality=99,
            ),
        ]

        # Should find T allele
        result_t = has_variant(variants, "chr10", 100, "A", "T")
        assert result_t is not None
        assert result_t.alt == "T"

        # Should find C allele (genotype has 2 which counts as alt present)
        result_c = has_variant(variants, "chr10", 100, "A", "C")
        assert result_c is not None
        assert result_c.alt == "C"
