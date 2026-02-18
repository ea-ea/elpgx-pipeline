"""Unit tests for pgx.confidence module."""

import pytest

from elpgx.confidence import (
    ConfidenceConfig,
    calculate_completeness_score,
    calculate_confidence,
    calculate_confidence_for_reference_call,
    calculate_depth_score,
    calculate_quality_score,
    determine_flags,
)
from elpgx.models import ConfidenceFlag, Variant


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
def default_config():
    """Create default confidence configuration."""
    return ConfidenceConfig()


class TestCalculateDepthScore:
    """Tests for calculate_depth_score function."""

    def test_high_depth(self, high_quality_variants):
        """Test depth score with high coverage."""
        score, mean, min_d = calculate_depth_score(high_quality_variants, target_depth=30)

        assert score > 0.9  # Should be close to 1.0
        assert mean == 42.5  # (45 + 40) / 2
        assert min_d == 40

    def test_low_depth(self, low_quality_variants):
        """Test depth score with low coverage."""
        score, mean, min_d = calculate_depth_score(low_quality_variants, target_depth=30)

        assert score < 0.5  # Should be low
        assert min_d == 5

    def test_empty_variants(self):
        """Test depth score with no variants."""
        score, mean, min_d = calculate_depth_score([])

        assert score == 0.0
        assert mean == 0.0
        assert min_d == 0

    def test_exceeds_target(self):
        """Test depth score capped at 1.0 when exceeding target."""
        high_depth_variant = [
            Variant(
                chrom="chr10",
                pos=100,
                ref="A",
                alt="T",
                genotype="0/1",
                quality=50.0,
                depth=100,
            )
        ]
        score, mean, min_d = calculate_depth_score(high_depth_variant, target_depth=30)

        assert score == 1.0  # Capped at 1.0


class TestCalculateQualityScore:
    """Tests for calculate_quality_score function."""

    def test_high_quality(self, high_quality_variants):
        """Test quality score with high-quality variants."""
        score = calculate_quality_score(high_quality_variants, max_quality=100.0)

        assert score > 0.4  # (50 + 45) / 2 / 100 = 0.475

    def test_low_quality(self, low_quality_variants):
        """Test quality score with low-quality variants."""
        score = calculate_quality_score(low_quality_variants, max_quality=100.0)

        assert score < 0.2  # 15 / 100 = 0.15

    def test_empty_variants(self):
        """Test quality score with no variants."""
        score = calculate_quality_score([])

        assert score == 0.0


class TestCalculateCompletenessScore:
    """Tests for calculate_completeness_score function."""

    def test_complete(self):
        """Test completeness score when all positions genotyped."""
        score = calculate_completeness_score(genotyped_positions=2, total_defining_positions=2)
        assert score == 1.0

    def test_partial(self):
        """Test completeness score when partially genotyped."""
        score = calculate_completeness_score(genotyped_positions=1, total_defining_positions=2)
        assert score == 0.5

    def test_no_defining_positions(self):
        """Test completeness score when no defining positions (reference allele)."""
        score = calculate_completeness_score(genotyped_positions=0, total_defining_positions=0)
        assert score == 1.0  # Reference allele case


class TestDetermineFlags:
    """Tests for determine_flags function."""

    def test_no_flags(self, high_quality_variants, default_config):
        """Test no flags for high-quality data."""
        flags = determine_flags(
            variants=high_quality_variants,
            completeness_score=1.0,
            call_quality="definitive",
            config=default_config,
        )

        assert len(flags) == 0

    def test_low_depth_flag(self, low_quality_variants, default_config):
        """Test LOW_DEPTH flag for low coverage."""
        flags = determine_flags(
            variants=low_quality_variants,
            completeness_score=1.0,
            call_quality="definitive",
            config=default_config,
        )

        assert ConfidenceFlag.LOW_DEPTH in flags

    def test_low_quality_flag(self, low_quality_variants, default_config):
        """Test LOW_QUALITY flag for low QUAL."""
        flags = determine_flags(
            variants=low_quality_variants,
            completeness_score=1.0,
            call_quality="definitive",
            config=default_config,
        )

        assert ConfidenceFlag.LOW_QUALITY in flags

    def test_incomplete_allele_flag(self, high_quality_variants, default_config):
        """Test INCOMPLETE_ALLELE flag."""
        flags = determine_flags(
            variants=high_quality_variants,
            completeness_score=0.5,  # Incomplete
            call_quality="definitive",
            config=default_config,
        )

        assert ConfidenceFlag.INCOMPLETE_ALLELE in flags

    def test_ambiguous_call_flag(self, high_quality_variants, default_config):
        """Test AMBIGUOUS_CALL flag."""
        flags = determine_flags(
            variants=high_quality_variants,
            completeness_score=1.0,
            call_quality="ambiguous",
            config=default_config,
        )

        assert ConfidenceFlag.AMBIGUOUS_CALL in flags

    def test_no_call_flag(self, high_quality_variants, default_config):
        """Test NO_CALL flag."""
        flags = determine_flags(
            variants=high_quality_variants,
            completeness_score=1.0,
            call_quality="no_call",
            config=default_config,
        )

        assert ConfidenceFlag.NO_CALL in flags


class TestCalculateConfidence:
    """Tests for calculate_confidence function."""

    def test_high_confidence(self, high_quality_variants):
        """Test high confidence calculation."""
        metrics = calculate_confidence(
            variants=high_quality_variants,
            total_defining_positions=2,
            call_quality="definitive",
        )

        assert metrics.score > 0.7
        assert len(metrics.flags) == 0
        assert metrics.variant_depth_stats is not None
        assert metrics.variant_depth_stats.mean_depth > 0

    def test_low_confidence(self, low_quality_variants):
        """Test low confidence calculation."""
        metrics = calculate_confidence(
            variants=low_quality_variants,
            total_defining_positions=1,
            call_quality="definitive",
        )

        assert metrics.score < 0.5
        assert len(metrics.flags) > 0

    def test_empty_variants(self):
        """Test confidence with no variants."""
        metrics = calculate_confidence(
            variants=[],
            total_defining_positions=0,
            call_quality="definitive",
        )

        # Empty variants should have moderate score
        assert 0.0 <= metrics.score <= 1.0

    def test_custom_config(self, high_quality_variants):
        """Test confidence with custom configuration."""
        config = ConfidenceConfig(
            depth_weight=0.5,
            quality_weight=0.3,
            completeness_weight=0.2,
            target_depth=50,
        )

        metrics = calculate_confidence(
            variants=high_quality_variants,
            total_defining_positions=2,
            call_quality="definitive",
            config=config,
        )

        assert 0.0 <= metrics.score <= 1.0


class TestCalculateConfidenceForReferenceCall:
    """Tests for calculate_confidence_for_reference_call function."""

    def test_high_coverage_reference(self):
        """Test reference call confidence with high coverage."""
        coverage = {94762706: 45, 94781859: 50}

        metrics = calculate_confidence_for_reference_call(coverage)

        assert metrics.score > 0.9
        assert len(metrics.flags) == 0

    def test_low_coverage_reference(self):
        """Test reference call confidence with low coverage."""
        coverage = {94762706: 5, 94781859: 8}

        metrics = calculate_confidence_for_reference_call(coverage)

        # With low depth (mean=6.5, target=30), depth_score=0.217
        # But quality=1.0 and completeness=1.0 for reference calls
        # So composite = 0.4*0.217 + 0.35*1.0 + 0.25*1.0 = 0.69
        assert metrics.score < 0.8
        assert ConfidenceFlag.LOW_DEPTH in metrics.flags

    def test_empty_coverage(self):
        """Test reference call confidence with no coverage data."""
        metrics = calculate_confidence_for_reference_call({})

        assert metrics.score == 0.5
        assert ConfidenceFlag.INCOMPLETE_ALLELE in metrics.flags


class TestConfidenceConfig:
    """Tests for ConfidenceConfig dataclass."""

    def test_default_values(self):
        """Test default configuration values."""
        config = ConfidenceConfig()

        assert config.depth_weight == 0.40
        assert config.quality_weight == 0.35
        assert config.completeness_weight == 0.25
        assert config.target_depth == 30
        assert config.min_depth_threshold == 10
        assert config.min_quality_threshold == 20.0

    def test_weights_sum_to_one(self):
        """Test that default weights sum to 1.0."""
        config = ConfidenceConfig()

        total = config.depth_weight + config.quality_weight + config.completeness_weight
        assert abs(total - 1.0) < 0.001

    def test_custom_values(self):
        """Test custom configuration values."""
        config = ConfidenceConfig(
            depth_weight=0.5,
            quality_weight=0.3,
            completeness_weight=0.2,
            target_depth=50,
        )

        assert config.depth_weight == 0.5
        assert config.target_depth == 50
