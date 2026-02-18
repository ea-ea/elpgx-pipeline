"""Unit tests for elpgx.models module."""

import json
from datetime import datetime

import pytest

from elpgx.models import (
    AlleleDefinition,
    AlleleVariant,
    BiomarkerType,
    ConfidenceFlag,
    ConfidenceMetrics,
    GeneAlleleDefinitions,
    GeneResult,
    GeneTier,
    PgxReport,
    PipelineMeta,
    Variant,
    VariantDepthStats,
)


class TestVariant:
    """Tests for Variant model."""

    def test_create_minimal_variant(self):
        """Test creating a variant with required fields only."""
        v = Variant(
            chrom="chr10",
            pos=94762706,
            ref="C",
            alt="T",
            genotype="0/1",
            quality=38.5,
            depth=42,
        )
        assert v.chrom == "chr10"
        assert v.pos == 94762706
        assert v.genotype == "0/1"
        assert v.rsid is None
        assert v.allele_defining is False

    def test_create_full_variant(self):
        """Test creating a variant with all fields."""
        v = Variant(
            rsid="rs12248560",
            chrom="chr10",
            pos=94762706,
            ref="C",
            alt="T",
            genotype="0/1",
            quality=38.5,
            depth=42,
            allele_defining=True,
            star_allele="*17",
            functional_impact="increased function",
        )
        assert v.rsid == "rs12248560"
        assert v.allele_defining is True
        assert v.star_allele == "*17"

    def test_variant_serialization(self):
        """Test variant JSON serialization."""
        v = Variant(
            chrom="chr10",
            pos=94762706,
            ref="C",
            alt="T",
            genotype="0/1",
            quality=38.5,
            depth=42,
        )
        json_str = v.model_dump_json()
        data = json.loads(json_str)
        assert data["chrom"] == "chr10"
        assert data["pos"] == 94762706


class TestConfidenceMetrics:
    """Tests for ConfidenceMetrics model."""

    def test_create_confidence_metrics(self):
        """Test creating confidence metrics."""
        cm = ConfidenceMetrics(
            score=0.95,
            flags=[],
            variant_depth_stats=VariantDepthStats(mean_depth=42.0, min_depth=28),
        )
        assert cm.score == 0.95
        assert cm.flags == []
        assert cm.variant_depth_stats.mean_depth == 42.0

    def test_confidence_score_bounds(self):
        """Test confidence score must be between 0 and 1."""
        with pytest.raises(ValueError):
            ConfidenceMetrics(score=1.5, flags=[])

        with pytest.raises(ValueError):
            ConfidenceMetrics(score=-0.1, flags=[])

    def test_confidence_with_flags(self):
        """Test confidence metrics with quality flags."""
        cm = ConfidenceMetrics(
            score=0.7,
            flags=[ConfidenceFlag.LOW_DEPTH, ConfidenceFlag.AMBIGUOUS_CALL],
        )
        assert len(cm.flags) == 2
        assert ConfidenceFlag.LOW_DEPTH in cm.flags


class TestGeneResult:
    """Tests for GeneResult model."""

    def test_create_gene_result(self):
        """Test creating a full gene result."""
        confidence = ConfidenceMetrics(score=0.95, flags=[])
        result = GeneResult(
            gene="CYP2C19",
            gene_tier=GeneTier.GROUP1,
            biomarker_type=BiomarkerType.METABOLIC_ENZYME,
            diplotype="*1/*17",
            functional_status="Normal Function / Increased Function",
            phenotype="Rapid Metabolizer",
            activity_score=2.5,
            prescribing_implication="Consider alternative.",
            confidence=confidence,
            variants=[],
        )
        assert result.gene == "CYP2C19"
        assert result.diplotype == "*1/*17"
        assert result.phenotype == "Rapid Metabolizer"

    def test_gene_tier_enum(self):
        """Test GeneTier enum values."""
        assert GeneTier.GROUP1.value == "group1"
        assert GeneTier.GROUP2.value == "group2"

    def test_biomarker_type_enum(self):
        """Test BiomarkerType enum values."""
        assert BiomarkerType.METABOLIC_ENZYME.value == "metabolic_enzyme"
        assert BiomarkerType.TRANSPORTER.value == "transporter"
        assert BiomarkerType.HLA.value == "hla"


class TestPgxReport:
    """Tests for PgxReport model."""

    def test_create_empty_report(self):
        """Test creating report with no results."""
        meta = PipelineMeta(
            sample_id="TEST_001",
            pipeline_version="1.0.0",
            pgx_db_version="CPIC v1.22",
        )
        report = PgxReport(meta=meta, results=[], warnings=[])
        assert report.meta.sample_id == "TEST_001"
        assert len(report.results) == 0

    def test_report_to_json(self):
        """Test report JSON serialization."""
        meta = PipelineMeta(
            sample_id="TEST_001",
            pipeline_version="1.0.0",
            pgx_db_version="CPIC v1.22",
            gene_tiers_included=[GeneTier.GROUP1],
        )
        report = PgxReport(meta=meta, results=[], warnings=[])

        json_str = report.to_json()
        data = json.loads(json_str)

        assert data["meta"]["sample_id"] == "TEST_001"
        assert "group1" in data["meta"]["gene_tiers_included"]

    def test_report_with_results(self):
        """Test report with gene results."""
        meta = PipelineMeta(
            sample_id="TEST_001",
            pipeline_version="1.0.0",
            pgx_db_version="CPIC v1.22",
        )

        confidence = ConfidenceMetrics(score=0.95, flags=[])
        result = GeneResult(
            gene="CYP2C19",
            gene_tier=GeneTier.GROUP1,
            biomarker_type=BiomarkerType.METABOLIC_ENZYME,
            diplotype="*1/*1",
            functional_status="Normal Function",
            phenotype="Normal Metabolizer",
            prescribing_implication="Standard dosing.",
            confidence=confidence,
        )

        report = PgxReport(meta=meta, results=[result], warnings=[])
        assert len(report.results) == 1
        assert report.results[0].gene == "CYP2C19"


class TestAlleleDefinition:
    """Tests for allele definition models."""

    def test_create_allele_definition(self):
        """Test creating an allele definition."""
        variant = AlleleVariant(
            rsid="rs12248560",
            chrom="chr10",
            pos=94762706,
            ref="C",
            alt="T",
            functional_impact="increased function",
        )
        allele = AlleleDefinition(
            allele="*17",
            function="Increased Function",
            activity_value=1.5,
            defining_variants=[variant],
        )
        assert allele.allele == "*17"
        assert len(allele.defining_variants) == 1

    def test_create_gene_allele_definitions(self):
        """Test creating gene-level allele definitions."""
        gene_def = GeneAlleleDefinitions(
            gene="CYP2C19",
            chromosome="chr10",
            gene_tier=GeneTier.GROUP1,
            biomarker_type=BiomarkerType.METABOLIC_ENZYME,
            reference_allele="*1",
            alleles=[],
        )
        assert gene_def.gene == "CYP2C19"
        assert gene_def.reference_allele == "*1"
