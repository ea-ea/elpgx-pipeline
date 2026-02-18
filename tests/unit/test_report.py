"""Unit tests for pgx.report module."""

import json
from datetime import datetime
from pathlib import Path

import pytest

from elpgx.allele_caller import DiplotypeCall
from elpgx.models import BiomarkerType, GeneTier, Variant
from elpgx.phenotype import PhenotypeResult
from elpgx.report import (
    build_gene_result,
    build_report,
    generate_report,
    load_gene_tier_info,
    write_report,
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
def gene_tiers_dir(tmp_path):
    """Create temporary gene tiers file."""
    gene_tiers = {
        "genes": {
            "CYP2C19": {
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

    return tmp_path


@pytest.fixture
def data_dir(tmp_path):
    """Create temporary data directory with all required files."""
    # Create gene_tiers.json
    gene_tiers = {
        "genes": {
            "CYP2C19": {
                "tier": "group1",
                "biomarker_type": "metabolic_enzyme",
            }
        }
    }
    tiers_file = tmp_path / "gene_tiers.json"
    with open(tiers_file, "w") as f:
        json.dump(gene_tiers, f)

    return tmp_path


class TestLoadGeneTierInfo:
    """Tests for load_gene_tier_info function."""

    def test_load_existing_file(self, gene_tiers_dir):
        """Test loading existing gene tiers file."""
        tiers_file = gene_tiers_dir / "gene_tiers.json"
        info = load_gene_tier_info(tiers_file)

        assert "CYP2C19" in info
        assert info["CYP2C19"]["tier"] == "group1"
        assert info["CYP2C19"]["biomarker_type"] == "metabolic_enzyme"

    def test_load_nonexistent_file(self, tmp_path):
        """Test loading non-existent file returns empty dict."""
        fake_path = tmp_path / "nonexistent.json"
        info = load_gene_tier_info(fake_path)

        assert info == {}


class TestBuildGeneResult:
    """Tests for build_gene_result function."""

    def test_build_gene_result(
        self, sample_diplotype_call, sample_phenotype_result, gene_tiers_dir
    ):
        """Test building a gene result."""
        tier_info = {"tier": "group1", "biomarker_type": "metabolic_enzyme"}

        result = build_gene_result(
            diplotype_call=sample_diplotype_call,
            phenotype_result=sample_phenotype_result,
            gene_tier_info=tier_info,
        )

        assert result.gene == "CYP2C19"
        assert result.diplotype == "*1/*17"
        assert result.phenotype == "Rapid Metabolizer"
        assert result.gene_tier == GeneTier.GROUP1
        assert result.biomarker_type == BiomarkerType.METABOLIC_ENZYME
        assert result.confidence.score > 0

    def test_build_gene_result_default_tier(
        self, sample_diplotype_call, sample_phenotype_result
    ):
        """Test building gene result with default tier info."""
        result = build_gene_result(
            diplotype_call=sample_diplotype_call,
            phenotype_result=sample_phenotype_result,
            gene_tier_info={},  # Empty tier info
        )

        assert result.gene_tier == GeneTier.GROUP1  # Default
        assert result.biomarker_type == BiomarkerType.METABOLIC_ENZYME  # Default


class TestBuildReport:
    """Tests for build_report function."""

    def test_build_report(
        self, sample_diplotype_call, sample_phenotype_result, gene_tiers_dir
    ):
        """Test building complete report."""
        tiers_file = gene_tiers_dir / "gene_tiers.json"

        report = build_report(
            sample_id="TEST_001",
            diplotype_calls=[sample_diplotype_call],
            phenotype_results=[sample_phenotype_result],
            gene_tiers_file=tiers_file,
            gene_tiers_included=["group1"],
            hla_enabled=False,
            pgx_db_version="CPIC v1.22",
            pipeline_version="1.0.0",
        )

        assert report.meta.sample_id == "TEST_001"
        assert report.meta.pipeline_version == "1.0.0"
        assert len(report.results) == 1
        assert report.results[0].gene == "CYP2C19"

    def test_build_report_empty_results(self, gene_tiers_dir):
        """Test building report with no results."""
        tiers_file = gene_tiers_dir / "gene_tiers.json"

        report = build_report(
            sample_id="TEST_002",
            diplotype_calls=[],
            phenotype_results=[],
            gene_tiers_file=tiers_file,
            gene_tiers_included=["group1"],
        )

        assert report.meta.sample_id == "TEST_002"
        assert len(report.results) == 0

    def test_build_report_with_hla(
        self, sample_diplotype_call, sample_phenotype_result, gene_tiers_dir
    ):
        """Test building report with HLA module enabled."""
        tiers_file = gene_tiers_dir / "gene_tiers.json"

        report = build_report(
            sample_id="TEST_003",
            diplotype_calls=[sample_diplotype_call],
            phenotype_results=[sample_phenotype_result],
            gene_tiers_file=tiers_file,
            gene_tiers_included=["group1"],
            hla_enabled=True,
        )

        assert report.meta.hla_module_enabled is True


class TestWriteReport:
    """Tests for write_report function."""

    def test_write_report(
        self, sample_diplotype_call, sample_phenotype_result, gene_tiers_dir, tmp_path
    ):
        """Test writing report to file."""
        tiers_file = gene_tiers_dir / "gene_tiers.json"

        report = build_report(
            sample_id="TEST_001",
            diplotype_calls=[sample_diplotype_call],
            phenotype_results=[sample_phenotype_result],
            gene_tiers_file=tiers_file,
            gene_tiers_included=["group1"],
        )

        output_path = tmp_path / "output" / "report.json"
        write_report(report, output_path)

        assert output_path.exists()

        with open(output_path) as f:
            data = json.load(f)

        assert data["meta"]["sample_id"] == "TEST_001"
        assert len(data["results"]) == 1


class TestGenerateReport:
    """Tests for generate_report function."""

    def test_generate_report_integration(
        self, sample_diplotype_call, sample_phenotype_result, data_dir, tmp_path
    ):
        """Test complete report generation integration."""
        output_path = tmp_path / "test_report.json"

        report = generate_report(
            sample_id="TEST_SAMPLE",
            diplotype_calls=[sample_diplotype_call],
            phenotype_results=[sample_phenotype_result],
            output_path=output_path,
            data_dir=data_dir,
            gene_tiers_included=["group1"],
            pgx_db_version="CPIC v1.22",
            pipeline_version="1.0.0",
        )

        assert output_path.exists()
        assert report.meta.sample_id == "TEST_SAMPLE"

        # Verify file contents
        with open(output_path) as f:
            data = json.load(f)

        assert "meta" in data
        assert "results" in data


class TestReportSerialization:
    """Tests for report JSON serialization."""

    def test_report_json_structure(
        self, sample_diplotype_call, sample_phenotype_result, gene_tiers_dir
    ):
        """Test report JSON structure matches expected schema."""
        tiers_file = gene_tiers_dir / "gene_tiers.json"

        report = build_report(
            sample_id="TEST_001",
            diplotype_calls=[sample_diplotype_call],
            phenotype_results=[sample_phenotype_result],
            gene_tiers_file=tiers_file,
            gene_tiers_included=["group1"],
        )

        json_str = report.to_json()
        data = json.loads(json_str)

        # Check required fields
        assert "meta" in data
        assert "results" in data
        assert "warnings" in data

        # Check meta fields
        assert "sample_id" in data["meta"]
        assert "pipeline_version" in data["meta"]
        assert "run_date" in data["meta"]
        assert "reference_build" in data["meta"]

        # Check results structure
        if data["results"]:
            result = data["results"][0]
            assert "gene" in result
            assert "diplotype" in result
            assert "phenotype" in result
            assert "confidence" in result

    def test_datetime_serialization(
        self, sample_diplotype_call, sample_phenotype_result, gene_tiers_dir
    ):
        """Test datetime is properly serialized."""
        tiers_file = gene_tiers_dir / "gene_tiers.json"

        report = build_report(
            sample_id="TEST_001",
            diplotype_calls=[sample_diplotype_call],
            phenotype_results=[sample_phenotype_result],
            gene_tiers_file=tiers_file,
            gene_tiers_included=["group1"],
        )

        json_str = report.to_json()
        data = json.loads(json_str)

        # run_date should be a valid ISO format string
        run_date = data["meta"]["run_date"]
        assert isinstance(run_date, str)
        # Should be parseable as datetime
        datetime.fromisoformat(run_date.replace("Z", "+00:00"))
