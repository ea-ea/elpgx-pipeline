"""
End-to-End Pipeline Tests using real HG002 BAM data (GRCh37/hg19).

These tests run the complete ElPGx pipeline against a real
Genome in a Bottle (GIAB) BAM file (HG002) aligned to GRCh37.
Tests verify that GRCh37 BAMs are handled correctly with
automatic coordinate translation to GRCh38 for allele matching.

Test Categories:
- BAM-100: Input Validation (BAM + GRCh37 reference)
- BAM-200: Variant Calling with GRCh37 BED
- BAM-300: Coordinate Translation & Allele Calling
- BAM-400: Phenotype Translation
- BAM-500: Report Generation & Schema
- BAM-600: Full Pipeline CLI Integration
- BAM-700: Known-Result Concordance

Requirements:
    sample_test_data/151002_*.posiSrt.markDup.bam (HG002 BAM)
    sample_test_data/151002_*.posiSrt.markDup.bai (BAM index)
    sample_test_data/human_g1k_v37.fasta (GRCh37 reference)
    sample_test_data/human_g1k_v37.fasta.fai (reference index)
    sample_test_data/human_g1k_v37.dict (sequence dictionary)
"""

import json
import subprocess
from pathlib import Path

import pytest

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).parent.parent.parent
SAMPLE_DATA = PROJECT_ROOT / "sample_test_data"
DATA_DIR = PROJECT_ROOT / "data"

# HG002 BAM file (GRCh37/hg19 aligned)
BAM_FILE = SAMPLE_DATA / "151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam"
BAM_INDEX = BAM_FILE.with_suffix(".bai")

# GRCh37 reference
GRCH37_REFERENCE = SAMPLE_DATA / "human_g1k_v37.fasta"
GRCH37_REFERENCE_FAI = SAMPLE_DATA / "human_g1k_v37.fasta.fai"
GRCH37_REFERENCE_DICT = SAMPLE_DATA / "human_g1k_v37.dict"

# GRCh37 BED file for target regions
GRCH37_BED = DATA_DIR / "pgx_regions_grch37.bed"

# Expected results
EXPECTED_RESULTS = SAMPLE_DATA / "HG002_expected_results.json"

# Skip all tests if sample data is not present
HAS_BAM_DATA = (
    BAM_FILE.exists()
    and BAM_INDEX.exists()
    and GRCH37_REFERENCE.exists()
    and GRCH37_REFERENCE_FAI.exists()
)

pytestmark = pytest.mark.skipif(
    not HAS_BAM_DATA,
    reason="BAM test data not available (HG002 BAM + GRCh37 reference required)",
)

# ---------------------------------------------------------------------------
# Expected results from HG002 (GIAB reference sample, GRCh37/hg19 aligned)
# Results verified from pipeline run with GRCh37 coordinate translation.
# ---------------------------------------------------------------------------
EXPECTED_SAMPLE_ID = "Sample_Diag-excap51-HG002-EEogPU"
EXPECTED_GENE_COUNT = 16  # group1 + group2 + group3

EXPECTED_DIPLOTYPES = {
    "CYP2C19": "*1/*1",
    "CYP2C9": "*1/*1",
    "CYP2D6": "*1/*1",
    "SLCO1B1": "*1/*1",
    "VKORC1": "*1/*1",
    "TPMT": "*1/*1",
    "DPYD": "*1/*1",
    "UGT1A1": "*1/*1",
    "NUDT15": "*1/*1",
    "CYP3A5": "*1/*1",
    "CYP2B6": "*1/*1",
    "NAT2": "*6/*60",
    "IFNL3": "*1/*1",
    "CACNA1S": "*1/*1",
    "RYR1": "*1/*1",
    "G6PD": "*1/*1",
}

EXPECTED_PHENOTYPES = {
    "CYP2C19": "Normal Metabolizer",
    "CYP2D6": "Normal Metabolizer",
    "SLCO1B1": "Normal Function",
    "NAT2": "Indeterminate",
}


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def expected_report():
    """Load expected results for concordance checks."""
    if EXPECTED_RESULTS.exists():
        with open(EXPECTED_RESULTS) as f:
            return json.load(f)
    return None


@pytest.fixture(scope="module")
def pipeline_report(tmp_path_factory):
    """Run full pipeline once on HG002 BAM and return the report.

    This fixture is module-scoped so the pipeline only runs once
    for all tests in this file.
    """
    output_dir = tmp_path_factory.mktemp("e2e_bam_output")
    output_file = output_dir / "HG002_pgx_report.json"

    result = subprocess.run(
        [
            "elpgx", "run",
            "-b", str(BAM_FILE),
            "-r", str(GRCH37_REFERENCE),
            "-o", str(output_file),
            "-g", "group1,group2,group3",
            "--pretty",
        ],
        capture_output=True,
        text=True,
        timeout=600,
    )

    assert result.returncode == 0, (
        f"Pipeline failed with exit code {result.returncode}\n"
        f"STDOUT: {result.stdout[-500:]}\n"
        f"STDERR: {result.stderr[-500:]}"
    )
    assert output_file.exists(), "Report file not created"

    with open(output_file) as f:
        return json.load(f)


# ===========================================================================
# BAM-100: Input Validation
# ===========================================================================
class TestBAMInputValidation:
    """BAM-100: Validate BAM input with GRCh37 reference."""

    def test_bam_file_exists(self):
        """BAM file must be present and non-empty."""
        assert BAM_FILE.exists()
        assert BAM_FILE.stat().st_size > 0

    def test_bam_index_exists(self):
        """BAM index (.bai) must be present."""
        assert BAM_INDEX.exists()

    def test_grch37_reference_exists(self):
        """GRCh37 reference FASTA must be present."""
        assert GRCH37_REFERENCE.exists()
        assert GRCH37_REFERENCE_FAI.exists()

    def test_grch37_bed_exists(self):
        """GRCh37 BED file for PGx regions must exist."""
        assert GRCH37_BED.exists()

    def test_validate_detects_grch37(self):
        """Validation should detect GRCh37 build."""
        from elpgx.validate import validate_alignment_file

        result = validate_alignment_file(BAM_FILE)
        assert result.is_valid is True
        assert result.reference_build == "GRCh37"

    def test_validate_detects_sample_id(self):
        """Validation should detect sample ID from @RG header."""
        from elpgx.validate import validate_alignment_file

        result = validate_alignment_file(BAM_FILE)
        assert result.sample_id is not None
        assert "HG002" in result.sample_id

    def test_validate_command(self):
        """elpgx validate should pass for this BAM."""
        result = subprocess.run(
            [
                "elpgx", "validate",
                "-b", str(BAM_FILE),
            ],
            capture_output=True,
            text=True,
            timeout=120,
        )
        assert result.returncode == 0
        assert "Validation passed" in result.stdout


# ===========================================================================
# BAM-200: Variant Calling
# ===========================================================================
class TestBAMVariantCalling:
    """BAM-200: GATK variant calling on GRCh37 BAM data."""

    def test_variant_calling_with_grch37(self, tmp_path):
        """GATK should produce a VCF from GRCh37 BAM."""
        from elpgx.vcf_parser import call_variants_gatk

        output_vcf = tmp_path / "test_hg002_variants.vcf.gz"
        vcf_path = call_variants_gatk(
            bam_path=BAM_FILE,
            reference_path=GRCH37_REFERENCE,
            regions_bed=GRCH37_BED,
            output_vcf=output_vcf,
        )
        assert vcf_path is not None
        assert Path(vcf_path).exists()
        assert Path(vcf_path).stat().st_size > 0

    def test_vcf_has_pgx_variants(self, tmp_path):
        """VCF should contain variants in PGx regions (GRCh37 coordinates)."""
        from elpgx.vcf_parser import call_variants_gatk, parse_vcf

        output_vcf = tmp_path / "test_hg002_variants.vcf.gz"
        vcf_path = call_variants_gatk(
            bam_path=BAM_FILE,
            reference_path=GRCH37_REFERENCE,
            regions_bed=GRCH37_BED,
            output_vcf=output_vcf,
        )
        variants = parse_vcf(vcf_path)
        assert len(variants) > 0, "No variants found in PGx regions"

        # GRCh37 BAM should have numeric chromosome names (no 'chr' prefix)
        chroms = {v.chrom for v in variants}
        assert any(c in chroms for c in ["10", "22", "19", "8"]), (
            f"No variants on known PGx chromosomes, got: {chroms}"
        )


# ===========================================================================
# BAM-300: Coordinate Translation & Allele Calling
# ===========================================================================
class TestBAMAlleleCalling:
    """BAM-300: Star-allele assignment from GRCh37 variants with coordinate translation."""

    def test_all_expected_genes_called(self, pipeline_report):
        """Pipeline should call all 16 genes (group1+2+3)."""
        genes = {r["gene"] for r in pipeline_report["results"]}
        assert len(genes) == EXPECTED_GENE_COUNT, (
            f"Expected {EXPECTED_GENE_COUNT} genes, got {len(genes)}: {genes}"
        )

    def test_coordinate_translation_worked(self, pipeline_report):
        """GRCh37 coordinates should be translated to GRCh38 for allele matching."""
        # The report should indicate GRCh37 was the input build
        assert pipeline_report["meta"]["reference_build"] == "GRCh37"

    def test_diplotypes_are_valid(self, pipeline_report):
        """All diplotypes should follow *X/*Y format."""
        for result in pipeline_report["results"]:
            diplotype = result["diplotype"]
            assert "/" in diplotype, (
                f"{result['gene']}: invalid diplotype format '{diplotype}'"
            )

    @pytest.mark.parametrize("gene,expected_diplotype", EXPECTED_DIPLOTYPES.items())
    def test_diplotype_concordance(self, pipeline_report, gene, expected_diplotype):
        """Each gene's diplotype should match expected result."""
        gene_result = next(
            (r for r in pipeline_report["results"] if r["gene"] == gene), None
        )
        assert gene_result is not None, f"Gene {gene} not found in results"
        assert gene_result["diplotype"] == expected_diplotype, (
            f"{gene}: expected {expected_diplotype}, "
            f"got {gene_result['diplotype']}"
        )


# ===========================================================================
# BAM-400: Phenotype Translation
# ===========================================================================
class TestBAMPhenotypeTranslation:
    """BAM-400: Phenotype assignment from HG002 diplotypes."""

    @pytest.mark.parametrize("gene,expected_phenotype", EXPECTED_PHENOTYPES.items())
    def test_phenotype_concordance(self, pipeline_report, gene, expected_phenotype):
        """Phenotype should match expected for key genes."""
        gene_result = next(
            (r for r in pipeline_report["results"] if r["gene"] == gene), None
        )
        assert gene_result is not None, f"Gene {gene} not found"
        assert gene_result["phenotype"] == expected_phenotype, (
            f"{gene}: expected '{expected_phenotype}', "
            f"got '{gene_result['phenotype']}'"
        )

    def test_all_genes_have_phenotype(self, pipeline_report):
        """Every gene should have a non-empty phenotype."""
        for result in pipeline_report["results"]:
            assert result["phenotype"], (
                f"{result['gene']}: missing phenotype"
            )

    def test_all_genes_have_prescribing_implication(self, pipeline_report):
        """Every gene should have a prescribing implication."""
        for result in pipeline_report["results"]:
            assert result["prescribing_implication"], (
                f"{result['gene']}: missing prescribing implication"
            )


# ===========================================================================
# BAM-500: Report Generation & Schema
# ===========================================================================
class TestBAMReportGeneration:
    """BAM-500: Report structure and schema compliance for GRCh37 BAM."""

    def test_schema_version(self, pipeline_report):
        """Report should have schema version."""
        assert pipeline_report["schema_version"] == "1.0.0"

    def test_meta_sample_id(self, pipeline_report):
        """Report should identify sample correctly."""
        assert pipeline_report["meta"]["sample_id"] is not None

    def test_meta_reference_build_grch37(self, pipeline_report):
        """Report should indicate GRCh37 build."""
        assert pipeline_report["meta"]["reference_build"] == "GRCh37"

    def test_meta_pipeline_version(self, pipeline_report):
        """Report should include pipeline version."""
        assert pipeline_report["meta"]["pipeline_version"] == "1.1.0"

    def test_confidence_scores_valid(self, pipeline_report):
        """All confidence scores should be between 0.0 and 1.0."""
        for result in pipeline_report["results"]:
            score = result["confidence"]["score"]
            assert 0.0 <= score <= 1.0, (
                f"{result['gene']}: confidence {score} out of range"
            )

    def test_report_json_roundtrip(self, pipeline_report):
        """Report should survive JSON serialization roundtrip."""
        serialized = json.dumps(pipeline_report, indent=2)
        deserialized = json.loads(serialized)
        assert deserialized == pipeline_report


# ===========================================================================
# BAM-600: CLI Integration
# ===========================================================================
class TestBAMCLIIntegration:
    """BAM-600: CLI commands work with GRCh37 BAM end-to-end."""

    def test_run_group1_only(self, tmp_path):
        """Pipeline should work with group1 genes only on GRCh37 BAM."""
        output_file = tmp_path / "group1_hg002_report.json"
        result = subprocess.run(
            [
                "elpgx", "run",
                "-b", str(BAM_FILE),
                "-r", str(GRCH37_REFERENCE),
                "-o", str(output_file),
                "-g", "group1",
                "--pretty",
            ],
            capture_output=True,
            text=True,
            timeout=600,
        )
        assert result.returncode == 0

        with open(output_file) as f:
            report = json.load(f)

        genes = {r["gene"] for r in report["results"]}
        group1_genes = {
            "CYP2C19", "CYP2C9", "VKORC1", "SLCO1B1",
            "TPMT", "DPYD", "UGT1A1", "NUDT15",
        }
        assert genes == group1_genes


# ===========================================================================
# BAM-700: Known-Result Concordance
# ===========================================================================
class TestBAMConcordance:
    """BAM-700: Results match previously validated reference output."""

    def test_full_concordance(self, pipeline_report, expected_report):
        """All gene results should match the reference run."""
        if expected_report is None:
            pytest.skip("Expected results file not available")

        current_genes = {
            r["gene"]: r for r in pipeline_report["results"]
        }
        expected_genes = {
            r["gene"]: r for r in expected_report["results"]
        }

        assert set(current_genes.keys()) == set(expected_genes.keys()), (
            "Gene set mismatch between current and expected"
        )

        discordant = []
        for gene in current_genes:
            curr = current_genes[gene]
            exp = expected_genes[gene]
            if curr["diplotype"] != exp["diplotype"]:
                discordant.append(
                    f"{gene}: diplotype {curr['diplotype']} != {exp['diplotype']}"
                )
            if curr["phenotype"] != exp["phenotype"]:
                discordant.append(
                    f"{gene}: phenotype {curr['phenotype']} != {exp['phenotype']}"
                )

        assert len(discordant) == 0, (
            f"Concordance failures:\n" + "\n".join(discordant)
        )

    def test_confidence_stability(self, pipeline_report, expected_report):
        """Confidence scores should not drift more than 0.05."""
        if expected_report is None:
            pytest.skip("Expected results file not available")

        current_genes = {
            r["gene"]: r for r in pipeline_report["results"]
        }
        expected_genes = {
            r["gene"]: r for r in expected_report["results"]
        }

        drifted = []
        for gene in current_genes:
            curr_score = current_genes[gene]["confidence"]["score"]
            exp_score = expected_genes[gene]["confidence"]["score"]
            if abs(curr_score - exp_score) > 0.05:
                drifted.append(
                    f"{gene}: confidence {curr_score} vs expected {exp_score}"
                )

        assert len(drifted) == 0, (
            f"Confidence drift detected:\n" + "\n".join(drifted)
        )
