"""
End-to-End Pipeline Tests using real HG00100 CRAM data.

These tests run the complete ElPGx pipeline against a real
1000 Genomes exome CRAM file (HG00100) with GRCh38 reference.

Test Categories:
- E2E-100: Input Validation (CRAM + reference)
- E2E-200: Variant Calling (GATK HaplotypeCaller)
- E2E-300: Allele Calling (star-allele assignment)
- E2E-400: Phenotype Translation
- E2E-500: Report Generation & Schema Compliance
- E2E-600: Full Pipeline CLI Integration
- E2E-700: Known-Result Concordance

Requirements:
    sample_test_data/HG00100.alt_bwamem_GRCh38DH.20150826.GBR.exome.cram
    sample_test_data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta
    sample_test_data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.fai
    sample_test_data/GCA_000001405.15_GRCh38_no_alt_analysis_set.dict
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
CRAM_FILE = SAMPLE_DATA / "HG00100.alt_bwamem_GRCh38DH.20150826.GBR.exome.cram"
CRAM_INDEX = CRAM_FILE.with_suffix(".cram.crai")
REFERENCE = SAMPLE_DATA / "GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"
REFERENCE_FAI = REFERENCE.with_suffix(".fasta.fai")
REFERENCE_DICT = SAMPLE_DATA / "GCA_000001405.15_GRCh38_no_alt_analysis_set.dict"
EXPECTED_RESULTS = SAMPLE_DATA / "HG00100_expected_results.json"
DATA_DIR = PROJECT_ROOT / "data"
BED_FILE = DATA_DIR / "pgx_regions.bed"

# Skip all tests if sample data is not present
HAS_SAMPLE_DATA = (
    CRAM_FILE.exists()
    and CRAM_INDEX.exists()
    and REFERENCE.exists()
    and REFERENCE_FAI.exists()
)

pytestmark = pytest.mark.skipif(
    not HAS_SAMPLE_DATA,
    reason="Sample test data not available (CRAM/FASTA files required)",
)

# ---------------------------------------------------------------------------
# Expected results from HG00100 (1000 Genomes, GBR population)
# ---------------------------------------------------------------------------
# These are verified reference calls from the pipeline run.
EXPECTED_SAMPLE_ID = "HG00100"
EXPECTED_GENE_COUNT = 16  # group1 + group2 + group3

EXPECTED_DIPLOTYPES = {
    "CYP2C19": "*1/*1",
    "CYP2C9": "*1/*1",
    "CYP2D6": "*13+*2/*13+*2",
    "SLCO1B1": "*1/*1",
    "VKORC1": "*1/*1",
    "TPMT": "*1/*1",
    "DPYD": "*1/*1",
    "UGT1A1": "*1/*1",
    "NUDT15": "*1/*1",
    "CYP3A5": "*1/*1",
    "CYP2B6": "*1/*6",
    "NAT2": "*6/*6",
    "IFNL3": "*1/*1",
    "CACNA1S": "*1/*1",
    "RYR1": "*1/*1",
    "G6PD": "*1/*1",
}

EXPECTED_PHENOTYPES = {
    "CYP2C19": "Normal Metabolizer",
    "CYP2B6": "Intermediate Metabolizer",
    "NAT2": "Poor Metabolizer",
    "SLCO1B1": "Normal Function",
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
    """Run full pipeline once and return the report.

    This fixture is module-scoped so the pipeline only runs once
    for all tests in this file.
    """
    output_dir = tmp_path_factory.mktemp("e2e_output")
    output_file = output_dir / "HG00100_pgx_report.json"

    result = subprocess.run(
        [
            "elpgx", "run",
            "-b", str(CRAM_FILE),
            "-r", str(REFERENCE),
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
# E2E-100: Input Validation
# ===========================================================================
class TestInputValidation:
    """E2E-100: Validate CRAM input with real data."""

    def test_cram_file_exists(self):
        """CRAM file must be present and non-empty."""
        assert CRAM_FILE.exists()
        assert CRAM_FILE.stat().st_size > 0

    def test_cram_index_exists(self):
        """CRAM index (.crai) must be present."""
        assert CRAM_INDEX.exists()

    def test_reference_exists(self):
        """Reference FASTA must be present."""
        assert REFERENCE.exists()
        assert REFERENCE_FAI.exists()

    def test_reference_dict_exists(self):
        """Reference .dict must be present for GATK."""
        assert REFERENCE_DICT.exists()

    def test_validate_command(self):
        """elpgx validate should pass for this CRAM."""
        result = subprocess.run(
            [
                "elpgx", "validate",
                "-b", str(CRAM_FILE),
                "-r", str(REFERENCE),
            ],
            capture_output=True,
            text=True,
            timeout=120,
        )
        assert result.returncode == 0
        assert "Validation passed" in result.stdout

    def test_validate_detects_sample_id(self):
        """Validation should detect sample ID HG00100."""
        result = subprocess.run(
            [
                "elpgx", "validate",
                "-b", str(CRAM_FILE),
                "-r", str(REFERENCE),
            ],
            capture_output=True,
            text=True,
            timeout=120,
        )
        assert "HG00100" in result.stdout

    def test_validate_api(self):
        """Direct API validation should work."""
        from elpgx.validate import validate_alignment_file

        result = validate_alignment_file(
            CRAM_FILE,
            reference_path=REFERENCE,
        )
        assert result.is_valid is True
        assert result.sample_id == "HG00100"


# ===========================================================================
# E2E-200: Variant Calling
# ===========================================================================
class TestVariantCalling:
    """E2E-200: GATK variant calling on real data."""

    def test_variant_calling_produces_vcf(self, tmp_path):
        """GATK should produce a VCF from the CRAM."""
        from elpgx.vcf_parser import call_variants_gatk

        output_vcf = tmp_path / "test_variants.vcf.gz"
        vcf_path = call_variants_gatk(
            bam_path=CRAM_FILE,
            reference_path=REFERENCE,
            regions_bed=BED_FILE,
            output_vcf=output_vcf,
        )
        assert vcf_path is not None
        assert Path(vcf_path).exists()
        assert Path(vcf_path).stat().st_size > 0

    def test_vcf_has_pgx_variants(self, tmp_path):
        """VCF should contain variants in PGx regions."""
        from elpgx.vcf_parser import call_variants_gatk, parse_vcf

        output_vcf = tmp_path / "test_variants.vcf.gz"
        vcf_path = call_variants_gatk(
            bam_path=CRAM_FILE,
            reference_path=REFERENCE,
            regions_bed=BED_FILE,
            output_vcf=output_vcf,
        )
        variants = parse_vcf(vcf_path)
        assert len(variants) > 0, "No variants found in PGx regions"

        # Should have variants on pharmacogene chromosomes
        chroms = {v.chrom for v in variants}
        assert any(c in chroms for c in ["chr10", "chr22", "chr19"]), (
            f"No variants on known PGx chromosomes, got: {chroms}"
        )


# ===========================================================================
# E2E-300: Allele Calling
# ===========================================================================
class TestAlleleCalling:
    """E2E-300: Star-allele assignment from real variants."""

    def test_all_expected_genes_called(self, pipeline_report):
        """Pipeline should call all 16 genes (group1+2+3)."""
        genes = {r["gene"] for r in pipeline_report["results"]}
        assert len(genes) == EXPECTED_GENE_COUNT, (
            f"Expected {EXPECTED_GENE_COUNT} genes, got {len(genes)}: {genes}"
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

    def test_cyp2b6_has_variant(self, pipeline_report):
        """CYP2B6 *6 should have at least one defining variant."""
        cyp2b6 = next(
            r for r in pipeline_report["results"] if r["gene"] == "CYP2B6"
        )
        assert len(cyp2b6["variants"]) > 0, "CYP2B6 *6 should have variants"

    def test_nat2_has_variants(self, pipeline_report):
        """NAT2 *6/*6 should have defining variants."""
        nat2 = next(
            r for r in pipeline_report["results"] if r["gene"] == "NAT2"
        )
        assert len(nat2["variants"]) > 0, "NAT2 *6/*6 should have variants"


# ===========================================================================
# E2E-400: Phenotype Translation
# ===========================================================================
class TestPhenotypeTranslation:
    """E2E-400: Phenotype assignment from diplotypes."""

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
# E2E-500: Report Generation & Schema
# ===========================================================================
class TestReportGeneration:
    """E2E-500: Report structure and schema compliance."""

    def test_schema_version(self, pipeline_report):
        """Report should have schema version."""
        assert pipeline_report["schema_version"] == "1.0.0"

    def test_meta_sample_id(self, pipeline_report):
        """Report should identify sample as HG00100."""
        assert pipeline_report["meta"]["sample_id"] == EXPECTED_SAMPLE_ID

    def test_meta_pipeline_version(self, pipeline_report):
        """Report should include pipeline version."""
        assert pipeline_report["meta"]["pipeline_version"] == "1.1.0"

    def test_meta_reference_build(self, pipeline_report):
        """Report should confirm GRCh38."""
        assert pipeline_report["meta"]["reference_build"] == "GRCh38"

    def test_meta_cpic_version(self, pipeline_report):
        """Report should reference CPIC v1.24."""
        assert pipeline_report["meta"]["pgx_db_version"] == "CPIC v1.24"

    def test_meta_gene_tiers(self, pipeline_report):
        """Report should include all requested tiers."""
        tiers = pipeline_report["meta"]["gene_tiers_included"]
        assert "group1" in tiers
        assert "group2" in tiers
        assert "group3" in tiers

    def test_provenance_present(self, pipeline_report):
        """Report should have provenance information."""
        prov = pipeline_report["meta"]["provenance"]
        assert prov["variant_caller"] == "GATK HaplotypeCaller"
        assert prov["pgx_db_hash"] is not None

    def test_confidence_scores_valid(self, pipeline_report):
        """All confidence scores should be between 0.0 and 1.0."""
        for result in pipeline_report["results"]:
            score = result["confidence"]["score"]
            assert 0.0 <= score <= 1.0, (
                f"{result['gene']}: confidence {score} out of range"
            )

    def test_schema_compliance(self, pipeline_report):
        """Report should validate against JSON schema."""
        schema_file = PROJECT_ROOT / "schemas" / "elpgx_report_schema.json"
        if not schema_file.exists():
            pytest.skip("Schema file not found")

        import jsonschema

        with open(schema_file) as f:
            schema = json.load(f)

        jsonschema.validate(pipeline_report, schema)

    def test_report_json_roundtrip(self, pipeline_report):
        """Report should survive JSON serialization roundtrip."""
        serialized = json.dumps(pipeline_report, indent=2)
        deserialized = json.loads(serialized)
        assert deserialized == pipeline_report


# ===========================================================================
# E2E-600: CLI Integration
# ===========================================================================
class TestCLIIntegration:
    """E2E-600: CLI commands work end-to-end."""

    def test_version_command(self):
        """elpgx version should work."""
        result = subprocess.run(
            ["elpgx", "version"],
            capture_output=True,
            text=True,
            timeout=30,
        )
        assert result.returncode == 0
        assert "1.1.0" in result.stdout

    def test_run_with_custom_sample_id(self, tmp_path):
        """Pipeline should accept custom sample ID."""
        output_file = tmp_path / "custom_id_report.json"
        result = subprocess.run(
            [
                "elpgx", "run",
                "-b", str(CRAM_FILE),
                "-r", str(REFERENCE),
                "-o", str(output_file),
                "-s", "CUSTOM_SAMPLE_001",
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
        assert report["meta"]["sample_id"] == "CUSTOM_SAMPLE_001"

    def test_run_group1_only(self, tmp_path):
        """Pipeline should work with group1 genes only."""
        output_file = tmp_path / "group1_report.json"
        result = subprocess.run(
            [
                "elpgx", "run",
                "-b", str(CRAM_FILE),
                "-r", str(REFERENCE),
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
        assert genes == group1_genes, (
            f"Group1 should have exactly {group1_genes}, got {genes}"
        )


# ===========================================================================
# E2E-700: Known-Result Concordance
# ===========================================================================
class TestConcordance:
    """E2E-700: Results match previously validated reference output."""

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
