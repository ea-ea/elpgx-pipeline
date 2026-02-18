"""Unit tests for pgx.phenotype module."""

import json
from pathlib import Path

import pytest

from elpgx.allele_caller import DiplotypeCall
from elpgx.phenotype import (
    PhenotypeResult,
    calculate_activity_score,
    get_functional_status,
    infer_phenotype_from_function,
    load_phenotype_map,
    load_prescribing_implications,
    lookup_phenotype_by_activity_score,
    lookup_phenotype_by_diplotype,
    translate_diplotype_to_phenotype,
)


@pytest.fixture
def sample_diplotype_call():
    """Create a sample diplotype call for testing."""
    return DiplotypeCall(
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


@pytest.fixture
def phenotype_map_dir(tmp_path):
    """Create temporary phenotype map files."""
    phenotype_map = {
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
                "functional_status": "Normal Function / Normal Function",
            },
            {
                "diplotype": "*1/*17",
                "phenotype": "Rapid Metabolizer",
                "activity_score": 2.5,
                "functional_status": "Normal Function / Increased Function",
            },
            {
                "diplotype": "*2/*2",
                "phenotype": "Poor Metabolizer",
                "activity_score": 0.0,
                "functional_status": "No Function / No Function",
            },
        ],
    }

    map_file = tmp_path / "CYP2C19_phenotype.json"
    with open(map_file, "w") as f:
        json.dump(phenotype_map, f)

    return tmp_path


@pytest.fixture
def implications_dir(tmp_path):
    """Create temporary prescribing implications files."""
    implications = {
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
            {
                "phenotype": "Poor Metabolizer",
                "implication": "Use alternative antiplatelet therapy.",
            },
        ],
    }

    impl_file = tmp_path / "CYP2C19_implications.json"
    with open(impl_file, "w") as f:
        json.dump(implications, f)

    return tmp_path


class TestCalculateActivityScore:
    """Tests for calculate_activity_score function."""

    def test_calculate_activity_score(self, sample_diplotype_call):
        """Test activity score calculation."""
        score = calculate_activity_score(sample_diplotype_call)
        assert score == 2.5  # 1.0 + 1.5

    def test_calculate_activity_score_none(self):
        """Test activity score with None values."""
        call = DiplotypeCall(
            gene="TEST",
            allele1="*1",
            allele2="*2",
            diplotype="*1/*2",
            variants=[],
            allele1_function="Normal",
            allele2_function="Decreased",
            allele1_activity=None,
            allele2_activity=None,
        )
        score = calculate_activity_score(call)
        assert score is None


class TestGetFunctionalStatus:
    """Tests for get_functional_status function."""

    def test_different_functions(self, sample_diplotype_call):
        """Test functional status with different functions."""
        status = get_functional_status(sample_diplotype_call)
        assert status == "Normal Function / Increased Function"

    def test_same_functions(self):
        """Test functional status with same functions."""
        call = DiplotypeCall(
            gene="TEST",
            allele1="*1",
            allele2="*1",
            diplotype="*1/*1",
            variants=[],
            allele1_function="Normal Function",
            allele2_function="Normal Function",
            allele1_activity=1.0,
            allele2_activity=1.0,
        )
        status = get_functional_status(call)
        assert status == "Normal Function"


class TestLoadPhenotypeMap:
    """Tests for load_phenotype_map function."""

    def test_load_existing_map(self, phenotype_map_dir):
        """Test loading existing phenotype map."""
        pmap = load_phenotype_map(phenotype_map_dir, "CYP2C19")
        assert pmap is not None
        assert pmap["gene"] == "CYP2C19"
        assert pmap["activity_score_based"] is True

    def test_load_nonexistent_map(self, tmp_path):
        """Test loading non-existent phenotype map."""
        pmap = load_phenotype_map(tmp_path, "NONEXISTENT")
        assert pmap is None


class TestLoadPrescribingImplications:
    """Tests for load_prescribing_implications function."""

    def test_load_existing_implications(self, implications_dir):
        """Test loading existing implications."""
        impl = load_prescribing_implications(implications_dir, "CYP2C19")
        assert "Normal Metabolizer" in impl
        assert "Rapid Metabolizer" in impl

    def test_load_nonexistent_implications(self, tmp_path):
        """Test loading non-existent implications."""
        impl = load_prescribing_implications(tmp_path, "NONEXISTENT")
        assert len(impl) == 0


class TestLookupPhenotypeByDiplotype:
    """Tests for lookup_phenotype_by_diplotype function."""

    def test_lookup_exact_match(self, phenotype_map_dir):
        """Test exact diplotype lookup."""
        pmap = load_phenotype_map(phenotype_map_dir, "CYP2C19")
        phenotype, score, status = lookup_phenotype_by_diplotype("*1/*17", pmap)

        assert phenotype == "Rapid Metabolizer"
        assert score == 2.5

    def test_lookup_reversed_diplotype(self, phenotype_map_dir):
        """Test reversed diplotype lookup (*17/*1 -> *1/*17)."""
        pmap = load_phenotype_map(phenotype_map_dir, "CYP2C19")
        phenotype, score, status = lookup_phenotype_by_diplotype("*17/*1", pmap)

        assert phenotype == "Rapid Metabolizer"

    def test_lookup_no_match(self, phenotype_map_dir):
        """Test diplotype with no match."""
        pmap = load_phenotype_map(phenotype_map_dir, "CYP2C19")
        phenotype, score, status = lookup_phenotype_by_diplotype("*99/*99", pmap)

        assert phenotype is None


class TestLookupPhenotypeByActivityScore:
    """Tests for lookup_phenotype_by_activity_score function."""

    def test_lookup_exact_score(self, phenotype_map_dir):
        """Test exact activity score lookup."""
        pmap = load_phenotype_map(phenotype_map_dir, "CYP2C19")
        phenotype = lookup_phenotype_by_activity_score(2.5, pmap)

        assert phenotype == "Rapid Metabolizer"

    def test_lookup_integer_score(self, phenotype_map_dir):
        """Test integer activity score lookup."""
        pmap = load_phenotype_map(phenotype_map_dir, "CYP2C19")
        phenotype = lookup_phenotype_by_activity_score(2.0, pmap)

        assert phenotype == "Normal Metabolizer"


class TestInferPhenotypeFromFunction:
    """Tests for infer_phenotype_from_function function."""

    def test_infer_poor_metabolizer(self):
        """Test inferring poor metabolizer."""
        phenotype = infer_phenotype_from_function(
            "No Function / No Function", "CYP2C19"
        )
        assert phenotype == "Poor Metabolizer"

    def test_infer_intermediate_metabolizer(self):
        """Test inferring intermediate metabolizer."""
        phenotype = infer_phenotype_from_function(
            "Normal Function / No Function", "CYP2C19"
        )
        assert phenotype == "Intermediate Metabolizer"

    def test_infer_normal_metabolizer(self):
        """Test inferring normal metabolizer."""
        phenotype = infer_phenotype_from_function(
            "Normal Function / Normal Function", "CYP2C19"
        )
        assert phenotype == "Normal Metabolizer"

    def test_infer_rapid_metabolizer(self):
        """Test inferring rapid metabolizer."""
        phenotype = infer_phenotype_from_function(
            "Normal Function / Increased Function", "CYP2C19"
        )
        assert phenotype == "Rapid Metabolizer"


class TestTranslateDiplotypeToPhenotype:
    """Tests for translate_diplotype_to_phenotype function."""

    def test_translate_complete(
        self, sample_diplotype_call, phenotype_map_dir, implications_dir
    ):
        """Test complete translation pipeline."""
        result = translate_diplotype_to_phenotype(
            sample_diplotype_call,
            phenotype_map_dir,
            implications_dir,
        )

        assert isinstance(result, PhenotypeResult)
        assert result.gene == "CYP2C19"
        assert result.diplotype == "*1/*17"
        assert result.phenotype == "Rapid Metabolizer"
        assert "clopidogrel" in result.prescribing_implication.lower()


class TestPhenotypeResult:
    """Tests for PhenotypeResult dataclass."""

    def test_phenotype_result_creation(self):
        """Test creating a PhenotypeResult."""
        result = PhenotypeResult(
            gene="CYP2C19",
            diplotype="*1/*17",
            activity_score=2.5,
            functional_status="Normal Function / Increased Function",
            phenotype="Rapid Metabolizer",
            prescribing_implication="Consider alternative.",
        )

        assert result.gene == "CYP2C19"
        assert result.activity_score == 2.5
        assert result.phenotype == "Rapid Metabolizer"
