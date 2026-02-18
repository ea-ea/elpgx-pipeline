"""Unit tests for elpgx.validate module."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from elpgx.validate import (
    ValidationResult,
    check_file_exists,
    check_index_exists,
    parse_bed_file,
)


@pytest.fixture
def sample_bam_file(tmp_path):
    """Create a sample BAM file (empty but exists)."""
    bam = tmp_path / "sample.bam"
    bam.write_bytes(b"BAM\x01")  # Minimal BAM magic bytes
    return bam


@pytest.fixture
def sample_bed_file(tmp_path):
    """Create a sample BED file with PGx regions."""
    bed_content = """# PGx target regions
chr10\t94760600\t94853000\tCYP2C19
chr10\t94938600\t94989500\tCYP2C9
chr12\t21130000\t21240000\tSLCO1B1
chr16\t88922500\t89004500\tVKORC1
"""
    bed_path = tmp_path / "pgx_regions.bed"
    bed_path.write_text(bed_content)
    return bed_path


class TestValidationResult:
    """Tests for ValidationResult dataclass."""

    def test_valid_result(self):
        """Test creating a valid result."""
        result = ValidationResult(
            is_valid=True,
            sample_id="SAMPLE_001",
            errors=[],
            warnings=[],
        )

        assert result.is_valid is True
        assert result.sample_id == "SAMPLE_001"
        assert len(result.errors) == 0

    def test_invalid_result_with_errors(self):
        """Test creating an invalid result with errors."""
        result = ValidationResult(
            is_valid=False,
            errors=["File not found", "Index missing"],
        )

        assert result.is_valid is False
        assert len(result.errors) == 2

    def test_result_with_warnings(self):
        """Test creating a result with warnings."""
        result = ValidationResult(
            is_valid=True,
            sample_id="SAMPLE_001",
            warnings=["Low coverage for CYP2C19: 5x"],
        )

        assert result.is_valid is True
        assert len(result.warnings) == 1

    def test_result_with_coverage(self):
        """Test creating a result with coverage stats."""
        coverage = {"CYP2C19": 45.5, "CYP2C9": 38.2}
        result = ValidationResult(
            is_valid=True,
            sample_id="SAMPLE_001",
            coverage_stats=coverage,
        )

        assert result.coverage_stats is not None
        assert result.coverage_stats["CYP2C19"] == 45.5


class TestCheckFileExists:
    """Tests for check_file_exists function."""

    def test_existing_file(self, sample_bam_file):
        """Test checking an existing file."""
        ok, err = check_file_exists(sample_bam_file)

        assert ok is True
        assert err is None

    def test_nonexistent_file(self, tmp_path):
        """Test checking a non-existent file."""
        fake_path = tmp_path / "nonexistent.bam"
        ok, err = check_file_exists(fake_path)

        assert ok is False
        assert "does not exist" in err

    def test_directory_not_file(self, tmp_path):
        """Test checking a directory instead of file."""
        ok, err = check_file_exists(tmp_path)

        assert ok is False
        assert "not a file" in err

    def test_empty_file(self, tmp_path):
        """Test checking an empty file."""
        empty_file = tmp_path / "empty.bam"
        empty_file.touch()

        ok, err = check_file_exists(empty_file)

        assert ok is False
        assert "empty" in err


class TestCheckIndexExists:
    """Tests for check_index_exists function."""

    def test_bam_with_bai_index(self, tmp_path):
        """Test BAM file with .bam.bai index."""
        bam = tmp_path / "sample.bam"
        bam.touch()
        index = tmp_path / "sample.bam.bai"
        index.touch()

        ok, err = check_index_exists(bam)

        assert ok is True
        assert err is None

    def test_bam_with_bai_suffix(self, tmp_path):
        """Test BAM file with .bai suffix index."""
        bam = tmp_path / "sample.bam"
        bam.touch()
        index = tmp_path / "sample.bai"
        index.touch()

        ok, err = check_index_exists(bam)

        assert ok is True

    def test_cram_with_crai_index(self, tmp_path):
        """Test CRAM file with .cram.crai index."""
        cram = tmp_path / "sample.cram"
        cram.touch()
        index = tmp_path / "sample.cram.crai"
        index.touch()

        ok, err = check_index_exists(cram)

        assert ok is True

    def test_missing_index(self, tmp_path):
        """Test BAM file without index."""
        bam = tmp_path / "sample.bam"
        bam.touch()

        ok, err = check_index_exists(bam)

        assert ok is False
        assert "Index file not found" in err

    def test_unsupported_format(self, tmp_path):
        """Test unsupported file format."""
        sam = tmp_path / "sample.sam"
        sam.touch()

        ok, err = check_index_exists(sam)

        assert ok is False
        assert "Unsupported file format" in err


class TestParseBedFile:
    """Tests for parse_bed_file function."""

    def test_parse_valid_bed(self, sample_bed_file):
        """Test parsing a valid BED file."""
        regions = parse_bed_file(sample_bed_file)

        assert len(regions) == 4

        # Check first region
        assert regions[0]["chrom"] == "chr10"
        assert regions[0]["start"] == 94760600
        assert regions[0]["end"] == 94853000
        assert regions[0]["name"] == "CYP2C19"

    def test_parse_bed_with_comments(self, tmp_path):
        """Test parsing BED with comment lines."""
        bed_content = """# Header comment
# Another comment
chr10\t100\t200\tgene1
chr10\t300\t400\tgene2
"""
        bed_path = tmp_path / "test.bed"
        bed_path.write_text(bed_content)

        regions = parse_bed_file(bed_path)

        assert len(regions) == 2

    def test_parse_bed_without_names(self, tmp_path):
        """Test parsing BED without 4th column."""
        bed_content = """chr10\t100\t200
chr10\t300\t400
"""
        bed_path = tmp_path / "test.bed"
        bed_path.write_text(bed_content)

        regions = parse_bed_file(bed_path)

        assert len(regions) == 2
        assert regions[0]["name"] is None

    def test_parse_empty_bed(self, tmp_path):
        """Test parsing empty BED file."""
        bed_path = tmp_path / "empty.bed"
        bed_path.write_text("")

        regions = parse_bed_file(bed_path)

        assert len(regions) == 0

    def test_parse_bed_skip_malformed(self, tmp_path):
        """Test parsing BED with malformed lines."""
        bed_content = """chr10\t100\t200\tgene1
invalid_line
chr10\t300
chr10\t400\t500\tgene2
"""
        bed_path = tmp_path / "test.bed"
        bed_path.write_text(bed_content)

        regions = parse_bed_file(bed_path)

        # Should only parse lines with at least 3 columns
        assert len(regions) == 2


class TestCheckReferenceBuild:
    """Tests for check_reference_build function (mocked)."""

    @patch("elpgx.validate.pysam.AlignmentFile")
    def test_valid_hg38_reference(self, mock_alignment_file):
        """Test validation of valid hg38 reference."""
        from elpgx.validate import check_reference_build

        # Mock BAM header with hg38 chromosomes (need 3+ to pass threshold)
        mock_bam = MagicMock()
        mock_bam.__enter__ = MagicMock(return_value=mock_bam)
        mock_bam.__exit__ = MagicMock(return_value=False)
        mock_bam.header.to_dict.return_value = {
            "SQ": [
                {"SN": "chr1", "LN": 248956422},
                {"SN": "chr2", "LN": 242193529},
                {"SN": "chr10", "LN": 133797422},
            ]
        }
        mock_alignment_file.return_value = mock_bam

        ok, err, chrom_style, detected_build = check_reference_build(Path("/fake/path.bam"))

        assert ok is True
        assert err is None
        assert detected_build == "GRCh38"

    @patch("elpgx.validate.pysam.AlignmentFile")
    def test_non_hg38_reference(self, mock_alignment_file):
        """Test validation rejects reference with too few matching chroms for any build."""
        from elpgx.validate import check_reference_build

        mock_bam = MagicMock()
        mock_bam.__enter__ = MagicMock(return_value=mock_bam)
        mock_bam.__exit__ = MagicMock(return_value=False)
        mock_bam.header.to_dict.return_value = {
            "SQ": [
                {"SN": "1", "LN": 248956422},  # Ensembl style
                {"SN": "10", "LN": 133797422},
            ]
        }
        mock_alignment_file.return_value = mock_bam

        ok, err, chrom_style, detected_build = check_reference_build(Path("/fake/path.bam"))

        # Only 2 matching chroms, needs 3 minimum
        assert ok is False

    @patch("elpgx.validate.pysam.AlignmentFile")
    def test_valid_hg19_reference(self, mock_alignment_file):
        """Test validation of valid hg19/GRCh37 reference."""
        from elpgx.validate import check_reference_build

        mock_bam = MagicMock()
        mock_bam.__enter__ = MagicMock(return_value=mock_bam)
        mock_bam.__exit__ = MagicMock(return_value=False)
        mock_bam.header.to_dict.return_value = {
            "SQ": [
                {"SN": "1", "LN": 249250621},
                {"SN": "2", "LN": 243199373},
                {"SN": "6", "LN": 171115067},
                {"SN": "10", "LN": 135534747},
            ]
        }
        mock_alignment_file.return_value = mock_bam

        ok, err, chrom_style, detected_build = check_reference_build(Path("/fake/path.bam"))

        assert ok is True
        assert err is None
        assert detected_build == "GRCh37"
        assert chrom_style == "ensembl"

    @patch("elpgx.validate.pysam.AlignmentFile")
    def test_wrong_chromosome_length(self, mock_alignment_file):
        """Test validation detects unknown chromosome lengths."""
        from elpgx.validate import check_reference_build

        mock_bam = MagicMock()
        mock_bam.__enter__ = MagicMock(return_value=mock_bam)
        mock_bam.__exit__ = MagicMock(return_value=False)
        mock_bam.header.to_dict.return_value = {
            "SQ": [
                {"SN": "chr1", "LN": 100000000},  # Invalid length
            ]
        }
        mock_alignment_file.return_value = mock_bam

        ok, err, chrom_style, detected_build = check_reference_build(Path("/fake/path.bam"))

        assert ok is False


class TestExtractSampleId:
    """Tests for extract_sample_id function (mocked)."""

    @patch("elpgx.validate.pysam.AlignmentFile")
    def test_extract_sample_id(self, mock_alignment_file):
        """Test extracting sample ID from @RG header."""
        from elpgx.validate import extract_sample_id

        mock_bam = MagicMock()
        mock_bam.__enter__ = MagicMock(return_value=mock_bam)
        mock_bam.__exit__ = MagicMock(return_value=False)
        mock_bam.header.to_dict.return_value = {
            "RG": [
                {"ID": "rg1", "SM": "SAMPLE_001", "PL": "ILLUMINA"},
            ]
        }
        mock_alignment_file.return_value = mock_bam

        sample_id = extract_sample_id(Path("/fake/path.bam"))

        assert sample_id == "SAMPLE_001"

    @patch("elpgx.validate.pysam.AlignmentFile")
    def test_no_rg_header(self, mock_alignment_file):
        """Test when no @RG header present."""
        from elpgx.validate import extract_sample_id

        mock_bam = MagicMock()
        mock_bam.__enter__ = MagicMock(return_value=mock_bam)
        mock_bam.__exit__ = MagicMock(return_value=False)
        mock_bam.header.to_dict.return_value = {"SQ": []}
        mock_alignment_file.return_value = mock_bam

        sample_id = extract_sample_id(Path("/fake/path.bam"))

        assert sample_id is None


class TestValidateAlignmentFile:
    """Tests for validate_alignment_file function (mocked)."""

    @patch("elpgx.validate.pysam.AlignmentFile")
    @patch("elpgx.validate.check_file_exists")
    @patch("elpgx.validate.check_index_exists")
    def test_valid_file(
        self, mock_check_index, mock_check_file, mock_alignment_file, tmp_path
    ):
        """Test validation of a valid file."""
        from elpgx.validate import validate_alignment_file

        mock_check_file.return_value = (True, None)
        mock_check_index.return_value = (True, None)

        mock_bam = MagicMock()
        mock_bam.__enter__ = MagicMock(return_value=mock_bam)
        mock_bam.__exit__ = MagicMock(return_value=False)
        mock_bam.header.to_dict.return_value = {
            "SQ": [
                {"SN": "chr1", "LN": 248956422},
                {"SN": "chr2", "LN": 242193529},
                {"SN": "chr10", "LN": 133797422},
            ],
            "RG": [{"SM": "SAMPLE_001"}],
        }
        mock_alignment_file.return_value = mock_bam

        bam_path = tmp_path / "sample.bam"

        result = validate_alignment_file(bam_path)

        assert result.is_valid is True
        assert result.sample_id == "SAMPLE_001"

    @patch("elpgx.validate.check_file_exists")
    def test_missing_file(self, mock_check_file, tmp_path):
        """Test validation when file is missing."""
        from elpgx.validate import validate_alignment_file

        mock_check_file.return_value = (False, "File does not exist")

        bam_path = tmp_path / "missing.bam"

        result = validate_alignment_file(bam_path)

        assert result.is_valid is False
        assert "File does not exist" in result.errors

    @patch("elpgx.validate.pysam.AlignmentFile")
    @patch("elpgx.validate.check_file_exists")
    @patch("elpgx.validate.check_index_exists")
    def test_missing_index(
        self, mock_check_index, mock_check_file, mock_alignment_file, tmp_path
    ):
        """Test validation when index is missing."""
        from elpgx.validate import validate_alignment_file

        mock_check_file.return_value = (True, None)
        mock_check_index.return_value = (False, "Index file not found")

        mock_bam = MagicMock()
        mock_bam.__enter__ = MagicMock(return_value=mock_bam)
        mock_bam.__exit__ = MagicMock(return_value=False)
        mock_bam.header.to_dict.return_value = {
            "SQ": [{"SN": "chr10", "LN": 133797422}],
            "RG": [{"SM": "SAMPLE_001"}],
        }
        mock_alignment_file.return_value = mock_bam

        bam_path = tmp_path / "sample.bam"

        result = validate_alignment_file(bam_path)

        assert result.is_valid is False
        assert any("Index" in e for e in result.errors)
