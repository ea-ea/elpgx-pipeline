"""
Input validation module for PGx Pipeline.

Validates BAM/CRAM files before pipeline processing:
- File existence and readability
- Index file presence
- hg38/GRCh38 reference verification
- Sample ID extraction
- Coverage assessment over PGx regions
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pysam


# Expected chromosome lengths for hg38/GRCh38 (subset for validation)
# Supports both UCSC (chr1) and Ensembl/NCBI (1) naming conventions
HG38_CHROM_LENGTHS: dict[str, int] = {
    # UCSC style (chr prefix)
    "chr1": 248956422,
    "chr2": 242193529,
    "chr6": 170805979,
    "chr10": 133797422,
    "chr12": 133275309,
    "chr13": 114364328,
    "chr16": 90338345,
    "chr22": 50818468,
    # Ensembl/NCBI style (no chr prefix)
    "1": 248956422,
    "2": 242193529,
    "6": 170805979,
    "10": 133797422,
    "12": 133275309,
    "13": 114364328,
    "16": 90338345,
    "22": 50818468,
}

# Expected chromosome lengths for hg19/GRCh37 (subset for validation)
HG37_CHROM_LENGTHS: dict[str, int] = {
    # UCSC style (chr prefix)
    "chr1": 249250621,
    "chr2": 243199373,
    "chr6": 171115067,
    "chr10": 135534747,
    "chr12": 133851895,
    "chr13": 115169878,
    "chr16": 90354753,
    "chr22": 51304566,
    # Ensembl/NCBI style (no chr prefix)
    "1": 249250621,
    "2": 243199373,
    "6": 171115067,
    "10": 135534747,
    "12": 133851895,
    "13": 115169878,
    "16": 90354753,
    "22": 51304566,
}


@dataclass
class ValidationResult:
    """Result of input file validation."""

    is_valid: bool
    sample_id: str | None = None
    reference_build: str | None = None
    chrom_style: str | None = None  # "ucsc" (chr1) or "ensembl" (1)
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    coverage_stats: dict[str, float] | None = None


def check_file_exists(path: Path) -> tuple[bool, str | None]:
    """
    Check if file exists and is readable.

    Args:
        path: Path to the file to check.

    Returns:
        Tuple of (success, error_message). error_message is None if success.
    """
    if not path.exists():
        return False, f"File does not exist: {path}"
    if not path.is_file():
        return False, f"Path is not a file: {path}"
    if not path.stat().st_size > 0:
        return False, f"File is empty: {path}"
    return True, None


def check_index_exists(bam_path: Path) -> tuple[bool, str | None]:
    """
    Check if BAM/CRAM index file exists.

    Looks for .bai, .crai, or .bam.bai / .cram.crai files.

    Args:
        bam_path: Path to the BAM/CRAM file.

    Returns:
        Tuple of (success, error_message). error_message is None if success.
    """
    suffix = bam_path.suffix.lower()

    if suffix == ".bam":
        index_paths = [
            bam_path.with_suffix(".bam.bai"),
            bam_path.with_suffix(".bai"),
        ]
    elif suffix == ".cram":
        index_paths = [
            bam_path.with_suffix(".cram.crai"),
            bam_path.with_suffix(".crai"),
        ]
    else:
        return False, f"Unsupported file format: {suffix}. Expected .bam or .cram"

    for idx_path in index_paths:
        if idx_path.exists():
            return True, None

    expected = " or ".join(str(p.name) for p in index_paths)
    return False, f"Index file not found. Expected: {expected}"


def check_reference_build(bam_path: Path) -> tuple[bool, str | None, str | None, str | None]:
    """
    Detect BAM/CRAM reference genome build (GRCh38 or GRCh37).

    Checks @SQ header lines and validates chromosome lengths against known values.
    Supports both UCSC (chr1, chr2...) and Ensembl/NCBI (1, 2...) naming conventions.

    Args:
        bam_path: Path to the BAM/CRAM file.

    Returns:
        Tuple of (is_valid, error_message, chrom_style, detected_build).
        - is_valid: True if reference is GRCh38 or GRCh37
        - error_message: None if valid, else error description
        - chrom_style: "ucsc" (chr prefix) or "ensembl" (no prefix) or None
        - detected_build: "GRCh38" or "GRCh37" or None
    """
    try:
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            header = bam.header.to_dict()
    except Exception as e:
        return False, f"Failed to read BAM header: {e}", None, None

    sq_records = header.get("SQ", [])
    if not sq_records:
        return False, "No @SQ records found in BAM header", None, None

    # Determine chromosome naming style
    chr_names = [rec.get("SN", "") for rec in sq_records]
    has_chr_prefix = any(name.startswith("chr") for name in chr_names)
    has_numeric = any(name in ["1", "2", "3", "10", "22", "X", "Y"] for name in chr_names)

    if has_chr_prefix:
        chrom_style = "ucsc"
    elif has_numeric:
        chrom_style = "ensembl"
    else:
        return False, "Could not determine chromosome naming convention", None, None

    # Try GRCh38 first
    matched_38 = 0
    for rec in sq_records:
        chrom = rec.get("SN", "")
        length = rec.get("LN", 0)
        if chrom in HG38_CHROM_LENGTHS and length == HG38_CHROM_LENGTHS[chrom]:
            matched_38 += 1

    if matched_38 >= 3:
        return True, None, chrom_style, "GRCh38"

    # Try GRCh37
    matched_37 = 0
    for rec in sq_records:
        chrom = rec.get("SN", "")
        length = rec.get("LN", 0)
        if chrom in HG37_CHROM_LENGTHS and length == HG37_CHROM_LENGTHS[chrom]:
            matched_37 += 1

    if matched_37 >= 3:
        return True, None, chrom_style, "GRCh37"

    return False, (
        f"Could not confirm reference build (matched {matched_38} GRCh38, "
        f"{matched_37} GRCh37 chromosomes; need at least 3)"
    ), chrom_style, None


def extract_sample_id(bam_path: Path) -> str | None:
    """
    Extract sample ID from BAM @RG SM tag.

    Args:
        bam_path: Path to the BAM/CRAM file.

    Returns:
        Sample ID string, or None if not found.
    """
    try:
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            header = bam.header.to_dict()
    except Exception:
        return None

    rg_records = header.get("RG", [])
    if not rg_records:
        return None

    # Return first sample ID found
    for rg in rg_records:
        sample_id = rg.get("SM")
        if sample_id:
            return sample_id

    return None


def parse_bed_file(bed_path: Path) -> list[dict[str, Any]]:
    """
    Parse BED file into list of regions.

    Args:
        bed_path: Path to BED file.

    Returns:
        List of dicts with keys: chrom, start, end, name (optional).
    """
    regions = []
    with open(bed_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            region = {
                "chrom": parts[0],
                "start": int(parts[1]),
                "end": int(parts[2]),
                "name": parts[3] if len(parts) > 3 else None,
            }
            regions.append(region)
    return regions


def calculate_region_coverage(
    bam_path: Path,
    regions_bed: Path,
    min_mapping_quality: int = 20,
) -> dict[str, float]:
    """
    Calculate mean coverage over PGx target regions.

    Args:
        bam_path: Path to BAM/CRAM file.
        regions_bed: Path to BED file with target regions.
        min_mapping_quality: Minimum mapping quality to count.

    Returns:
        Dict mapping gene name (or region) to mean depth.
    """
    regions = parse_bed_file(regions_bed)
    coverage_stats: dict[str, float] = {}

    try:
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            for region in regions:
                chrom = region["chrom"]
                start = region["start"]
                end = region["end"]
                name = region["name"] or f"{chrom}:{start}-{end}"

                try:
                    # Get pileup coverage over region
                    depths = []
                    for pileup_col in bam.pileup(
                        chrom,
                        start,
                        end,
                        min_mapping_quality=min_mapping_quality,
                        truncate=True,
                    ):
                        depths.append(pileup_col.n)

                    if depths:
                        mean_depth = sum(depths) / len(depths)
                    else:
                        mean_depth = 0.0

                    coverage_stats[name] = round(mean_depth, 2)

                except ValueError:
                    # Region not in BAM (chromosome not found)
                    coverage_stats[name] = 0.0

    except Exception:
        # Return empty stats on error
        return {}

    return coverage_stats


def validate_alignment_file(
    bam_path: Path,
    reference_path: Path | None = None,
    pgx_regions_bed: Path | None = None,
    min_mapping_quality: int = 20,
) -> ValidationResult:
    """
    Validate a BAM/CRAM file for PGx analysis.

    Performs the following checks:
    1. File exists and is readable
    2. Index file exists (.bai/.crai)
    3. Header confirms hg38/GRCh38 reference
    4. Sample ID extracted from @RG header
    5. (Optional) Coverage over PGx regions

    Args:
        bam_path: Path to BAM/CRAM file.
        reference_path: Path to reference FASTA (for CRAM, optional).
        pgx_regions_bed: Path to BED file with PGx regions (optional).
        min_mapping_quality: Minimum mapping quality for coverage calc.

    Returns:
        ValidationResult with is_valid, sample_id, errors, warnings.
    """
    errors: list[str] = []
    warnings: list[str] = []
    sample_id: str | None = None
    chrom_style: str | None = None
    coverage_stats: dict[str, float] | None = None

    # 1. Check file exists
    ok, err = check_file_exists(bam_path)
    if not ok:
        errors.append(err)  # type: ignore
        return ValidationResult(
            is_valid=False,
            errors=errors,
            warnings=warnings,
        )

    # 2. Check index exists
    ok, err = check_index_exists(bam_path)
    if not ok:
        errors.append(err)  # type: ignore

    # 3. Check reference build (supports both GRCh38 and GRCh37)
    ok, err, chrom_style, detected_build = check_reference_build(bam_path)
    if not ok:
        errors.append(err)  # type: ignore
    else:
        if chrom_style == "ensembl":
            warnings.append(
                "BAM uses Ensembl/NCBI chromosome naming (1, 2, 3...). "
                "Pipeline will handle conversion if needed."
            )
        if detected_build == "GRCh37":
            warnings.append(
                "BAM is aligned to GRCh37/hg19. Pipeline will translate "
                "coordinates to GRCh38 for allele matching."
            )

    # 4. Extract sample ID
    sample_id = extract_sample_id(bam_path)
    if sample_id is None:
        warnings.append("Could not extract sample ID from @RG header")
        sample_id = bam_path.stem  # Fallback to filename

    # 5. Calculate coverage (if BED provided)
    if pgx_regions_bed is not None:
        bed_ok, bed_err = check_file_exists(pgx_regions_bed)
        if bed_ok:
            coverage_stats = calculate_region_coverage(
                bam_path,
                pgx_regions_bed,
                min_mapping_quality,
            )
            # Warn about low coverage regions
            for gene, depth in coverage_stats.items():
                if depth < 10:
                    warnings.append(f"Low coverage for {gene}: {depth}x")
        else:
            warnings.append(f"Could not read BED file: {bed_err}")

    is_valid = len(errors) == 0

    return ValidationResult(
        is_valid=is_valid,
        sample_id=sample_id,
        reference_build=detected_build if is_valid else None,
        chrom_style=chrom_style if is_valid else None,
        errors=errors,
        warnings=warnings,
        coverage_stats=coverage_stats,
    )


def validate_bam_file(
    bam_path: Path,
    reference_path: Path | None = None,
    pgx_regions_bed: Path | None = None,
) -> ValidationResult:
    """
    Validate a BAM/CRAM file for PGx analysis.

    This is a convenience wrapper around validate_alignment_file.

    Args:
        bam_path: Path to BAM/CRAM file.
        reference_path: Path to reference FASTA (required for variant calling).
        pgx_regions_bed: Path to BED file with PGx regions (optional).

    Returns:
        ValidationResult with is_valid, sample_id, reference_build, errors, warnings.
    """
    result = validate_alignment_file(
        bam_path=bam_path,
        reference_path=reference_path,
        pgx_regions_bed=pgx_regions_bed,
    )

    # Validate reference if provided
    if reference_path is not None:
        ref_ok, ref_err = validate_reference(reference_path)
        if not ref_ok:
            result.errors.append(f"Reference: {ref_err}")
            result.is_valid = False

    return result


def is_bgzipped(file_path: Path) -> bool:
    """
    Check if a file is bgzipped (block gzip) vs regular gzip.

    BGZip files have a specific magic number and block structure.
    Regular gzip: 1f 8b
    BGZip: 1f 8b 08 04 ... with BGZF extra field

    Args:
        file_path: Path to the gzipped file.

    Returns:
        True if file is bgzipped, False if regular gzip or other.
    """
    try:
        with open(file_path, "rb") as f:
            # Read first 18 bytes to check for BGZF signature
            header = f.read(18)
            if len(header) < 18:
                return False
            # Check gzip magic (1f 8b), compression method (08), flags (04 = extra field)
            if header[0:2] != b'\x1f\x8b':
                return False
            if header[2] != 8:  # deflate compression
                return False
            if not (header[3] & 0x04):  # FEXTRA flag not set
                return False
            # Check for BGZF extra subfield (SI1=66='B', SI2=67='C')
            # Extra field starts at byte 10, length at bytes 10-11
            xlen = header[10] + (header[11] << 8)
            if xlen < 6:
                return False
            # Read extra field
            extra = header[12:12 + min(xlen, 6)]
            if len(extra) >= 4 and extra[0:2] == b'BC':
                return True
            return False
    except Exception:
        return False


def validate_reference(reference_path: Path) -> tuple[bool, str | None]:
    """
    Validate reference genome file.

    Supports both plain (.fa, .fasta) and bgzipped (.fa.gz, .fasta.gz) FASTA files.
    Note: Regular gzip is NOT supported - must use bgzip for compressed files.

    Args:
        reference_path: Path to reference FASTA.

    Returns:
        Tuple of (is_valid, error_message).
    """
    ok, err = check_file_exists(reference_path)
    if not ok:
        return False, err

    ref_str = str(reference_path)
    is_gzipped = ref_str.endswith(".gz")

    # Check for .fai index
    fai_path = Path(ref_str + ".fai")
    if not fai_path.exists():
        return False, (
            f"Reference index not found: {fai_path.name}. "
            f"Run: samtools faidx {reference_path.name}"
        )

    # For gzipped files, verify it's bgzipped (not regular gzip)
    if is_gzipped:
        if not is_bgzipped(reference_path):
            return False, (
                f"Reference file is compressed with gzip, but bgzip is required.\n"
                f"Please decompress and optionally re-compress with bgzip:\n"
                f"  gunzip {reference_path.name}\n"
                f"  samtools faidx {reference_path.stem}\n"
                f"Or to re-compress with bgzip:\n"
                f"  gunzip {reference_path.name} && bgzip {reference_path.stem}"
            )
        # Check for .gzi index
        gzi_path = Path(ref_str + ".gzi")
        if not gzi_path.exists():
            return False, (
                f"BGZip index not found: {gzi_path.name}. "
                f"Run: samtools faidx {reference_path.name}"
            )

    return True, None
