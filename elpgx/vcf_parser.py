"""
VCF parsing module for PGx Pipeline.

Parses VCF files and extracts variants relevant to pharmacogenetic analysis.
Supports both cyvcf2 and pysam backends.

Also provides variant calling from BAM/CRAM files using GATK HaplotypeCaller.
"""

import json
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

from elpgx.models import Variant


@dataclass
class ParsedVariant:
    """Parsed variant from VCF before conversion to model."""

    chrom: str
    pos: int
    ref: str
    alt: str
    rsid: str | None
    genotype: str  # "0/0", "0/1", "1/1", "./."
    quality: float
    depth: int
    genotype_quality: int


def _try_import_cyvcf2() -> Any:
    """Try to import cyvcf2, return module or None."""
    try:
        import cyvcf2

        return cyvcf2
    except ImportError:
        return None


def _try_import_pysam() -> Any:
    """Try to import pysam.VariantFile, return module or None."""
    try:
        import pysam

        return pysam
    except ImportError:
        return None


def extract_genotype_cyvcf2(variant: Any, sample_idx: int = 0) -> str:
    """
    Extract genotype string from cyvcf2 variant.

    Args:
        variant: cyvcf2 Variant object.
        sample_idx: Sample index (default 0 for single-sample VCF).

    Returns:
        Genotype string (e.g., "0/1", "1/1", "./.").
    """
    gt = variant.genotypes[sample_idx]
    allele1, allele2, phased = gt[0], gt[1], gt[2]

    if allele1 == -1 or allele2 == -1:
        return "./."

    sep = "|" if phased else "/"
    return f"{allele1}{sep}{allele2}"


def extract_genotype_pysam(record: Any, sample_idx: int = 0) -> str:
    """
    Extract genotype string from pysam VariantRecord.

    Args:
        record: pysam VariantRecord object.
        sample_idx: Sample index (default 0).

    Returns:
        Genotype string (e.g., "0/1", "1/1", "./.").
    """
    sample = record.samples[sample_idx]
    gt = sample.get("GT")

    if gt is None:
        return "./."

    alleles = list(gt)
    if None in alleles:
        return "./."

    # Check if phased
    phased = sample.phased if hasattr(sample, "phased") else False
    sep = "|" if phased else "/"

    return sep.join(str(a) for a in alleles)


def extract_depth_cyvcf2(variant: Any, sample_idx: int = 0) -> int:
    """
    Extract read depth from cyvcf2 variant.

    Tries DP from FORMAT, then from INFO.

    Args:
        variant: cyvcf2 Variant object.
        sample_idx: Sample index.

    Returns:
        Read depth, or 0 if not found.
    """
    # Try FORMAT DP
    try:
        dp = variant.format("DP")
        if dp is not None and len(dp) > sample_idx:
            val = dp[sample_idx]
            if val is not None and val >= 0:
                return int(val)
    except (KeyError, TypeError):
        pass

    # Try INFO DP
    try:
        dp = variant.INFO.get("DP")
        if dp is not None:
            return int(dp)
    except (KeyError, TypeError, AttributeError):
        pass

    # Try AD (allelic depths) sum
    try:
        ad = variant.format("AD")
        if ad is not None and len(ad) > sample_idx:
            return int(sum(ad[sample_idx]))
    except (KeyError, TypeError):
        pass

    return 0


def extract_depth_pysam(record: Any, sample_idx: int = 0) -> int:
    """
    Extract read depth from pysam VariantRecord.

    Args:
        record: pysam VariantRecord object.
        sample_idx: Sample index.

    Returns:
        Read depth, or 0 if not found.
    """
    sample = record.samples[sample_idx]

    # Try FORMAT DP
    dp = sample.get("DP")
    if dp is not None:
        return int(dp)

    # Try INFO DP
    dp = record.info.get("DP")
    if dp is not None:
        return int(dp)

    # Try AD sum
    ad = sample.get("AD")
    if ad is not None:
        return sum(int(x) for x in ad if x is not None)

    return 0


def extract_gq_cyvcf2(variant: Any, sample_idx: int = 0) -> int:
    """Extract genotype quality from cyvcf2 variant."""
    try:
        gq = variant.format("GQ")
        if gq is not None and len(gq) > sample_idx:
            val = gq[sample_idx]
            if val is not None and val >= 0:
                return int(val)
    except (KeyError, TypeError):
        pass
    return 0


def extract_gq_pysam(record: Any, sample_idx: int = 0) -> int:
    """Extract genotype quality from pysam VariantRecord."""
    sample = record.samples[sample_idx]
    gq = sample.get("GQ")
    if gq is not None:
        return int(gq)
    return 0


def parse_vcf_cyvcf2(vcf_path: Path) -> list[ParsedVariant]:
    """
    Parse VCF using cyvcf2.

    Args:
        vcf_path: Path to VCF file.

    Returns:
        List of ParsedVariant objects.
    """
    cyvcf2 = _try_import_cyvcf2()
    if cyvcf2 is None:
        raise ImportError("cyvcf2 is required for VCF parsing")

    variants = []
    vcf = cyvcf2.VCF(str(vcf_path))

    for record in vcf:
        # Handle multi-allelic sites
        for alt_idx, alt in enumerate(record.ALT):
            rsid = record.ID if record.ID and record.ID != "." else None

            parsed = ParsedVariant(
                chrom=record.CHROM,
                pos=record.POS,
                ref=record.REF,
                alt=alt,
                rsid=rsid,
                genotype=extract_genotype_cyvcf2(record),
                quality=record.QUAL if record.QUAL is not None else 0.0,
                depth=extract_depth_cyvcf2(record),
                genotype_quality=extract_gq_cyvcf2(record),
            )
            variants.append(parsed)

    vcf.close()
    return variants


def parse_vcf_pysam(vcf_path: Path) -> list[ParsedVariant]:
    """
    Parse VCF using pysam.

    Args:
        vcf_path: Path to VCF file.

    Returns:
        List of ParsedVariant objects.
    """
    pysam = _try_import_pysam()
    if pysam is None:
        raise ImportError("pysam is required for VCF parsing")

    variants = []

    with pysam.VariantFile(str(vcf_path)) as vcf:
        for record in vcf:
            for alt_idx, alt in enumerate(record.alts or []):
                rsid = record.id if record.id and record.id != "." else None

                parsed = ParsedVariant(
                    chrom=record.chrom,
                    pos=record.pos,
                    ref=record.ref,
                    alt=alt,
                    rsid=rsid,
                    genotype=extract_genotype_pysam(record),
                    quality=record.qual if record.qual is not None else 0.0,
                    depth=extract_depth_pysam(record),
                    genotype_quality=extract_gq_pysam(record),
                )
                variants.append(parsed)

    return variants


def parse_vcf(vcf_path: Path) -> list[ParsedVariant]:
    """
    Parse VCF file and return all variants.

    Automatically selects cyvcf2 or pysam backend based on availability.

    Args:
        vcf_path: Path to VCF file (.vcf or .vcf.gz).

    Returns:
        List of ParsedVariant objects.

    Raises:
        ImportError: If neither cyvcf2 nor pysam is available.
        FileNotFoundError: If VCF file does not exist.
    """
    if not vcf_path.exists():
        raise FileNotFoundError(f"VCF file not found: {vcf_path}")

    # Prefer cyvcf2 for performance
    if _try_import_cyvcf2() is not None:
        return parse_vcf_cyvcf2(vcf_path)
    elif _try_import_pysam() is not None:
        return parse_vcf_pysam(vcf_path)
    else:
        raise ImportError("Either cyvcf2 or pysam is required for VCF parsing")


def load_grch37_coordinate_map(data_dir: Path) -> dict[str, dict[int, tuple[str, int]]]:
    """
    Load GRCh37-to-GRCh38 coordinate mapping.

    The mapping file maps GRCh38 positions to GRCh37. This function inverts it
    to provide GRCh37-to-GRCh38 translation keyed by GRCh37 chrom and position.

    Args:
        data_dir: Path to data directory containing grch37_coordinate_map.json.

    Returns:
        Dict mapping GRCh37 chrom -> {GRCh37 pos -> (GRCh38 chrom, GRCh38 pos)}.
    """
    map_file = data_dir / "grch37_coordinate_map.json"
    if not map_file.exists():
        return {}

    with open(map_file) as f:
        raw_map = json.load(f)

    # Build inverted map: grch37_chrom -> {grch37_pos -> (grch38_chrom, grch38_pos)}
    coord_map: dict[str, dict[int, tuple[str, int]]] = {}

    # Also need the GRCh38 chromosome from allele definitions
    for gene_name, positions in raw_map.items():
        for grch38_pos_str, grch37_info in positions.items():
            grch37_chrom = grch37_info["grch37_chrom"]
            grch37_pos = grch37_info["grch37_pos"]
            grch38_pos = int(grch38_pos_str)
            # GRCh38 chrom has "chr" prefix
            grch38_chrom = f"chr{grch37_chrom}"

            if grch37_chrom not in coord_map:
                coord_map[grch37_chrom] = {}
            coord_map[grch37_chrom][grch37_pos] = (grch38_chrom, grch38_pos)

    return coord_map


def translate_variants_grch37_to_grch38(
    variants: list[ParsedVariant],
    coord_map: dict[str, dict[int, tuple[str, int]]],
) -> list[ParsedVariant]:
    """
    Translate variant coordinates from GRCh37 to GRCh38.

    Creates new ParsedVariant objects with GRCh38 coordinates for variants
    that have a known mapping. Variants without a mapping are dropped.

    Args:
        variants: List of variants with GRCh37 coordinates.
        coord_map: GRCh37-to-GRCh38 coordinate mapping.

    Returns:
        List of variants with GRCh38 coordinates.
    """
    translated = []
    for v in variants:
        chrom = v.chrom
        if chrom in coord_map and v.pos in coord_map[chrom]:
            grch38_chrom, grch38_pos = coord_map[chrom][v.pos]
            translated.append(ParsedVariant(
                chrom=grch38_chrom,
                pos=grch38_pos,
                ref=v.ref,
                alt=v.alt,
                rsid=v.rsid,
                genotype=v.genotype,
                quality=v.quality,
                depth=v.depth,
                genotype_quality=v.genotype_quality,
            ))
        else:
            # Keep variant with "chr" prefix added for potential matching
            chr_chrom = chrom if chrom.startswith("chr") else f"chr{chrom}"
            translated.append(ParsedVariant(
                chrom=chr_chrom,
                pos=v.pos,
                ref=v.ref,
                alt=v.alt,
                rsid=v.rsid,
                genotype=v.genotype,
                quality=v.quality,
                depth=v.depth,
                genotype_quality=v.genotype_quality,
            ))
    return translated


def load_target_positions(allele_defs_dir: Path) -> dict[str, set[int]]:
    """
    Load all target positions from allele definition files.

    Scans all JSON files in directory, extracts variant positions.

    Args:
        allele_defs_dir: Path to directory containing allele definition JSONs.

    Returns:
        Dict mapping chromosome to set of positions.
    """
    target_positions: dict[str, set[int]] = {}

    for json_file in allele_defs_dir.glob("*.json"):
        with open(json_file) as f:
            data = json.load(f)

        for allele in data.get("alleles", []):
            for variant in allele.get("defining_variants", []):
                chrom = variant.get("chrom")
                pos = variant.get("pos")
                if chrom and pos:
                    if chrom not in target_positions:
                        target_positions[chrom] = set()
                    target_positions[chrom].add(pos)

    return target_positions


def filter_by_positions(
    variants: list[ParsedVariant],
    target_positions: dict[str, set[int]],
) -> list[ParsedVariant]:
    """
    Filter variants to only those at target positions.

    Args:
        variants: List of parsed variants.
        target_positions: Dict of chrom -> set of positions.

    Returns:
        Filtered list of variants at target positions.
    """
    filtered = []
    for v in variants:
        if v.chrom in target_positions:
            if v.pos in target_positions[v.chrom]:
                filtered.append(v)
    return filtered


def to_variant_model(
    parsed: ParsedVariant,
    allele_defining: bool = False,
    star_allele: str | None = None,
    functional_impact: str | None = None,
) -> Variant:
    """
    Convert ParsedVariant to elpgx.models.Variant.

    Args:
        parsed: ParsedVariant from VCF parsing.
        allele_defining: Whether this variant defines a star allele.
        star_allele: Star allele name if defining.
        functional_impact: Functional impact description.

    Returns:
        Variant model instance.
    """
    return Variant(
        rsid=parsed.rsid,
        chrom=parsed.chrom,
        pos=parsed.pos,
        ref=parsed.ref,
        alt=parsed.alt,
        genotype=parsed.genotype,
        quality=parsed.quality,
        depth=parsed.depth,
        allele_defining=allele_defining,
        star_allele=star_allele,
        functional_impact=functional_impact,
    )


def variants_at_position(
    variants: list[ParsedVariant],
    chrom: str,
    pos: int,
) -> list[ParsedVariant]:
    """
    Get all variants at a specific position.

    Args:
        variants: List of parsed variants.
        chrom: Chromosome name.
        pos: Position.

    Returns:
        List of variants at that position.
    """
    return [v for v in variants if v.chrom == chrom and v.pos == pos]


def has_variant(
    variants: list[ParsedVariant],
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
) -> ParsedVariant | None:
    """
    Check if a specific variant exists and is present in genotype.

    Args:
        variants: List of parsed variants.
        chrom, pos, ref, alt: Variant coordinates.

    Returns:
        The matching ParsedVariant if found and has alt allele, else None.
    """
    for v in variants:
        if v.chrom == chrom and v.pos == pos and v.ref == ref and v.alt == alt:
            # Check if alt allele is present in genotype
            if "1" in v.genotype or "2" in v.genotype:
                return v
    return None


def call_variants_gatk(
    bam_path: Path,
    reference_path: Path,
    regions_bed: Path,
    output_vcf: Path,
    log_callback: Callable[[str], None] | None = None,
    java_options: str = "-Xmx4g",
) -> Path:
    """
    Call variants from BAM/CRAM using GATK HaplotypeCaller.

    Calls variants only in the regions specified by the BED file (PGx genes).

    Args:
        bam_path: Path to input BAM/CRAM file.
        reference_path: Path to reference FASTA (hg38).
        regions_bed: Path to BED file with target regions.
        output_vcf: Path for output VCF file.
        log_callback: Optional callback for logging messages.
        java_options: Java memory options for GATK.

    Returns:
        Path to the output VCF file.

    Raises:
        RuntimeError: If GATK is not available or variant calling fails.
        FileNotFoundError: If input files are missing.
    """

    def log(msg: str) -> None:
        if log_callback:
            log_callback(msg)

    # Find GATK JAR file (call Java directly to avoid Python wrapper issues)
    gatk_jar = None
    java_cmd = shutil.which("java")

    # Determine base path for bundled resources
    import sys
    if getattr(sys, 'frozen', False):
        # Running as PyInstaller bundle
        bundle_dir = Path(sys._MEIPASS)
    else:
        # Running as script - check project lib directory
        bundle_dir = Path(__file__).parent.parent

    # GATK JAR locations - prioritize bundled JAR
    gatk_jar_paths: list[Path] = [
        # 1. Bundled with application (highest priority)
        bundle_dir / "lib" / "gatk-package-4.6.2.0-local.jar",
        # 2. Conda installations
        Path.home() / "miniconda3" / "share" / "gatk4-4.6.2.0-0" / "gatk-package-4.6.2.0-local.jar",
        Path.home() / "anaconda3" / "share" / "gatk4-4.6.2.0-0" / "gatk-package-4.6.2.0-local.jar",
        # 3. Manual installations
        Path("/opt/gatk/gatk-package-4.6.2.0-local.jar"),
        Path.home() / "gatk" / "gatk-package-4.6.2.0-local.jar",
    ]

    # Add any gatk4 version found in conda share directories
    miniconda_share = Path.home() / "miniconda3" / "share"
    anaconda_share = Path.home() / "anaconda3" / "share"
    if miniconda_share.exists():
        gatk_jar_paths.extend(miniconda_share.glob("gatk4-*/gatk-package-*-local.jar"))
    if anaconda_share.exists():
        gatk_jar_paths.extend(anaconda_share.glob("gatk4-*/gatk-package-*-local.jar"))

    for jar_path in gatk_jar_paths:
        if isinstance(jar_path, Path) and jar_path.exists():
            gatk_jar = str(jar_path)
            log(f"  Using GATK: {jar_path.name}")
            break

    # Find Java executable - prioritize conda installations over system stubs
    # Note: /usr/bin/java on macOS is often a stub that fails without system JRE
    import os
    java_paths = [
        # 1. Conda installations (most reliable)
        Path.home() / "miniconda3" / "bin" / "java",
        Path.home() / "anaconda3" / "bin" / "java",
        Path("/usr/local/opt/openjdk/bin/java"),  # Homebrew OpenJDK
        # 2. Homebrew installations
        Path("/opt/homebrew/bin/java"),
        Path("/usr/local/bin/java"),
        # 3. System Java (check last - macOS stub may exist without JRE)
        Path("/usr/bin/java"),
    ]

    # Also check JAVA_HOME environment variable
    java_home = os.environ.get("JAVA_HOME")
    if java_home:
        java_paths.insert(0, Path(java_home) / "bin" / "java")

    # Find first working Java
    java_cmd = None
    for jp in java_paths:
        if not jp.exists():
            continue
        # Verify this Java actually works
        try:
            result = subprocess.run(
                [str(jp), "-version"],
                capture_output=True,
                timeout=10,
            )
            if result.returncode == 0:
                java_cmd = str(jp)
                break
        except (subprocess.TimeoutExpired, FileNotFoundError, OSError):
            continue

    if java_cmd is None:
        log("Java Runtime not available - using bcftools fallback")
        log("  Note: Install Java (JRE 8+) to use GATK HaplotypeCaller")
        return call_variants_bcftools(
            bam_path=bam_path,
            reference_path=reference_path,
            regions_bed=regions_bed,
            output_vcf=output_vcf,
            log_callback=log_callback,
        )

    if gatk_jar is None:
        log("GATK JAR not found - using bcftools fallback")
        return call_variants_bcftools(
            bam_path=bam_path,
            reference_path=reference_path,
            regions_bed=regions_bed,
            output_vcf=output_vcf,
            log_callback=log_callback,
        )

    # Validate inputs
    if not bam_path.exists():
        raise FileNotFoundError(f"BAM file not found: {bam_path}")
    if not reference_path.exists():
        raise FileNotFoundError(f"Reference not found: {reference_path}")
    if not regions_bed.exists():
        raise FileNotFoundError(f"Regions BED not found: {regions_bed}")

    # Create output directory if needed
    output_vcf.parent.mkdir(parents=True, exist_ok=True)

    # Build GATK HaplotypeCaller command with best practice parameters
    # Reference: GATK Best Practices for Germline Short Variant Discovery
    # Using Java directly to avoid Python wrapper issues in bundled apps
    cmd = [
        java_cmd,
        "-Xmx4g",  # Memory allocation
        "-Dsamjdk.use_async_io_read_samtools=false",
        "-Dsamjdk.use_async_io_write_samtools=true",
        "-Dsamjdk.use_async_io_write_tribble=false",
        "-Dsamjdk.compression_level=2",
        "-jar", gatk_jar,
        "HaplotypeCaller",
        # Input/Output
        "-R", str(reference_path),
        "-I", str(bam_path),
        "-O", str(output_vcf),
        "-L", str(regions_bed),
        # Output mode - emit all confident variants
        "--output-mode", "EMIT_VARIANTS_ONLY",
        # Quality thresholds (GATK best practices)
        "--standard-min-confidence-threshold-for-calling", "30.0",
        "--min-base-quality-score", "20",
        # Annotation for downstream filtering
        "-A", "Coverage",
        "-A", "DepthPerAlleleBySample",
        "-A", "QualByDepth",
        "-A", "FisherStrand",
        "-A", "StrandOddsRatio",
        "-A", "MappingQualityRankSumTest",
        "-A", "ReadPosRankSumTest",
        # Performance optimization
        "--native-pair-hmm-threads", "4",
        # Accuracy improvements
        "--dont-use-soft-clipped-bases", "true",
        # For pharmacogenetic regions - be more sensitive
        "--minimum-mapping-quality", "20",
        # Handle reference/BAM contig mismatches (e.g., HLA contigs in BAM but not in reference)
        "--disable-sequence-dictionary-validation", "true",
    ]

    log(f"Running GATK HaplotypeCaller (Best Practices)...")
    log(f"  Input: {bam_path.name}")
    log(f"  Reference: {reference_path.name}")
    log(f"  Regions: {regions_bed.name}")
    log(f"  Confidence threshold: 30.0")

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=1800,  # 30 minute timeout
        )

        if result.returncode != 0:
            error_msg = result.stderr[-2000:] if result.stderr else "Unknown error"
            raise RuntimeError(f"GATK HaplotypeCaller failed: {error_msg}")

        log(f"  Output: {output_vcf.name}")
        return output_vcf

    except subprocess.TimeoutExpired:
        raise RuntimeError("GATK HaplotypeCaller timed out after 30 minutes")
    except FileNotFoundError:
        log("GATK execution failed - trying bcftools fallback")
        return call_variants_bcftools(
            bam_path=bam_path,
            reference_path=reference_path,
            regions_bed=regions_bed,
            output_vcf=output_vcf,
            log_callback=log_callback,
        )


def call_variants_bcftools(
    bam_path: Path,
    reference_path: Path,
    regions_bed: Path,
    output_vcf: Path,
    log_callback: Callable[[str], None] | None = None,
) -> Path:
    """
    Call variants using bcftools mpileup + call (fallback when GATK unavailable).

    Args:
        bam_path: Path to input BAM/CRAM file.
        reference_path: Path to reference FASTA.
        regions_bed: Path to BED file with target regions.
        output_vcf: Path for output VCF file.
        log_callback: Optional callback for logging messages.

    Returns:
        Path to the output VCF file.

    Raises:
        RuntimeError: If bcftools is not available or fails.
    """

    def log(msg: str) -> None:
        if log_callback:
            log_callback(msg)

    # Check bcftools availability (check common paths for bundled apps)
    bcftools_cmd = shutil.which("bcftools")
    if bcftools_cmd is None:
        # Check common installation paths (Homebrew, conda, system)
        common_paths = [
            "/usr/local/bin/bcftools",      # Homebrew (Intel Mac)
            "/opt/homebrew/bin/bcftools",   # Homebrew (Apple Silicon)
            "/usr/bin/bcftools",            # System
            "/opt/conda/bin/bcftools",      # Conda
            Path.home() / "miniconda3" / "bin" / "bcftools",  # Miniconda
            Path.home() / "anaconda3" / "bin" / "bcftools",   # Anaconda
        ]
        for path in common_paths:
            path = Path(path)
            if path.exists() and path.is_file():
                bcftools_cmd = str(path)
                break

    if bcftools_cmd is None:
        raise RuntimeError(
            "Neither GATK nor bcftools found. Please install one of:\n"
            "  - GATK: https://gatk.broadinstitute.org\n"
            "  - bcftools: https://samtools.github.io/bcftools/"
        )

    # Create output directory
    output_vcf.parent.mkdir(parents=True, exist_ok=True)

    # Determine output format (bgzipped if .gz extension)
    output_format = "z" if str(output_vcf).endswith(".gz") else "v"

    log("Running bcftools mpileup + call (GATK not available)...")
    log(f"  Input: {bam_path.name}")
    log(f"  Note: Install GATK for best practice variant calling")

    # Build bcftools mpileup | bcftools call pipeline
    # Using parameters aligned with GATK best practices where possible
    mpileup_cmd = [
        bcftools_cmd,
        "mpileup",
        "-f", str(reference_path),
        "-R", str(regions_bed),
        # Annotations for quality assessment
        "-a", "FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR,INFO/AD",
        # Quality thresholds (matching GATK)
        "-q", "20",   # Minimum mapping quality
        "-Q", "20",   # Minimum base quality
        # Additional options for accuracy
        "--max-depth", "1000",
        "--min-MQ", "20",
        str(bam_path),
    ]

    call_cmd = [
        bcftools_cmd,
        "call",
        "-m",           # Multiallelic caller
        "-v",           # Output variants only
        "--ploidy", "2",  # Diploid
        f"-O{output_format}",
        "-o", str(output_vcf),
    ]

    try:
        # Run mpileup
        mpileup_proc = subprocess.Popen(
            mpileup_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

        # Pipe to call
        call_proc = subprocess.Popen(
            call_cmd,
            stdin=mpileup_proc.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

        # Wait for completion
        mpileup_proc.stdout.close()  # type: ignore
        stdout, stderr = call_proc.communicate(timeout=1800)

        if call_proc.returncode != 0:
            error_msg = stderr.decode()[:500] if stderr else "Unknown error"
            raise RuntimeError(f"bcftools call failed: {error_msg}")

        log(f"  Output: {output_vcf.name}")

        # Index if bgzipped
        if str(output_vcf).endswith(".gz"):
            subprocess.run(
                [bcftools_cmd, "index", "-t", str(output_vcf)],
                check=True,
            )

        return output_vcf

    except subprocess.TimeoutExpired:
        raise RuntimeError("bcftools timed out after 30 minutes")
