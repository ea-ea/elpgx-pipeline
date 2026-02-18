"""
Provenance collection module for PGx Pipeline.

Collects all information needed for full reproducibility:
- Git commit and branch
- Nextflow version
- Docker image digests
- File checksums
- Tool versions
"""

import hashlib
import os
import subprocess
from pathlib import Path
from typing import Optional

from elpgx.models import ProvenanceInfo


def get_git_commit() -> tuple[Optional[str], Optional[str]]:
    """
    Get current git commit hash and branch.

    Returns:
        Tuple of (commit_hash, branch_name) or (None, None) if not in git repo.
    """
    try:
        commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            stderr=subprocess.DEVNULL,
            text=True,
        ).strip()

        branch = subprocess.check_output(
            ["git", "rev-parse", "--abbrev-ref", "HEAD"],
            stderr=subprocess.DEVNULL,
            text=True,
        ).strip()

        return commit, branch
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None, None


def get_nextflow_version() -> Optional[str]:
    """
    Get Nextflow version if available.

    Returns:
        Version string or None.
    """
    # Check environment variable first (set by Nextflow)
    version = os.environ.get("NXF_VER")
    if version:
        return version

    try:
        output = subprocess.check_output(
            ["nextflow", "-version"],
            stderr=subprocess.DEVNULL,
            text=True,
        )
        # Parse "nextflow version 23.10.0.5889"
        for line in output.split("\n"):
            if "version" in line.lower():
                parts = line.split()
                for i, part in enumerate(parts):
                    if part.lower() == "version" and i + 1 < len(parts):
                        return parts[i + 1]
        return None
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None


def get_docker_image_digest(image_name: str) -> Optional[str]:
    """
    Get Docker image digest (sha256).

    Args:
        image_name: Full image name with tag.

    Returns:
        Digest string or None.
    """
    try:
        output = subprocess.check_output(
            ["docker", "inspect", "--format", "{{index .RepoDigests 0}}", image_name],
            stderr=subprocess.DEVNULL,
            text=True,
        ).strip()
        # Extract digest from "image@sha256:abc123..."
        if "@" in output:
            return output.split("@")[1]
        return None
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None


def calculate_file_md5(file_path: Path) -> Optional[str]:
    """
    Calculate MD5 checksum of a file.

    Args:
        file_path: Path to file.

    Returns:
        MD5 hex string or None if file doesn't exist.
    """
    if not file_path.exists():
        return None

    md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            md5.update(chunk)
    return md5.hexdigest()


def calculate_file_sha256(file_path: Path) -> Optional[str]:
    """
    Calculate SHA256 checksum of a file.

    Args:
        file_path: Path to file.

    Returns:
        SHA256 hex string or None if file doesn't exist.
    """
    if not file_path.exists():
        return None

    sha256 = hashlib.sha256()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            sha256.update(chunk)
    return sha256.hexdigest()


def calculate_directory_hash(dir_path: Path, pattern: str = "*.json") -> Optional[str]:
    """
    Calculate combined SHA256 hash of all files in directory matching pattern.

    Args:
        dir_path: Path to directory.
        pattern: Glob pattern for files to include.

    Returns:
        Combined SHA256 hex string or None.
    """
    if not dir_path.exists():
        return None

    sha256 = hashlib.sha256()

    # Sort files for deterministic ordering
    files = sorted(dir_path.glob(pattern))

    for file_path in files:
        # Include filename in hash for integrity
        sha256.update(file_path.name.encode())
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(8192), b""):
                sha256.update(chunk)

    return sha256.hexdigest() if files else None


def get_gatk_version() -> Optional[str]:
    """
    Get GATK version if available.

    Returns:
        Version string or None.
    """
    try:
        output = subprocess.check_output(
            ["gatk", "--version"],
            stderr=subprocess.DEVNULL,
            text=True,
        )
        # Parse first line which contains version
        for line in output.split("\n"):
            if line.strip():
                # "The Genome Analysis Toolkit (GATK) v4.5.0.0"
                parts = line.split()
                for part in parts:
                    if part.startswith("v") and "." in part:
                        return part
                    if part[0].isdigit() and "." in part:
                        return part
        return None
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None


def collect_provenance(
    reference_path: Optional[Path] = None,
    input_bam_path: Optional[Path] = None,
    data_dir: Optional[Path] = None,
    docker_images: Optional[list[str]] = None,
) -> ProvenanceInfo:
    """
    Collect all provenance information.

    Args:
        reference_path: Path to reference FASTA.
        input_bam_path: Path to input BAM/CRAM.
        data_dir: Path to PGx data directory.
        docker_images: List of Docker image names to get digests for.

    Returns:
        ProvenanceInfo with all available information.
    """
    # Git info
    git_commit, git_branch = get_git_commit()

    # Nextflow version
    nextflow_version = get_nextflow_version()

    # Docker digests
    docker_digests = {}
    if docker_images:
        for image in docker_images:
            digest = get_docker_image_digest(image)
            if digest:
                docker_digests[image] = digest

    # Reference FASTA checksum (MD5 for compatibility with GATK dict)
    reference_md5 = None
    reference_path_str = None
    if reference_path:
        reference_path_str = str(reference_path)
        # Try to read from .md5 file first (common convention)
        md5_file = Path(str(reference_path) + ".md5")
        if md5_file.exists():
            reference_md5 = md5_file.read_text().strip().split()[0]
        else:
            # Calculate if file is small enough (< 1GB)
            if reference_path.exists() and reference_path.stat().st_size < 1_000_000_000:
                reference_md5 = calculate_file_md5(reference_path)

    # Input BAM checksum
    input_bam_md5 = None
    if input_bam_path and input_bam_path.exists():
        # For large BAMs, check for pre-computed md5
        md5_file = Path(str(input_bam_path) + ".md5")
        if md5_file.exists():
            input_bam_md5 = md5_file.read_text().strip().split()[0]
        # Don't compute on the fly for large files

    # PGx database hash
    pgx_db_hash = None
    if data_dir:
        allele_defs_dir = data_dir / "allele_definitions"
        if allele_defs_dir.exists():
            pgx_db_hash = calculate_directory_hash(allele_defs_dir, "*.json")

    # GATK version
    gatk_version = get_gatk_version()

    return ProvenanceInfo(
        git_commit=git_commit,
        git_branch=git_branch,
        nextflow_version=nextflow_version,
        docker_images=docker_digests,
        reference_fasta_md5=reference_md5,
        reference_fasta_path=reference_path_str,
        pgx_db_hash=pgx_db_hash,
        input_bam_md5=input_bam_md5,
        variant_caller="GATK HaplotypeCaller",
        variant_caller_version=gatk_version,
    )


def validate_provenance(provenance: ProvenanceInfo) -> list[str]:
    """
    Validate provenance information completeness.

    Returns:
        List of missing or invalid fields.
    """
    issues = []

    if not provenance.git_commit:
        issues.append("git_commit: Not available (not in git repository?)")

    if not provenance.nextflow_version:
        issues.append("nextflow_version: Not available")

    if not provenance.docker_images:
        issues.append("docker_images: No Docker digests collected")

    if not provenance.reference_fasta_md5:
        issues.append("reference_fasta_md5: Reference checksum not available")

    if not provenance.pgx_db_hash:
        issues.append("pgx_db_hash: PGx database hash not available")

    if not provenance.variant_caller_version:
        issues.append("variant_caller_version: GATK version not available")

    return issues
