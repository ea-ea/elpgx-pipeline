"""
Pydantic data models for the PGx Pipeline.

These models represent the full interpretation chain:
raw variants -> star alleles -> diplotype -> functional status -> phenotype -> prescribing implication -> confidence

All models are designed to serialize directly to the JSON output format.
"""

from datetime import datetime, timezone
from enum import Enum
from typing import Literal

from pydantic import BaseModel, Field


class BiomarkerType(str, Enum):
    """Categories of PGx biomarkers."""
    METABOLIC_ENZYME = "metabolic_enzyme"
    TRANSPORTER = "transporter"
    DRUG_TARGET = "drug_target"
    HLA = "hla"
    OTHER = "other"


class GeneTier(str, Enum):
    """Gene tier classification based on clinical evidence and complexity."""
    GROUP1 = "group1"
    GROUP2 = "group2"
    GROUP3 = "group3"
    GROUP4 = "group4"


class ConfidenceFlag(str, Enum):
    """Quality flags that may affect confidence."""
    LOW_DEPTH = "LOW_DEPTH"
    LOW_QUALITY = "LOW_QUALITY"
    INCOMPLETE_ALLELE = "INCOMPLETE_ALLELE"
    AMBIGUOUS_CALL = "AMBIGUOUS_CALL"
    NO_CALL = "NO_CALL"
    NO_DEFINING_VARIANTS = "NO_DEFINING_VARIANTS"
    COMPLEX_LOCUS_SHORT_READ = "COMPLEX_LOCUS_SHORT_READ"


class Variant(BaseModel):
    """
    A single PGx-relevant variant.

    Represents one observed variant at a pharmacogene locus,
    including its contribution to star-allele assignment.
    """
    rsid: str | None = Field(default=None, description="dbSNP rsID if available")
    chrom: str = Field(..., description="Chromosome (e.g., 'chr10')")
    pos: int = Field(..., description="1-based genomic position")
    ref: str = Field(..., description="Reference allele")
    alt: str = Field(..., description="Alternate allele")
    genotype: str = Field(..., description="Genotype call (e.g., '0/1', '1/1')")
    quality: float = Field(..., description="Variant call quality score")
    depth: int = Field(..., description="Read depth at this position")
    allele_defining: bool = Field(
        default=False,
        description="Whether this variant defines a star allele"
    )
    star_allele: str | None = Field(
        default=None,
        description="Which star allele this variant defines (e.g., '*17')"
    )
    functional_impact: str | None = Field(
        default=None,
        description="Functional consequence (e.g., 'increased function', 'no function')"
    )


class VariantDepthStats(BaseModel):
    """Depth statistics across all variants for a gene."""
    mean_depth: float = Field(..., description="Mean read depth across defining positions")
    min_depth: int = Field(..., description="Minimum read depth at any defining position")


class ConfidenceMetrics(BaseModel):
    """
    Confidence assessment for a gene call.

    Composite score based on depth, quality, and completeness.
    """
    score: float = Field(
        ...,
        ge=0.0,
        le=1.0,
        description="Overall confidence score (0.0-1.0)"
    )
    flags: list[ConfidenceFlag] = Field(
        default_factory=list,
        description="Warning flags affecting this call"
    )
    variant_depth_stats: VariantDepthStats | None = Field(
        default=None,
        description="Depth statistics for variants"
    )


class GeneResult(BaseModel):
    """
    Complete result for a single pharmacogene.

    Contains the full interpretation chain from variants through
    to phenotype and prescribing implications.
    """
    gene: str = Field(..., description="HGNC gene symbol (e.g., 'CYP2C19')")
    gene_tier: GeneTier = Field(..., description="Gene tier classification")
    biomarker_type: BiomarkerType = Field(..., description="Biomarker category")
    diplotype: str = Field(..., description="Star-allele diplotype (e.g., '*1/*17')")
    functional_status: str = Field(
        ...,
        description="Per-allele functional classification"
    )
    phenotype: str = Field(
        ...,
        description="CPIC-defined phenotype (e.g., 'Rapid Metabolizer')"
    )
    activity_score: float | None = Field(
        default=None,
        description="Sum of allele activity values (null if not applicable)"
    )
    prescribing_implication: str = Field(
        ...,
        description="Brief clinical guidance based on CPIC guidelines"
    )
    confidence: ConfidenceMetrics = Field(
        ...,
        description="Confidence metrics for this call"
    )
    variants: list[Variant] = Field(
        default_factory=list,
        description="All PGx-relevant variants detected for this gene"
    )


class ProvenanceInfo(BaseModel):
    """
    Complete provenance information for reproducibility.

    All fields required for full audit trail and result reproduction.
    """
    git_commit: str | None = Field(
        default=None,
        description="Git commit hash of pipeline code"
    )
    git_branch: str | None = Field(
        default=None,
        description="Git branch name"
    )
    nextflow_version: str | None = Field(
        default=None,
        description="Nextflow version used"
    )
    docker_images: dict[str, str] = Field(
        default_factory=dict,
        description="Docker image names with digests (sha256)"
    )
    reference_fasta_md5: str | None = Field(
        default=None,
        description="MD5 checksum of reference FASTA"
    )
    reference_fasta_path: str | None = Field(
        default=None,
        description="Path to reference FASTA used"
    )
    pgx_db_hash: str | None = Field(
        default=None,
        description="SHA256 hash of combined PGx database files"
    )
    input_bam_md5: str | None = Field(
        default=None,
        description="MD5 checksum of input BAM/CRAM"
    )
    variant_caller: str = Field(
        default="GATK HaplotypeCaller",
        description="Variant caller used"
    )
    variant_caller_version: str | None = Field(
        default=None,
        description="Version of variant caller"
    )


class PipelineMeta(BaseModel):
    """Metadata about the pipeline run."""
    sample_id: str = Field(..., description="Sample identifier")
    pipeline_version: str = Field(..., description="Pipeline version string")
    run_date: datetime = Field(
        default_factory=lambda: datetime.now(timezone.utc),
        description="ISO 8601 timestamp of the run"
    )
    reference_build: Literal["GRCh38", "GRCh37"] = Field(
        default="GRCh38",
        description="Reference genome build"
    )
    pgx_db_version: str = Field(..., description="PGx database version (e.g., 'CPIC v1.24')")
    gene_tiers_included: list[GeneTier] = Field(
        default_factory=list,
        description="Which gene tiers were analyzed"
    )
    hla_module_enabled: bool = Field(
        default=False,
        description="Whether HLA calling was enabled"
    )
    provenance: ProvenanceInfo = Field(
        default_factory=ProvenanceInfo,
        description="Complete provenance for reproducibility"
    )


class PgxReport(BaseModel):
    """
    Top-level PGx report output.

    This is the complete JSON output of the pipeline,
    containing metadata, results for all genes, and any warnings.

    Output Structure:
    - schema_version: Version of the output schema format
    - meta: Pipeline run metadata (sample_id, timestamps, versions)
    - results: Per-gene results with:
        - diplotype: Star-allele nomenclature (e.g., "*1/*17")
        - phenotype: Predicted phenotype (e.g., "Rapid Metabolizer")
        - confidence: Quality score (0.0-1.0) with flags
    - warnings: Any pipeline-level warnings
    """
    schema_version: str = Field(
        default="1.0.0",
        description="Schema version for this report format"
    )
    meta: PipelineMeta = Field(..., description="Pipeline run metadata")
    results: list[GeneResult] = Field(
        default_factory=list,
        description="Results for each analyzed gene"
    )
    warnings: list[str] = Field(
        default_factory=list,
        description="Pipeline-level warnings"
    )

    def to_json(self, indent: int = 2) -> str:
        """Serialize report to JSON string."""
        return self.model_dump_json(indent=indent)


# --- Allele Definition Models (for loading reference data) ---

class AlleleVariant(BaseModel):
    """A variant that defines a star allele."""
    rsid: str | None = None
    chrom: str
    pos: int
    ref: str
    alt: str
    functional_impact: str | None = None


class AlleleDefinition(BaseModel):
    """Definition of a single star allele for a gene."""
    allele: str = Field(..., description="Star allele name (e.g., '*17')")
    function: str | None = Field(
        default="Unknown",
        description="Functional category (e.g., 'Increased Function')"
    )
    activity_value: float | None = Field(
        default=None,
        description="Activity score value for this allele"
    )
    defining_variants: list[AlleleVariant] = Field(
        default_factory=list,
        description="Variants that define this allele"
    )
    evidence_strength: str | None = Field(
        default=None,
        description="CPIC evidence strength (e.g., 'Definitive', 'Limited')"
    )
    citations: list[str] | None = Field(
        default=None,
        description="PubMed IDs for supporting evidence"
    )


class GeneAlleleDefinitions(BaseModel):
    """Complete allele definitions for a single gene."""
    gene: str = Field(..., description="HGNC gene symbol")
    chromosome: str = Field(..., description="Chromosome location")
    gene_tier: GeneTier = Field(..., description="Gene tier")
    biomarker_type: BiomarkerType = Field(..., description="Biomarker category")
    reference_allele: str = Field(
        default="*1",
        description="Default reference allele when no variants found"
    )
    alleles: list[AlleleDefinition] = Field(
        default_factory=list,
        description="All defined star alleles for this gene"
    )


# --- Phenotype Mapping Models ---

class DiplotypePhenotype(BaseModel):
    """Mapping from a diplotype to its phenotype."""
    diplotype: str = Field(..., description="Diplotype (e.g., '*1/*17')")
    phenotype: str = Field(..., description="Phenotype label")
    activity_score: float | None = None
    functional_status: str | None = None


class GenePhenotypeMap(BaseModel):
    """Phenotype mappings for a single gene."""
    gene: str
    mappings: list[DiplotypePhenotype] = Field(default_factory=list)


# --- Prescribing Implication Models ---

class PhenotypeImplication(BaseModel):
    """Prescribing implication for a phenotype."""
    phenotype: str
    implication: str


class GeneImplications(BaseModel):
    """Prescribing implications for a single gene."""
    gene: str
    implications: list[PhenotypeImplication] = Field(default_factory=list)
