# Changelog

All notable changes to ElPGx will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.0] - 2026-02-18

### Added

#### GRCh37/hg19 Support
- **Automatic build detection**: Pipeline detects GRCh38 vs GRCh37 from BAM/CRAM header chromosome names and lengths
- **Coordinate translation**: GRCh37 variant positions are automatically translated to GRCh38 coordinates for allele matching
- **GRCh37 BED regions**: Separate target regions file (`pgx_regions_grch37.bed`) for variant calling with GRCh37 references
- **CLI `--build` option**: Manual build override for `elpgx analyze` command

#### Expanded Gene Panel (17 genes)
- **Group 2**: Added CYP2D6 (with phenotype blocking strategy for low-confidence calls), CYP2B6, G6PD
- **Group 3**: Added CACNA1S, RYR1, IFNL3, NAT2
- **Group 4**: Added HLA-A, HLA-B (requires `--enable_hla`)
- **PGx database**: Updated to CPIC v1.24 with official allele definitions for all 17 genes

#### CYP2D6 Strategy Engine
- **Phenotype blocking rules**: Gene-specific strategy files that block phenotype reporting when confidence is below threshold
- **CYP2D6 strategy**: Blocks phenotype assignment for low-confidence calls due to CNV complexity

#### Clinical Validation Framework
- **FDA/CAP compliance**: Validation system for clinical-grade reporting
- **Critical gene completeness checks**: Ensures required alleles present for high-priority genes
- **Allele definition auditing**: Validates CPIC naming conventions in allele definition files

#### Provenance Tracking
- **Full audit trail**: Git commit, Nextflow version, Docker image digests, reference FASTA MD5
- **PGx database hash**: SHA-256 hash of gene tier configuration for reproducibility
- **Schema version**: Report schema updated to v1.1.0 with provenance fields

#### Nextflow Pipeline Improvements
- **GRCh37 support in all stages**: VALIDATE_INPUT detects build, VARIANT_CALL selects correct BED, ALLELE_CALL translates coordinates
- **CRAM alt-contig support**: `--disable-sequence-dictionary-validation` allows CRAMs aligned against GRCh38DH references
- **Docker entrypoint fix**: Replaced deprecated `docker.userEmulation` with `--entrypoint ""` override
- **Container compatibility**: Added `procps` to Docker image for Nextflow process monitoring

#### Docker Image
- **Self-contained image**: GATK 4.5.0.0, samtools, bcftools, Java 21 all included in `elifffarslan/elpgx-python:1.1.0`
- **No separate GATK image needed** for basic CLI usage

#### End-to-End Tests
- **HG002 BAM test** (GRCh37): Validated full pipeline with GIAB HG002 sample
- **HG00100 CRAM test** (GRCh38): Validated CRAM input with 1000 Genomes sample

### Fixed
- **Output directory auto-creation**: Pipeline now creates output directory if it doesn't exist
- **`elpgx analyze` GRCh37 bug**: VCF analysis command was not translating GRCh37 coordinates, causing incorrect diplotype calls (e.g., NAT2 *1/*1 instead of *6/*60)
- **Nextflow `false` rendering**: Fixed Groovy `false` being passed as Python literal instead of string comparison for `hla_module_enabled`
- **Hardcoded reference build**: Report now uses detected build instead of hardcoded "GRCh38"
- **`datetime.utcnow()` deprecation**: Updated to `datetime.now(timezone.utc)`
- **Nextflow `.fai` staging**: Reference index now explicitly staged in VALIDATE_INPUT

### Changed
- **Report schema**: Updated from 1.0.0 to 1.1.0 with provenance and expanded gene support
- **PGx database**: Upgraded from CPIC v1.22 to CPIC v1.24
- **Gene tier assignment**: Reorganized genes across tiers to match CPIC clinical evidence levels

---

## [1.0.1] - 2025-02-10

### Added
- **BED file input**: Target regions BED file support for WES/panel data (CLI: `--regions-bed`, GUI: Target Regions field)
- **Gzipped FASTA support**: Reference genome can now be `.fa.gz` or `.fasta.gz` (with `.fai` and `.gzi` indices)

### Fixed
- **Chromosome naming**: Support both UCSC (`chr1`) and Ensembl/NCBI (`1`) naming conventions
- **macOS file browser**: Fixed file dialog to properly show `.gz` files
- **PyInstaller bundling**: Added pysam and cyvcf2 to hidden imports for standalone builds

### Changed
- Improved validation error messages with clearer guidance
- Updated help text in CLI and GUI to show supported file formats

---

## [1.0.0] - 2025-02-09

### Added

#### Core Pipeline
- **VCF Parser** (`vcf_parser.py`): Robust parsing of VCF/VCF.gz files with support for multi-sample files, phased/unphased genotypes, and quality filtering
- **Allele Caller** (`allele_caller.py`): Star allele calling engine implementing CPIC-compliant algorithms for diplotype determination
- **Phenotype Translator** (`phenotype.py`): Activity score calculation and phenotype assignment following CPIC clinical guidelines
- **Confidence Scoring** (`confidence.py`): Multi-factor confidence assessment considering coverage, phasing, and variant quality
- **Report Generator** (`report.py`): Structured JSON output with complete audit trail and clinical annotations

#### Supported Genes (Initial Release - 8 genes)
- **Group 1**: CYP2C19, CYP2C9, DPYD, NUDT15, SLCO1B1, TPMT, UGT1A1, VKORC1

#### User Interfaces
- **Command-Line Interface** (`cli.py`): Full-featured CLI with rich terminal output, progress indicators, and batch processing support
- **Desktop GUI** (`gui.py`): Cross-platform graphical interface built with tkinter

#### Infrastructure
- **Nextflow Pipeline** (`main.nf`): Scalable workflow for HPC/cloud execution
- **Docker Support**: Containerized deployment with pre-configured environment
- **PyInstaller Build**: Standalone desktop executable generation

---

## Version History

| Version | Date | Highlights |
|---------|------|------------|
| 1.1.0 | 2026-02-18 | GRCh37 support, 17 genes, CPIC v1.24, CRAM support, provenance tracking |
| 1.0.1 | 2025-02-10 | BED file input, gzipped FASTA, chromosome naming fix |
| 1.0.0 | 2025-02-09 | Initial release with 8-gene pipeline |

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
