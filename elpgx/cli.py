"""
ElPGx Command Line Interface

Provides standalone CLI for pharmacogenetic analysis without Nextflow.
Suitable for single-sample analysis and integration into custom pipelines.

Usage:
    elpgx analyze --vcf <input.vcf> --output <report.json>
    elpgx validate --bam <input.bam>
    elpgx version
"""

import json
import sys
from datetime import datetime
from pathlib import Path
from typing import Optional

try:
    import click
    from rich.console import Console
    from rich.panel import Panel
    from rich.progress import Progress, SpinnerColumn, TextColumn
    from rich.table import Table
except ImportError:
    print("Error: Required packages not installed.")
    print("Please install with: pip install click rich")
    sys.exit(1)

from elpgx import __version__
from elpgx.clinical_validation import (
    run_clinical_validation,
    ClinicalValidationError,
    ValidationSeverity,
)

console = Console()


def _run_clinical_validation(strict: bool = False):
    """
    Run clinical-grade validation checks.

    For production clinical use, set strict=True to enforce all checks.
    """
    try:
        project_root = Path(__file__).parent.parent
        report = run_clinical_validation(
            project_root,
            strict=strict,
            raise_on_failure=False,
        )

        # Show critical/error findings
        critical = [f for f in report.findings if f.severity == ValidationSeverity.CRITICAL]
        errors = [f for f in report.findings if f.severity == ValidationSeverity.ERROR]
        warnings = [f for f in report.findings if f.severity == ValidationSeverity.WARNING]

        if critical:
            console.print("[bold red]CLINICAL VALIDATION FAILED - CRITICAL ISSUES:[/bold red]")
            for f in critical:
                console.print(f"  [{f.code}] {f.message}")
            sys.exit(1)

        if errors and strict:
            console.print("[bold red]CLINICAL VALIDATION FAILED - ERRORS:[/bold red]")
            for f in errors:
                console.print(f"  [{f.code}] {f.message}")
            sys.exit(1)

        if warnings:
            for f in warnings:
                console.print(f"[yellow]Warning [{f.code}]:[/yellow] {f.message}", style="dim")

    except Exception as e:
        console.print(f"[yellow]Warning:[/yellow] Clinical validation skipped: {e}", style="dim")


def print_banner():
    """Print ElPGx banner."""
    banner = """
╔═══════════════════════════════════════════════════════════════════╗
║                     ElPGx v{version}                                  ║
║          Clinical Pharmacogenetics Analysis Pipeline              ║
║                       by Elif Arslan                              ║
╚═══════════════════════════════════════════════════════════════════╝
    """.format(version=__version__)
    console.print(Panel(banner.strip(), style="blue"))


@click.group()
@click.version_option(version=__version__, prog_name="ElPGx")
def main():
    """ElPGx - Clinical Pharmacogenetics Analysis Pipeline

    A comprehensive toolkit for pharmacogenetic analysis from NGS data.
    Predicts drug response phenotypes based on genetic variants.
    """
    # Run clinical-grade validation on CLI initialization
    _run_clinical_validation(strict=False)


@main.command()
@click.option(
    "--bam", "-b",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Input BAM/CRAM file (aligned to hg38 or hg19)"
)
@click.option(
    "--reference", "-r",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Reference genome FASTA (must match BAM alignment build) with .fai index"
)
@click.option(
    "--regions-bed", "-L",
    type=click.Path(exists=True, path_type=Path),
    default=None,
    help="Target regions BED file (for WES/panel data; default: PGx regions only)"
)
@click.option(
    "--output", "-o",
    type=click.Path(path_type=Path),
    default=None,
    help="Output JSON report path (default: <input_name>_pgx_report.json)"
)
@click.option(
    "--sample-id", "-s",
    type=str,
    default=None,
    help="Sample ID (default: derived from BAM @RG header or filename)"
)
@click.option(
    "--data-dir", "-d",
    type=click.Path(exists=True, path_type=Path),
    default=None,
    help="Data directory with allele definitions (default: bundled data)"
)
@click.option(
    "--gene-tiers", "-g",
    type=str,
    default="group1",
    help="Gene tiers to analyze (comma-separated, default: group1)"
)
@click.option(
    "--enable-hla/--no-hla",
    default=False,
    help="Enable HLA typing module (requires arcasHLA)"
)
@click.option(
    "--pretty/--compact",
    default=True,
    help="Pretty-print JSON output"
)
def run(
    bam: Path,
    reference: Path,
    regions_bed: Optional[Path],
    output: Optional[Path],
    sample_id: Optional[str],
    data_dir: Optional[Path],
    gene_tiers: str,
    enable_hla: bool,
    pretty: bool,
):
    """Run full PGx analysis pipeline from BAM/CRAM file.

    This is the primary command for pharmacogenetic analysis.
    It runs the complete pipeline:
      BAM -> Validate -> Variant Call -> Allele Call -> Phenotype -> Report

    For WES/panel data, use --regions-bed to specify the capture kit's target BED file.
    For WGS data, you can omit --regions-bed to use default PGx regions only.

    Gene tiers available:
      - group1: Core panel (CYP2C19, CYP2C9, VKORC1, SLCO1B1, TPMT, DPYD, UGT1A1, NUDT15)
      - group2: Adds CYP2D6 (requires CNV analysis)
      - group3: Adds CYP3A5, NAT2, IFNL3
      - group4: Adds HLA markers (requires --enable-hla)
    """
    print_banner()

    # Derive output filename from input if not specified
    if output is None:
        # Remove all extensions from BAM filename
        bam_stem = bam.name
        for ext in ['.cram', '.bam', '.gz']:
            if bam_stem.endswith(ext):
                bam_stem = bam_stem[:-len(ext)]
        output = bam.parent / f"{bam_stem}_pgx_report.json"

    # Determine data directory
    if data_dir is None:
        data_dir = Path(__file__).parent.parent / "data"
        if not data_dir.exists():
            data_dir = Path.cwd() / "data"

    if not data_dir.exists():
        console.print(f"[red]Error:[/red] Data directory not found: {data_dir}")
        sys.exit(1)

    # Determine target regions BED file (will be overridden for GRCh37)
    if regions_bed is None:
        regions_bed = data_dir / "pgx_regions.bed"

    console.print(f"\n[bold]Input BAM:[/bold] {bam}")
    console.print(f"[bold]Reference:[/bold] {reference}")
    console.print(f"[bold]Target Regions:[/bold] {regions_bed}")
    console.print(f"[bold]Output:[/bold] {output}")
    console.print(f"[bold]Gene Tiers:[/bold] {gene_tiers}")
    console.print()

    detected_build = "GRCh38"  # default

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        console=console,
    ) as progress:
        # Stage 1: Validate input
        task1 = progress.add_task("Validating input files...", total=None)
        try:
            from elpgx.validate import validate_bam_file

            result = validate_bam_file(bam, reference)
            if not result.is_valid:
                progress.update(task1, description=f"[red]✗[/red] Validation failed")
                for err in result.errors:
                    console.print(f"  [red]{err}[/red]")
                sys.exit(1)

            # Use sample ID from BAM if not provided
            if sample_id is None:
                sample_id = result.sample_id or bam.stem

            # Detect genome build
            detected_build = result.reference_build or "GRCh38"

            # For GRCh37, use GRCh37 BED file for variant calling
            if detected_build == "GRCh37":
                grch37_bed = data_dir / "pgx_regions_grch37.bed"
                if grch37_bed.exists():
                    regions_bed = grch37_bed

            build_label = f", build: {detected_build}" if detected_build != "GRCh38" else ""
            progress.update(task1, description=f"[green]✓[/green] Validation passed (sample: {sample_id}{build_label})")
        except Exception as e:
            progress.update(task1, description=f"[red]✗[/red] Validation error: {e}")
            sys.exit(1)

        # Stage 2: Call variants
        task2 = progress.add_task("Calling variants (GATK/bcftools)...", total=None)
        try:
            from elpgx.vcf_parser import call_variants_gatk

            output.parent.mkdir(parents=True, exist_ok=True)
            temp_dir = output.parent / ".elpgx_temp"
            temp_dir.mkdir(exist_ok=True)
            vcf_path = temp_dir / f"{sample_id}_pgx_variants.vcf.gz"

            vcf_path = call_variants_gatk(
                bam_path=bam,
                reference_path=reference,
                regions_bed=regions_bed,
                output_vcf=vcf_path,
            )
            progress.update(task2, description=f"[green]✓[/green] Variant calling complete")
        except Exception as e:
            progress.update(task2, description=f"[red]✗[/red] Variant calling failed: {e}")
            sys.exit(1)

        # Stage 3: Parse VCF and call alleles
        task3 = progress.add_task("Calling star alleles...", total=None)
        try:
            from elpgx.vcf_parser import (
                parse_vcf, load_target_positions, filter_by_positions,
                load_grch37_coordinate_map, translate_variants_grch37_to_grch38,
            )
            from elpgx.allele_caller import load_gene_definitions, determine_diplotype

            all_variants = parse_vcf(vcf_path)

            # If GRCh37, translate coordinates to GRCh38 before allele matching
            if detected_build == "GRCh37":
                coord_map = load_grch37_coordinate_map(data_dir)
                all_variants = translate_variants_grch37_to_grch38(all_variants, coord_map)

            target_positions = load_target_positions(data_dir / "allele_definitions")
            variants = filter_by_positions(all_variants, target_positions)

            gene_defs = load_gene_definitions(data_dir / "allele_definitions")
            diplotype_calls = []
            tier_list = [t.strip() for t in gene_tiers.split(",")]

            for gene_name, gene_def in gene_defs.items():
                if gene_def.gene_tier.value in tier_list:
                    gene_variants = [v for v in variants if v.chrom == gene_def.chromosome]
                    call = determine_diplotype(gene_def, gene_variants)
                    diplotype_calls.append(call)

            progress.update(task3, description=f"[green]✓[/green] Called {len(diplotype_calls)} genes")
        except Exception as e:
            progress.update(task3, description=f"[red]✗[/red] Allele calling failed: {e}")
            sys.exit(1)

        # Stage 4: Translate phenotypes
        task4 = progress.add_task("Translating phenotypes...", total=None)
        try:
            from elpgx.phenotype import translate_diplotype_to_phenotype

            phenotype_results = []
            for call in diplotype_calls:
                result = translate_diplotype_to_phenotype(
                    call,
                    data_dir / "phenotype_maps",
                    data_dir / "prescribing_implications",
                )
                phenotype_results.append(result)

            progress.update(task4, description=f"[green]✓[/green] Translated phenotypes")
        except Exception as e:
            progress.update(task4, description=f"[red]✗[/red] Phenotype translation failed: {e}")
            sys.exit(1)

        # Stage 5: Generate report
        task5 = progress.add_task("Generating report...", total=None)
        try:
            from elpgx.report import generate_report

            report = generate_report(
                sample_id=sample_id,
                diplotype_calls=diplotype_calls,
                phenotype_results=phenotype_results,
                output_path=output,
                data_dir=data_dir,
                gene_tiers_included=tier_list,
                hla_enabled=enable_hla,
                pipeline_version=__version__,
                pgx_db_version="CPIC v1.24",
                reference_build=detected_build,
            )
            progress.update(task5, description=f"[green]✓[/green] Report generated")
        except Exception as e:
            progress.update(task5, description=f"[red]✗[/red] Report generation failed: {e}")
            sys.exit(1)

    # Print summary
    console.print("\n")
    table = Table(title="Pharmacogenetic Results", show_header=True, header_style="bold blue")
    table.add_column("Gene", style="cyan")
    table.add_column("Diplotype", style="green")
    table.add_column("Phenotype", style="yellow")
    table.add_column("Confidence", justify="right")

    for res in report.results:
        conf_str = f"{res.confidence.score:.0%}"
        style = "green" if res.confidence.score >= 0.9 else "yellow" if res.confidence.score >= 0.7 else "red"
        table.add_row(res.gene, res.diplotype, res.phenotype, f"[{style}]{conf_str}[/{style}]")

    console.print(table)
    console.print(f"\n[bold green]Analysis complete![/bold green] Report: {output}")


@main.command()
@click.option(
    "--vcf", "-v",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Input VCF file with genotyped variants (pre-called)"
)
@click.option(
    "--output", "-o",
    type=click.Path(path_type=Path),
    default=None,
    help="Output JSON report path (default: <input_name>_pgx_report.json)"
)
@click.option(
    "--sample-id", "-s",
    type=str,
    default=None,
    help="Sample ID (default: derived from VCF filename)"
)
@click.option(
    "--data-dir", "-d",
    type=click.Path(exists=True, path_type=Path),
    default=None,
    help="Data directory with allele definitions (default: bundled data)"
)
@click.option(
    "--gene-tiers", "-g",
    type=str,
    default="group1",
    help="Gene tiers to analyze (comma-separated, default: group1)"
)
@click.option(
    "--enable-hla/--no-hla",
    default=False,
    help="Enable HLA typing module"
)
@click.option(
    "--pretty/--compact",
    default=True,
    help="Pretty-print JSON output"
)
@click.option(
    "--build",
    type=click.Choice(["GRCh38", "GRCh37"], case_sensitive=False),
    default=None,
    help="Reference build of the VCF (auto-detected from chromosome names if not set)"
)
def analyze(
    vcf: Path,
    output: Optional[Path],
    sample_id: Optional[str],
    data_dir: Optional[Path],
    gene_tiers: str,
    enable_hla: bool,
    pretty: bool,
    build: Optional[str],
):
    """Analyze pre-called VCF file (variant calling already done).

    Use this command if you already have a VCF file with called variants.
    For full pipeline from BAM/CRAM, use 'elpgx run' instead.

    Processes variants from the VCF file, calls star alleles,
    translates to phenotypes, and generates a clinical report.
    """
    print_banner()

    # Derive sample ID from filename if not provided
    if sample_id is None:
        sample_id = vcf.stem.replace(".vcf", "").replace(".gz", "")

    # Derive output filename from input if not specified
    if output is None:
        vcf_stem = vcf.name
        for ext in ['.vcf', '.gz']:
            if vcf_stem.endswith(ext):
                vcf_stem = vcf_stem[:-len(ext)]
        output = vcf.parent / f"{vcf_stem}_pgx_report.json"

    console.print(f"\n[bold]Sample ID:[/bold] {sample_id}")
    console.print(f"[bold]Input VCF:[/bold] {vcf}")
    console.print(f"[bold]Output:[/bold] {output}")
    console.print(f"[bold]Gene Tiers:[/bold] {gene_tiers}")
    console.print()

    # Determine data directory
    if data_dir is None:
        # Use bundled data directory
        data_dir = Path(__file__).parent.parent / "data"
        if not data_dir.exists():
            data_dir = Path.cwd() / "data"

    if not data_dir.exists():
        console.print(f"[red]Error:[/red] Data directory not found: {data_dir}")
        console.print("Please specify --data-dir or ensure data/ exists")
        sys.exit(1)

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        console=console,
    ) as progress:
        # Stage 1: Parse VCF
        task = progress.add_task("Parsing VCF file...", total=None)
        try:
            from elpgx.vcf_parser import parse_vcf, load_target_positions, filter_by_positions

            all_variants = parse_vcf(vcf)

            # Auto-detect or apply build, translate GRCh37 if needed
            detected_build = build or "GRCh38"
            if build is None and all_variants:
                chroms = {v.chrom for v in all_variants}
                if any(not c.startswith("chr") for c in chroms):
                    detected_build = "GRCh37"

            if detected_build == "GRCh37":
                from elpgx.vcf_parser import (
                    load_grch37_coordinate_map, translate_variants_grch37_to_grch38,
                )
                coord_map = load_grch37_coordinate_map(data_dir)
                all_variants = translate_variants_grch37_to_grch38(all_variants, coord_map)

            target_positions = load_target_positions(data_dir / "allele_definitions")
            variants = filter_by_positions(all_variants, target_positions)
            build_label = f" (build: {detected_build})" if detected_build != "GRCh38" else ""
            progress.update(task, description=f"[green]✓[/green] Parsed {len(variants)} pharmacogenetic variants{build_label}")
        except ImportError as e:
            progress.update(task, description=f"[red]✗[/red] VCF parsing failed: {e}")
            console.print("\n[yellow]Note:[/yellow] Install pysam or cyvcf2 for VCF parsing")
            sys.exit(1)
        except Exception as e:
            progress.update(task, description=f"[red]✗[/red] VCF parsing failed: {e}")
            sys.exit(1)

        # Stage 2: Load gene definitions
        task2 = progress.add_task("Loading allele definitions...", total=None)
        try:
            from elpgx.allele_caller import load_gene_definitions

            gene_defs = load_gene_definitions(data_dir / "allele_definitions")
            progress.update(task2, description=f"[green]✓[/green] Loaded {len(gene_defs)} gene definitions")
        except Exception as e:
            progress.update(task2, description=f"[red]✗[/red] Failed to load definitions: {e}")
            sys.exit(1)

        # Stage 3: Call alleles
        task3 = progress.add_task("Calling star alleles...", total=None)
        try:
            from elpgx.allele_caller import determine_diplotype

            diplotype_calls = []
            tier_list = [t.strip() for t in gene_tiers.split(",")]

            for gene_name, gene_def in gene_defs.items():
                if gene_def.gene_tier.value in tier_list:
                    gene_variants = [v for v in variants if v.chrom == gene_def.chromosome]
                    call = determine_diplotype(gene_def, gene_variants)
                    diplotype_calls.append(call)

            progress.update(task3, description=f"[green]✓[/green] Called {len(diplotype_calls)} genes")
        except Exception as e:
            progress.update(task3, description=f"[red]✗[/red] Allele calling failed: {e}")
            sys.exit(1)

        # Stage 4: Translate phenotypes
        task4 = progress.add_task("Translating phenotypes...", total=None)
        try:
            from elpgx.phenotype import translate_diplotype_to_phenotype

            phenotype_results = []
            for call in diplotype_calls:
                result = translate_diplotype_to_phenotype(
                    call,
                    data_dir / "phenotype_maps",
                    data_dir / "prescribing_implications",
                )
                phenotype_results.append(result)

            progress.update(task4, description=f"[green]✓[/green] Translated {len(phenotype_results)} phenotypes")
        except Exception as e:
            progress.update(task4, description=f"[red]✗[/red] Phenotype translation failed: {e}")
            sys.exit(1)

        # Stage 5: Generate report
        task5 = progress.add_task("Generating report...", total=None)
        try:
            from elpgx.report import generate_report

            gene_tiers_file = data_dir / "gene_tiers.json"
            tier_list = [t.strip() for t in gene_tiers.split(",")]

            report = generate_report(
                sample_id=sample_id,
                diplotype_calls=diplotype_calls,
                phenotype_results=phenotype_results,
                output_path=output,
                data_dir=data_dir,
                gene_tiers_included=tier_list,
                hla_enabled=enable_hla,
                pipeline_version=__version__,
                pgx_db_version="CPIC v1.24",
                reference_build=detected_build,
            )

            progress.update(task5, description=f"[green]✓[/green] Report saved to {output}")
        except Exception as e:
            progress.update(task5, description=f"[red]✗[/red] Report generation failed: {e}")
            sys.exit(1)

    # Print summary table
    console.print("\n")
    table = Table(title="Analysis Summary", show_header=True, header_style="bold blue")
    table.add_column("Gene", style="cyan")
    table.add_column("Diplotype", style="green")
    table.add_column("Phenotype", style="yellow")
    table.add_column("Confidence", justify="right")

    for result in report.results:
        confidence_str = f"{result.confidence.score:.0%}"
        if result.confidence.score >= 0.9:
            confidence_style = "green"
        elif result.confidence.score >= 0.7:
            confidence_style = "yellow"
        else:
            confidence_style = "red"

        table.add_row(
            result.gene,
            result.diplotype,
            result.phenotype,
            f"[{confidence_style}]{confidence_str}[/{confidence_style}]",
        )

    console.print(table)
    console.print(f"\n[bold green]Analysis complete![/bold green] Report saved to: {output}")


@main.command()
@click.option(
    "--bam", "-b",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Input BAM/CRAM file to validate"
)
@click.option(
    "--reference", "-r",
    type=click.Path(exists=True, path_type=Path),
    default=None,
    help="Reference FASTA file (.fa/.fa.gz) required for CRAM files"
)
def validate(bam: Path, reference: Optional[Path]):
    """Validate input BAM/CRAM file for pharmacogenetic analysis.

    Checks file format, index presence, reference build (hg38),
    and optionally coverage over pharmacogene regions.
    """
    print_banner()

    console.print(f"\n[bold]Validating:[/bold] {bam}")

    try:
        from elpgx.validate import validate_alignment_file

        result = validate_alignment_file(bam, reference)

        if result.is_valid:
            console.print(f"\n[bold green]✓ Validation passed[/bold green]")
            console.print(f"  Sample ID: {result.sample_id}")
        else:
            console.print(f"\n[bold red]✗ Validation failed[/bold red]")
            for error in result.errors:
                console.print(f"  [red]• {error}[/red]")

        if result.warnings:
            console.print("\n[yellow]Warnings:[/yellow]")
            for warning in result.warnings:
                console.print(f"  [yellow]• {warning}[/yellow]")

        if result.coverage_stats:
            console.print("\n[bold]Coverage Statistics:[/bold]")
            for gene, depth in result.coverage_stats.items():
                if depth >= 30:
                    style = "green"
                elif depth >= 10:
                    style = "yellow"
                else:
                    style = "red"
                console.print(f"  {gene}: [{style}]{depth:.1f}x[/{style}]")

        sys.exit(0 if result.is_valid else 1)

    except ImportError as e:
        console.print(f"\n[red]Error:[/red] {e}")
        console.print("[yellow]Note:[/yellow] Install pysam for BAM validation")
        sys.exit(1)


@main.command()
def version():
    """Display ElPGx version and system information."""
    print_banner()

    console.print(f"\n[bold]Version:[/bold] {__version__}")
    console.print(f"[bold]Python:[/bold] {sys.version}")
    console.print(f"[bold]Platform:[/bold] {sys.platform}")

    # Check optional dependencies
    console.print("\n[bold]Optional Dependencies:[/bold]")

    deps = [
        ("pysam", "BAM/CRAM handling"),
        ("cyvcf2", "VCF parsing (preferred)"),
        ("pydantic", "Data validation"),
        ("jsonschema", "Schema validation"),
    ]

    for pkg, purpose in deps:
        try:
            mod = __import__(pkg)
            ver = getattr(mod, "__version__", "installed")
            console.print(f"  [green]✓[/green] {pkg}: {ver} ({purpose})")
        except ImportError:
            console.print(f"  [red]✗[/red] {pkg}: not installed ({purpose})")


@main.command()
@click.option(
    "--output-dir", "-o",
    type=click.Path(path_type=Path),
    default=Path.cwd(),
    help="Output directory for data files"
)
def init(output_dir: Path):
    """Initialize ElPGx data directory with required files.

    Downloads or copies allele definitions, phenotype maps,
    and prescribing implications to the specified directory.
    """
    print_banner()

    console.print(f"\n[bold]Initializing data directory:[/bold] {output_dir}")

    # Create directories
    data_dir = output_dir / "data"
    dirs_to_create = [
        data_dir / "allele_definitions",
        data_dir / "phenotype_maps",
        data_dir / "prescribing_implications",
    ]

    for d in dirs_to_create:
        d.mkdir(parents=True, exist_ok=True)
        console.print(f"  [green]✓[/green] Created {d}")

    # Copy bundled data if available
    bundled_data = Path(__file__).parent.parent / "data"
    if bundled_data.exists():
        import shutil

        for subdir in ["allele_definitions", "phenotype_maps", "prescribing_implications"]:
            src = bundled_data / subdir
            dst = data_dir / subdir
            if src.exists():
                for f in src.glob("*.json"):
                    shutil.copy(f, dst / f.name)
                console.print(f"  [green]✓[/green] Copied {subdir}")

        # Copy other files
        for f in ["gene_tiers.json", "pgx_regions.bed"]:
            src = bundled_data / f
            if src.exists():
                shutil.copy(src, data_dir / f)
                console.print(f"  [green]✓[/green] Copied {f}")

    console.print(f"\n[bold green]Data directory initialized![/bold green]")
    console.print(f"You can now run: elpgx analyze --vcf <file.vcf> --output report.json")


if __name__ == "__main__":
    main()
