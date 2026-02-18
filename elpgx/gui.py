"""
ElPGx Desktop GUI

Graphical interface for clinical pharmacogenetic analysis from BAM/CRAM files.
Processes aligned NGS data (hg38) to predict drug response phenotypes.

Pipeline Flow:
    BAM/CRAM → Validate → Variant Call → Allele Call → Phenotype → Report

Usage:
    python -m elpgx.gui
    or
    elpgx-gui
"""

import json
import subprocess
import sys
import threading
from datetime import datetime
from pathlib import Path
from typing import Optional

try:
    import tkinter as tk
    from tkinter import filedialog, messagebox, scrolledtext, ttk
except ImportError:
    print("Error: tkinter not available.")
    print("Please install Python with tkinter support.")
    sys.exit(1)

from elpgx import __version__
from elpgx.clinical_validation import (
    run_clinical_validation,
    ValidationSeverity,
    print_validation_report,
)


def _run_gui_clinical_validation():
    """Run clinical validation on GUI startup."""
    try:
        project_root = Path(__file__).parent.parent
        report = run_clinical_validation(
            project_root,
            strict=False,
            raise_on_failure=False,
        )

        critical = [f for f in report.findings if f.severity == ValidationSeverity.CRITICAL]
        errors = [f for f in report.findings if f.severity == ValidationSeverity.ERROR]

        if critical:
            msg = "CRITICAL VALIDATION ISSUES:\n\n"
            for f in critical:
                msg += f"[{f.code}] {f.message}\n"
            msg += "\nThe pipeline cannot be used for clinical purposes until these are resolved."
            messagebox.showerror("Clinical Validation Failed", msg)
            return False

        if errors:
            msg = "Validation warnings detected:\n\n"
            for f in errors[:5]:  # Show first 5
                msg += f"[{f.code}] {f.message}\n"
            if len(errors) > 5:
                msg += f"\n... and {len(errors) - 5} more issues"
            msg += "\n\nReview these issues for clinical use."
            messagebox.showwarning("Validation Warnings", msg)

        return True
    except Exception as e:
        print(f"Warning: Clinical validation skipped: {e}")
        return True


# Gene tier definitions with full details
GENE_TIERS = {
    "group1": {
        "name": "Group 1 - Core PGx Panel (v1.0)",
        "description": "CPIC Level A genes, low-moderate complexity",
        "genes": [
            ("CYP2C19", "Clopidogrel, PPIs, Voriconazole"),
            ("CYP2C9", "Warfarin, Phenytoin, NSAIDs"),
            ("VKORC1", "Warfarin sensitivity"),
            ("SLCO1B1", "Simvastatin, statins (myopathy risk)"),
            ("TPMT", "Azathioprine, 6-Mercaptopurine"),
            ("DPYD", "5-FU, Capecitabine (safety-critical)"),
            ("UGT1A1", "Irinotecan, Atazanavir"),
            ("NUDT15", "Azathioprine, 6-MP"),
        ],
    },
    "group2": {
        "name": "Group 2 - Complex Genes (v1.1+)",
        "description": "CPIC Level A genes requiring CNV detection",
        "genes": [
            ("CYP2D6", "Codeine, Tamoxifen, SSRIs - requires CNV analysis"),
        ],
    },
    "group3": {
        "name": "Group 3 - Extended Panel (v1.1+)",
        "description": "CPIC Level B/C genes, standard complexity",
        "genes": [
            ("CYP3A5", "Tacrolimus"),
            ("NAT2", "Isoniazid, Hydralazine (acetylator status)"),
            ("IFNL3", "Peginterferon (historical)"),
        ],
    },
    "group4": {
        "name": "Group 4 - HLA Markers (v2.0+)",
        "description": "HLA typing required - immunogenetic markers",
        "genes": [
            ("HLA-B*57:01", "Abacavir hypersensitivity"),
            ("HLA-B*15:02", "Carbamazepine SJS/TEN risk"),
            ("HLA-A*31:01", "Carbamazepine DRESS risk"),
        ],
    },
}


class ElPGxGUI:
    """Main ElPGx Desktop Application."""

    def __init__(self):
        self.root = tk.Tk()
        self.root.title(f"ElPGx - Clinical Pharmacogenetics Analysis v{__version__}")
        self.root.geometry("1000x800")
        self.root.minsize(900, 700)

        # Variables
        self.bam_path = tk.StringVar()
        self.reference_path = tk.StringVar()
        self.bed_path = tk.StringVar()  # Target regions BED file
        self.output_path = tk.StringVar()
        self.sample_id = tk.StringVar()
        self.gene_tiers = tk.StringVar(value="group1")
        self.enable_hla = tk.BooleanVar(value=False)
        self.data_dir = tk.StringVar()

        # Set default output
        self.output_path.set(str(Path.cwd() / "elpgx_report.json"))

        self._create_widgets()
        self._create_menu()

    def _create_menu(self):
        """Create application menu bar."""
        menubar = tk.Menu(self.root)
        self.root.config(menu=menubar)

        # File menu
        file_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="File", menu=file_menu)
        file_menu.add_command(label="Open BAM/CRAM...", command=self._browse_bam)
        file_menu.add_command(label="Set Reference...", command=self._browse_reference)
        file_menu.add_command(label="Set Output...", command=self._browse_output)
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=self.root.quit)

        # Help menu
        help_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Help", menu=help_menu)
        help_menu.add_command(label="Gene Information", command=self._show_gene_info)
        help_menu.add_command(label="About", command=self._show_about)
        help_menu.add_command(label="Documentation", command=self._open_docs)

    def _create_widgets(self):
        """Create main application widgets."""
        # Main container with scrollbar
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.pack(fill=tk.BOTH, expand=True)

        # Header
        header_frame = ttk.Frame(main_frame)
        header_frame.pack(fill=tk.X, pady=(0, 10))

        title_label = ttk.Label(
            header_frame,
            text="ElPGx - Clinical Pharmacogenetics Analysis",
            font=("Helvetica", 16, "bold"),
        )
        title_label.pack(side=tk.LEFT)

        version_label = ttk.Label(
            header_frame,
            text=f"v{__version__}",
            font=("Helvetica", 10),
        )
        version_label.pack(side=tk.RIGHT)

        # Input section
        input_frame = ttk.LabelFrame(main_frame, text="Input Files (Required)", padding="10")
        input_frame.pack(fill=tk.X, pady=5)

        # BAM/CRAM file
        ttk.Label(input_frame, text="BAM/CRAM File:").grid(row=0, column=0, sticky=tk.W, pady=2)
        bam_entry = ttk.Entry(input_frame, textvariable=self.bam_path, width=60)
        bam_entry.grid(row=0, column=1, padx=5, pady=2)
        ttk.Button(input_frame, text="Browse...", command=self._browse_bam).grid(
            row=0, column=2, pady=2
        )
        ttk.Label(
            input_frame,
            text="Aligned reads (hg38), coordinate-sorted with index (.bai/.crai)",
            font=("Helvetica", 9),
            foreground="gray",
        ).grid(row=1, column=1, sticky=tk.W, padx=5)

        # Reference genome
        ttk.Label(input_frame, text="Reference Genome:").grid(row=2, column=0, sticky=tk.W, pady=2)
        ref_entry = ttk.Entry(input_frame, textvariable=self.reference_path, width=60)
        ref_entry.grid(row=2, column=1, padx=5, pady=2)
        ttk.Button(input_frame, text="Browse...", command=self._browse_reference).grid(
            row=2, column=2, pady=2
        )
        ttk.Label(
            input_frame,
            text="hg38/GRCh38 FASTA (.fa, .fasta, .fa.gz, .fasta.gz) with .fai index",
            font=("Helvetica", 9),
            foreground="gray",
        ).grid(row=3, column=1, sticky=tk.W, padx=5)

        # Target regions BED file
        ttk.Label(input_frame, text="Target Regions (BED):").grid(row=4, column=0, sticky=tk.W, pady=2)
        bed_entry = ttk.Entry(input_frame, textvariable=self.bed_path, width=60)
        bed_entry.grid(row=4, column=1, padx=5, pady=2)
        ttk.Button(input_frame, text="Browse...", command=self._browse_bed).grid(
            row=4, column=2, pady=2
        )
        ttk.Label(
            input_frame,
            text="Optional: BED file for WES/panel capture regions (leave empty to use default PGx regions)",
            font=("Helvetica", 9),
            foreground="gray",
        ).grid(row=5, column=1, sticky=tk.W, padx=5)

        # Sample ID
        ttk.Label(input_frame, text="Sample ID:").grid(row=6, column=0, sticky=tk.W, pady=2)
        sample_entry = ttk.Entry(input_frame, textvariable=self.sample_id, width=40)
        sample_entry.grid(row=6, column=1, sticky=tk.W, padx=5, pady=2)
        ttk.Label(
            input_frame,
            text="(optional - auto-detected from BAM @RG header or filename)",
            font=("Helvetica", 9),
            foreground="gray",
        ).grid(row=6, column=2, sticky=tk.W)

        # Output section
        output_frame = ttk.LabelFrame(main_frame, text="Output", padding="10")
        output_frame.pack(fill=tk.X, pady=5)

        ttk.Label(output_frame, text="Report File:").grid(row=0, column=0, sticky=tk.W, pady=2)
        output_entry = ttk.Entry(output_frame, textvariable=self.output_path, width=60)
        output_entry.grid(row=0, column=1, padx=5, pady=2)
        ttk.Button(output_frame, text="Browse...", command=self._browse_output).grid(
            row=0, column=2, pady=2
        )
        ttk.Label(
            output_frame,
            text="Structured JSON report with diplotypes, phenotypes, and prescribing implications",
            font=("Helvetica", 9),
            foreground="gray",
        ).grid(row=1, column=1, sticky=tk.W, padx=5)

        # Gene Tiers section
        tiers_frame = ttk.LabelFrame(
            main_frame, text="Gene Tiers - Select which genes to analyze", padding="10"
        )
        tiers_frame.pack(fill=tk.X, pady=5)

        # Gene tier selection with descriptions
        tier_options = [
            ("group1", "Group 1: Core Panel (8 genes) - CYP2C19, CYP2C9, VKORC1, SLCO1B1, TPMT, DPYD, UGT1A1, NUDT15"),
            ("group1,group2", "Groups 1-2: Core + CYP2D6 (9 genes) - includes CNV-complex gene"),
            ("group1,group2,group3", "Groups 1-3: Extended Panel (12 genes) - adds CYP3A5, NAT2, IFNL3"),
            ("group1,group2,group3,group4", "Groups 1-4: Full Panel + HLA (15 markers) - requires HLA typing"),
        ]

        ttk.Label(tiers_frame, text="Select Gene Tier:").grid(row=0, column=0, sticky=tk.W, pady=2)
        tier_combo = ttk.Combobox(
            tiers_frame,
            textvariable=self.gene_tiers,
            values=[t[0] for t in tier_options],
            width=80,
            state="readonly",
        )
        tier_combo.grid(row=0, column=1, sticky=tk.W, padx=5, pady=2)
        tier_combo.bind("<<ComboboxSelected>>", self._update_tier_description)

        # Tier description label
        self.tier_desc_var = tk.StringVar(value=tier_options[0][1])
        tier_desc_label = ttk.Label(
            tiers_frame,
            textvariable=self.tier_desc_var,
            font=("Helvetica", 9),
            foreground="blue",
            wraplength=700,
        )
        tier_desc_label.grid(row=1, column=1, sticky=tk.W, padx=5, pady=2)

        # Show gene details button
        ttk.Button(tiers_frame, text="View Gene Details...", command=self._show_gene_info).grid(
            row=0, column=2, padx=10, pady=2
        )

        # HLA checkbox
        hla_frame = ttk.Frame(tiers_frame)
        hla_frame.grid(row=2, column=1, sticky=tk.W, pady=5)
        hla_check = ttk.Checkbutton(
            hla_frame,
            text="Enable HLA Typing Module",
            variable=self.enable_hla,
        )
        hla_check.pack(side=tk.LEFT)
        ttk.Label(
            hla_frame,
            text="(requires arcasHLA - for Group 4 markers: HLA-B*57:01, HLA-B*15:02, HLA-A*31:01)",
            font=("Helvetica", 9),
            foreground="gray",
        ).pack(side=tk.LEFT, padx=10)

        # Advanced Options section
        advanced_frame = ttk.LabelFrame(main_frame, text="Advanced Options", padding="10")
        advanced_frame.pack(fill=tk.X, pady=5)

        # Data directory
        ttk.Label(advanced_frame, text="Data Directory:").grid(row=0, column=0, sticky=tk.W, pady=2)
        data_entry = ttk.Entry(advanced_frame, textvariable=self.data_dir, width=50)
        data_entry.grid(row=0, column=1, padx=5, pady=2)
        ttk.Button(advanced_frame, text="Browse...", command=self._browse_data_dir).grid(
            row=0, column=2, pady=2
        )

        # Data directory description
        data_desc = ttk.Label(
            advanced_frame,
            text=(
                "PGx reference data location (leave empty for bundled data).\n"
                "Contains: allele_definitions/ (star allele JSON), phenotype_maps/ (diplotype→phenotype),\n"
                "prescribing_implications/ (clinical guidance), pgx_regions.bed (target regions)"
            ),
            font=("Helvetica", 9),
            foreground="gray",
            justify=tk.LEFT,
        )
        data_desc.grid(row=1, column=1, sticky=tk.W, padx=5, pady=2)

        # Action buttons
        button_frame = ttk.Frame(main_frame)
        button_frame.pack(fill=tk.X, pady=10)

        self.analyze_btn = ttk.Button(
            button_frame,
            text="Run Analysis",
            command=self._run_analysis,
        )
        self.analyze_btn.pack(side=tk.LEFT, padx=5)

        ttk.Button(button_frame, text="View Report", command=self._view_report).pack(
            side=tk.LEFT, padx=5
        )

        ttk.Button(button_frame, text="Clear Log", command=self._clear_log).pack(
            side=tk.RIGHT, padx=5
        )

        # Progress bar
        self.progress = ttk.Progressbar(main_frame, mode="indeterminate")
        self.progress.pack(fill=tk.X, pady=5)

        # Log output
        log_frame = ttk.LabelFrame(main_frame, text="Pipeline Log", padding="5")
        log_frame.pack(fill=tk.BOTH, expand=True, pady=5)

        self.log_text = scrolledtext.ScrolledText(
            log_frame, height=12, font=("Courier", 10), wrap=tk.WORD
        )
        self.log_text.pack(fill=tk.BOTH, expand=True)

        # Status bar
        self.status_var = tk.StringVar(value="Ready - Load BAM/CRAM file and reference genome to begin")
        status_bar = ttk.Label(
            main_frame, textvariable=self.status_var, relief=tk.SUNKEN, anchor=tk.W
        )
        status_bar.pack(fill=tk.X, pady=(5, 0))

    def _update_tier_description(self, event=None):
        """Update tier description when selection changes."""
        tier_descriptions = {
            "group1": "Group 1: Core Panel (8 genes) - CYP2C19, CYP2C9, VKORC1, SLCO1B1, TPMT, DPYD, UGT1A1, NUDT15",
            "group1,group2": "Groups 1-2: Core + CYP2D6 (9 genes) - includes CNV-complex gene",
            "group1,group2,group3": "Groups 1-3: Extended Panel (12 genes) - adds CYP3A5, NAT2, IFNL3",
            "group1,group2,group3,group4": "Groups 1-4: Full Panel + HLA (15 markers) - requires HLA typing",
        }
        self.tier_desc_var.set(tier_descriptions.get(self.gene_tiers.get(), ""))

    def _browse_bam(self):
        """Browse for BAM/CRAM file."""
        path = filedialog.askopenfilename(
            title="Select BAM/CRAM File",
            filetypes=[
                ("Alignment files", "*.bam *.cram"),
                ("BAM files", "*.bam"),
                ("CRAM files", "*.cram"),
                ("All files", "*.*"),
            ],
        )
        if path:
            self.bam_path.set(path)
            # Auto-set output path
            bam_file = Path(path)
            self.output_path.set(str(bam_file.parent / f"{bam_file.stem}_elpgx_report.json"))
            # Try to extract sample ID from filename if not set
            if not self.sample_id.get():
                stem = bam_file.stem
                # Remove common suffixes
                for suffix in [".sorted", ".aligned", ".dedup", ".recal"]:
                    stem = stem.replace(suffix, "")
                self.sample_id.set(stem)

    def _browse_reference(self):
        """Browse for reference genome."""
        path = filedialog.askopenfilename(
            title="Select Reference Genome (hg38)",
            filetypes=[
                ("FASTA files", ("*.fa", "*.fasta")),
                ("Gzipped FASTA", ("*.gz",)),
                ("All files", "*.*"),
            ],
        )
        if path:
            self.reference_path.set(path)

    def _browse_bed(self):
        """Browse for target regions BED file."""
        path = filedialog.askopenfilename(
            title="Select Target Regions BED File",
            filetypes=[
                ("BED files", ("*.bed",)),
                ("Gzipped BED", ("*.gz",)),
                ("All files", "*.*"),
            ],
        )
        if path:
            self.bed_path.set(path)

    def _browse_output(self):
        """Browse for output file location."""
        path = filedialog.asksaveasfilename(
            title="Save Report As",
            defaultextension=".json",
            filetypes=[
                ("JSON files", "*.json"),
                ("All files", "*.*"),
            ],
        )
        if path:
            self.output_path.set(path)

    def _browse_data_dir(self):
        """Browse for data directory."""
        path = filedialog.askdirectory(title="Select PGx Data Directory")
        if path:
            self.data_dir.set(path)

    def _log(self, message: str, level: str = "info"):
        """Add message to log output."""
        timestamp = datetime.now().strftime("%H:%M:%S")
        prefix = {"info": "[INFO]", "success": "[OK]", "error": "[ERROR]", "warning": "[WARN]"}.get(
            level, "[...]"
        )
        self.log_text.insert(tk.END, f"[{timestamp}] {prefix} {message}\n")
        self.log_text.see(tk.END)
        self.root.update_idletasks()

    def _clear_log(self):
        """Clear log output."""
        self.log_text.delete(1.0, tk.END)

    def _run_analysis(self):
        """Run pharmacogenetic analysis pipeline."""
        # Validate inputs
        if not self.bam_path.get():
            messagebox.showerror("Error", "Please select a BAM/CRAM file.")
            return

        if not self.reference_path.get():
            messagebox.showerror(
                "Error",
                "Please select a reference genome (hg38 FASTA).\n\n"
                "The reference is required for variant calling with GATK.",
            )
            return

        if not self.output_path.get():
            messagebox.showerror("Error", "Please specify an output file.")
            return

        bam_file = Path(self.bam_path.get())
        if not bam_file.exists():
            messagebox.showerror("Error", f"BAM/CRAM file not found: {bam_file}")
            return

        ref_file = Path(self.reference_path.get())
        if not ref_file.exists():
            messagebox.showerror("Error", f"Reference genome not found: {ref_file}")
            return

        # Check for index files
        bam_index = bam_file.with_suffix(bam_file.suffix + ".bai")
        if not bam_index.exists():
            bam_index = bam_file.with_suffix(".bai")
        if not bam_index.exists() and bam_file.suffix == ".cram":
            bam_index = bam_file.with_suffix(".crai")

        if not bam_index.exists():
            messagebox.showwarning(
                "Warning",
                f"Index file not found for {bam_file.name}.\n\n"
                "Expected: .bai (for BAM) or .crai (for CRAM)\n"
                "The pipeline may fail without an index.",
            )

        # Disable button and show progress
        self.analyze_btn.config(state=tk.DISABLED)
        self.progress.start()
        self.status_var.set("Running analysis pipeline...")

        # Run in thread to keep GUI responsive
        thread = threading.Thread(target=self._analysis_thread, daemon=True)
        thread.start()

    def _analysis_thread(self):
        """Analysis thread - runs the full pipeline."""
        try:
            self._log("=" * 60)
            self._log("ElPGx Clinical Pharmacogenetics Analysis Pipeline")
            self._log("=" * 60)

            bam_path = Path(self.bam_path.get())
            ref_path = Path(self.reference_path.get())
            output_path = Path(self.output_path.get())
            sample_id = self.sample_id.get() or bam_path.stem

            # Determine data directory
            data_dir_str = self.data_dir.get()
            if data_dir_str:
                data_dir = Path(data_dir_str)
            else:
                # Try bundled data first
                data_dir = Path(__file__).parent.parent / "data"
                if not data_dir.exists():
                    data_dir = Path.cwd() / "data"

            if not data_dir.exists():
                raise FileNotFoundError(
                    f"Data directory not found: {data_dir}\n"
                    "Please specify the PGx data directory in Advanced Options."
                )

            # Determine target regions BED file
            bed_path_str = self.bed_path.get()
            if bed_path_str:
                regions_bed = Path(bed_path_str)
                if not regions_bed.exists():
                    raise FileNotFoundError(f"BED file not found: {regions_bed}")
            else:
                regions_bed = data_dir / "pgx_regions.bed"

            self._log(f"Sample ID: {sample_id}")
            self._log(f"Input BAM: {bam_path}")
            self._log(f"Reference: {ref_path}")
            self._log(f"Target regions: {regions_bed}")
            self._log(f"Data directory: {data_dir}")
            self._log(f"Gene tiers: {self.gene_tiers.get()}")
            self._log("")

            # STAGE 1: Input Validation
            self._log("STAGE 1: Validating input files...")
            from elpgx.validate import validate_bam_file, validate_reference

            validation_result = validate_bam_file(bam_path, ref_path)
            if not validation_result.is_valid:
                raise ValueError(f"Validation failed: {validation_result.errors}")
            self._log("  Input validation passed", "success")
            self._log(f"  Reference build: {validation_result.reference_build}")
            self._log(f"  Sample from BAM: {validation_result.sample_id or 'not found'}")

            # Use sample ID from BAM if available
            if validation_result.sample_id and not self.sample_id.get():
                sample_id = validation_result.sample_id
                self._log(f"  Using sample ID from BAM: {sample_id}")

            # STAGE 2: Variant Calling (PGx regions only)
            self._log("")
            self._log("STAGE 2: Calling variants in PGx regions...")
            self._log("  This may take a few minutes depending on BAM size...")

            # Create temp directory for intermediate files
            output_path.parent.mkdir(parents=True, exist_ok=True)
            temp_dir = output_path.parent / ".elpgx_temp"
            temp_dir.mkdir(exist_ok=True)
            vcf_path = temp_dir / f"{sample_id}_pgx_variants.vcf.gz"

            from elpgx.vcf_parser import call_variants_gatk

            vcf_path = call_variants_gatk(
                bam_path=bam_path,
                reference_path=ref_path,
                regions_bed=regions_bed,
                output_vcf=vcf_path,
                log_callback=lambda msg: self._log(f"  {msg}"),
            )
            self._log("  Variant calling complete", "success")

            # STAGE 3: Parse variants and call alleles
            self._log("")
            self._log("STAGE 3: Parsing variants and calling star alleles...")
            from elpgx.vcf_parser import parse_vcf, load_target_positions, filter_by_positions
            from elpgx.allele_caller import load_gene_definitions, determine_diplotype

            all_variants = parse_vcf(vcf_path)
            target_positions = load_target_positions(data_dir / "allele_definitions")
            variants = filter_by_positions(all_variants, target_positions)
            self._log(f"  Found {len(variants)} PGx-relevant variants", "success")

            # Load gene definitions
            gene_defs = load_gene_definitions(data_dir / "allele_definitions")
            self._log(f"  Loaded {len(gene_defs)} gene definitions")

            # Call diplotypes for selected tiers
            diplotype_calls = []
            tier_list = [t.strip() for t in self.gene_tiers.get().split(",")]

            self._log("")
            self._log("  Star allele calls:")
            for gene_name, gene_def in sorted(gene_defs.items()):
                if gene_def.gene_tier.value in tier_list:
                    gene_variants = [v for v in variants if v.chrom == gene_def.chromosome]
                    call = determine_diplotype(gene_def, gene_variants)
                    diplotype_calls.append(call)
                    confidence_str = f"{call.confidence.score:.0%}"
                    self._log(f"    {gene_name}: {call.diplotype} (confidence: {confidence_str})")

            self._log(f"  Called {len(diplotype_calls)} genes", "success")

            # STAGE 4: Phenotype Translation
            self._log("")
            self._log("STAGE 4: Translating diplotypes to phenotypes...")
            from elpgx.phenotype import translate_diplotype_to_phenotype

            phenotype_results = []
            for call in diplotype_calls:
                result = translate_diplotype_to_phenotype(
                    call,
                    data_dir / "phenotype_maps",
                    data_dir / "prescribing_implications",
                )
                phenotype_results.append(result)

            self._log(f"  Translated {len(phenotype_results)} phenotypes", "success")

            # STAGE 5: Generate Report
            self._log("")
            self._log("STAGE 5: Generating clinical report...")
            from elpgx.report import generate_report

            report = generate_report(
                sample_id=sample_id,
                diplotype_calls=diplotype_calls,
                phenotype_results=phenotype_results,
                output_path=output_path,
                data_dir=data_dir,
                gene_tiers_included=tier_list,
                hla_enabled=self.enable_hla.get(),
                pipeline_version=__version__,
                pgx_db_version="CPIC v1.24",
                reference_build="GRCh38",
            )

            self._log(f"  Report saved: {output_path}", "success")

            # Summary
            self._log("")
            self._log("=" * 60)
            self._log("ANALYSIS COMPLETE", "success")
            self._log("=" * 60)
            self._log("")
            self._log("Results Summary:")
            self._log("-" * 40)

            for result in report.results:
                confidence_pct = f"{result.confidence.score:.0%}"
                self._log(
                    f"  {result.gene}: {result.diplotype} -> {result.phenotype} ({confidence_pct})"
                )

            self.root.after(
                0,
                lambda: messagebox.showinfo(
                    "Analysis Complete",
                    f"PGx analysis completed successfully!\n\n"
                    f"Sample: {sample_id}\n"
                    f"Genes analyzed: {len(report.results)}\n\n"
                    f"Report saved to:\n{output_path}",
                ),
            )

        except Exception as e:
            self._log(f"Pipeline error: {e}", "error")
            import traceback

            self._log(traceback.format_exc(), "error")
            self.root.after(0, lambda: messagebox.showerror("Pipeline Error", str(e)))

        finally:
            self.root.after(0, self._analysis_complete)

    def _analysis_complete(self):
        """Called when analysis completes."""
        self.progress.stop()
        self.analyze_btn.config(state=tk.NORMAL)
        self.status_var.set("Ready")

    def _view_report(self):
        """View generated report."""
        output_path = Path(self.output_path.get())
        if not output_path.exists():
            messagebox.showwarning("Warning", "Report file not found. Run analysis first.")
            return

        # Open in new window
        viewer = tk.Toplevel(self.root)
        viewer.title(f"ElPGx Report - {output_path.name}")
        viewer.geometry("900x700")

        # Create notebook for tabs
        notebook = ttk.Notebook(viewer)
        notebook.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # Load report data
        with open(output_path) as f:
            data = json.load(f)

        # Summary tab
        summary_frame = ttk.Frame(notebook)
        notebook.add(summary_frame, text="Summary")

        summary_text = scrolledtext.ScrolledText(summary_frame, font=("Helvetica", 11), wrap=tk.WORD)
        summary_text.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Build summary
        summary_lines = [
            f"Sample ID: {data.get('meta', {}).get('sample_id', 'Unknown')}",
            f"Analysis Date: {data.get('meta', {}).get('run_date', 'Unknown')}",
            f"Pipeline Version: {data.get('meta', {}).get('pipeline_version', 'Unknown')}",
            f"Reference Build: {data.get('meta', {}).get('reference_build', 'Unknown')}",
            "",
            "=" * 50,
            "PHARMACOGENETIC RESULTS",
            "=" * 50,
            "",
        ]

        for result in data.get("results", []):
            summary_lines.extend(
                [
                    f"Gene: {result.get('gene', 'Unknown')}",
                    f"  Diplotype: {result.get('diplotype', 'Unknown')}",
                    f"  Phenotype: {result.get('phenotype', 'Unknown')}",
                    f"  Confidence: {result.get('confidence', {}).get('score', 0):.0%}",
                    f"  Implication: {result.get('prescribing_implication', 'N/A')[:100]}...",
                    "",
                ]
            )

        summary_text.insert(tk.END, "\n".join(summary_lines))
        summary_text.config(state=tk.DISABLED)

        # Raw JSON tab
        json_frame = ttk.Frame(notebook)
        notebook.add(json_frame, text="Raw JSON")

        json_text = scrolledtext.ScrolledText(json_frame, font=("Courier", 10))
        json_text.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        json_text.insert(tk.END, json.dumps(data, indent=2))
        json_text.config(state=tk.DISABLED)

    def _show_gene_info(self):
        """Show detailed gene tier information."""
        info_window = tk.Toplevel(self.root)
        info_window.title("ElPGx Gene Tier Information")
        info_window.geometry("800x600")

        text = scrolledtext.ScrolledText(info_window, font=("Courier", 10), wrap=tk.WORD)
        text.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        lines = [
            "ELPGX GENE TIER INFORMATION",
            "=" * 60,
            "",
            "Genes are organized into tiers based on clinical evidence",
            "strength and analytical complexity.",
            "",
        ]

        for tier_id, tier_info in GENE_TIERS.items():
            lines.extend(
                [
                    f"{tier_info['name']}",
                    "-" * 50,
                    f"{tier_info['description']}",
                    "",
                    "Genes:",
                ]
            )
            for gene, drugs in tier_info["genes"]:
                lines.append(f"  - {gene}: {drugs}")
            lines.append("")

        lines.extend(
            [
                "=" * 60,
                "CPIC LEVELS:",
                "  Level A: Strong evidence, actionable prescribing guidance",
                "  Level B: Moderate evidence, some prescribing guidance",
                "  Level C: Limited evidence, informational",
                "",
                "For more information: https://cpicpgx.org",
            ]
        )

        text.insert(tk.END, "\n".join(lines))
        text.config(state=tk.DISABLED)

    def _show_about(self):
        """Show about dialog."""
        messagebox.showinfo(
            "About ElPGx",
            f"ElPGx - Clinical Pharmacogenetics Analysis Pipeline\n\n"
            f"Version: {__version__}\n"
            f"Author: Elif Arslan\n\n"
            f"A comprehensive toolkit for pharmacogenetic analysis\n"
            f"from NGS data (BAM/CRAM). Predicts drug response\n"
            f"phenotypes based on genetic variants following\n"
            f"CPIC guidelines.\n\n"
            f"Pipeline Flow:\n"
            f"BAM -> Validate -> Variant Call -> Allele Call -> Phenotype -> Report\n\n"
            f"Supported: 17+ pharmacogenes across 4 tiers",
        )

    def _open_docs(self):
        """Open documentation."""
        import webbrowser

        webbrowser.open("https://github.com/ea-ea/elpgx-pipeline")

    def run(self):
        """Run the application."""
        self.root.mainloop()


def main():
    """Entry point for GUI application."""
    # Run clinical validation before starting GUI
    if not _run_gui_clinical_validation():
        sys.exit(1)

    app = ElPGxGUI()
    app.run()


if __name__ == "__main__":
    main()
