"""
Workflow Optimizer Agent

Identifies bottlenecks or unnecessary steps in the workflow.
Suggests simplifications and validates that behavior remains correct.

Rules:
- No refactors "just because"
- Any change must trace directly to a measurable improvement or user requirement

Deliverables:
- Optimization proposals
- Before/after measurement results
- Tests proving correctness remains unchanged
"""

import json
import sys
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path


@dataclass
class WorkflowStep:
    """Represents a step in the workflow."""

    name: str
    duration_seconds: float
    memory_mb: float
    cpu_percent: float
    input_size_bytes: int
    output_size_bytes: int


@dataclass
class OptimizationProposal:
    """A proposed optimization with justification."""

    step_name: str
    proposal: str
    expected_improvement: str
    risk_level: str  # "low", "medium", "high"
    measurement_required: str
    implemented: bool = False


@dataclass
class WorkflowAnalysis:
    """Complete workflow analysis results."""

    total_duration: float
    bottlenecks: list[str]
    proposals: list[OptimizationProposal]
    baseline_captured: bool


class WorkflowOptimizer:
    """
    Workflow optimization analyzer for PGx Pipeline.

    Analyzes execution traces to identify bottlenecks and propose improvements.
    """

    def __init__(self, trace_file: Path | None = None):
        self.trace_file = trace_file
        self.steps: list[WorkflowStep] = []
        self.proposals: list[OptimizationProposal] = []

    def parse_nextflow_trace(self, trace_path: Path) -> list[WorkflowStep]:
        """
        Parse Nextflow execution trace file.

        Expected format: TSV with columns including name, duration, realtime, etc.
        """
        steps = []

        if not trace_path.exists():
            return steps

        with open(trace_path) as f:
            lines = f.readlines()

        if len(lines) < 2:
            return steps

        header = lines[0].strip().split("\t")
        name_idx = header.index("name") if "name" in header else 0
        duration_idx = header.index("duration") if "duration" in header else -1
        realtime_idx = header.index("realtime") if "realtime" in header else -1
        memory_idx = header.index("peak_rss") if "peak_rss" in header else -1
        cpu_idx = header.index("%cpu") if "%cpu" in header else -1

        for line in lines[1:]:
            cols = line.strip().split("\t")
            if len(cols) <= name_idx:
                continue

            name = cols[name_idx]

            # Parse duration (format: 1m 30s or 90s)
            duration = 0.0
            if duration_idx >= 0 and duration_idx < len(cols):
                duration = self._parse_duration(cols[duration_idx])

            # Parse memory
            memory = 0.0
            if memory_idx >= 0 and memory_idx < len(cols):
                memory = self._parse_memory(cols[memory_idx])

            # Parse CPU
            cpu = 0.0
            if cpu_idx >= 0 and cpu_idx < len(cols):
                try:
                    cpu = float(cols[cpu_idx].replace("%", ""))
                except ValueError:
                    pass

            steps.append(
                WorkflowStep(
                    name=name,
                    duration_seconds=duration,
                    memory_mb=memory,
                    cpu_percent=cpu,
                    input_size_bytes=0,
                    output_size_bytes=0,
                )
            )

        return steps

    def _parse_duration(self, duration_str: str) -> float:
        """Parse duration string to seconds."""
        duration_str = duration_str.strip()
        total_seconds = 0.0

        # Handle formats: "1h 30m 45s", "30m 45s", "45s", "45.5s"
        parts = duration_str.split()
        for part in parts:
            if part.endswith("h"):
                total_seconds += float(part[:-1]) * 3600
            elif part.endswith("m"):
                total_seconds += float(part[:-1]) * 60
            elif part.endswith("s"):
                total_seconds += float(part[:-1])
            else:
                try:
                    total_seconds += float(part)
                except ValueError:
                    pass

        return total_seconds

    def _parse_memory(self, memory_str: str) -> float:
        """Parse memory string to MB."""
        memory_str = memory_str.strip().upper()

        try:
            if memory_str.endswith("GB"):
                return float(memory_str[:-2]) * 1024
            elif memory_str.endswith("MB"):
                return float(memory_str[:-2])
            elif memory_str.endswith("KB"):
                return float(memory_str[:-2]) / 1024
            else:
                return float(memory_str) / (1024 * 1024)  # Assume bytes
        except ValueError:
            return 0.0

    def identify_bottlenecks(
        self,
        threshold_percent: float = 30.0,
    ) -> list[str]:
        """
        Identify steps that consume significant portion of total time.

        Args:
            threshold_percent: Minimum percentage of total time to flag.

        Returns:
            List of bottleneck step names.
        """
        if not self.steps:
            return []

        total_duration = sum(s.duration_seconds for s in self.steps)
        if total_duration == 0:
            return []

        bottlenecks = []
        for step in self.steps:
            percent = (step.duration_seconds / total_duration) * 100
            if percent >= threshold_percent:
                bottlenecks.append(f"{step.name} ({percent:.1f}%)")

        return bottlenecks

    def generate_proposals(self) -> list[OptimizationProposal]:
        """Generate optimization proposals based on analysis."""
        proposals = []

        # Check for variant calling bottleneck
        for step in self.steps:
            if "VARIANT_CALL" in step.name and step.duration_seconds > 300:
                proposals.append(
                    OptimizationProposal(
                        step_name=step.name,
                        proposal="Consider using bcftools mpileup instead of GATK for faster variant calling on targeted regions",
                        expected_improvement="30-50% reduction in variant calling time",
                        risk_level="medium",
                        measurement_required="Compare VCF outputs between GATK and bcftools",
                    )
                )

            # Check for high memory usage
            if step.memory_mb > 8000:
                proposals.append(
                    OptimizationProposal(
                        step_name=step.name,
                        proposal=f"Reduce memory allocation for {step.name} (currently {step.memory_mb:.0f} MB)",
                        expected_improvement="Resource cost reduction",
                        risk_level="low",
                        measurement_required="Monitor peak memory during execution",
                    )
                )

        # General proposals
        if len(self.steps) > 5:
            proposals.append(
                OptimizationProposal(
                    step_name="workflow",
                    proposal="Consider parallelizing allele calling across genes",
                    expected_improvement="Reduced wall-clock time for multi-gene analysis",
                    risk_level="low",
                    measurement_required="Compare single vs parallel execution time",
                )
            )

        return proposals

    def capture_baseline(self, output_file: Path) -> None:
        """Capture current workflow metrics as baseline."""
        baseline = {
            "captured_at": datetime.utcnow().isoformat(),
            "steps": [
                {
                    "name": s.name,
                    "duration_seconds": s.duration_seconds,
                    "memory_mb": s.memory_mb,
                    "cpu_percent": s.cpu_percent,
                }
                for s in self.steps
            ],
            "total_duration": sum(s.duration_seconds for s in self.steps),
        }

        output_file.parent.mkdir(parents=True, exist_ok=True)
        with open(output_file, "w") as f:
            json.dump(baseline, f, indent=2)

    def analyze(self) -> WorkflowAnalysis:
        """Run complete workflow analysis."""
        if self.trace_file and self.trace_file.exists():
            self.steps = self.parse_nextflow_trace(self.trace_file)

        bottlenecks = self.identify_bottlenecks()
        proposals = self.generate_proposals()

        return WorkflowAnalysis(
            total_duration=sum(s.duration_seconds for s in self.steps),
            bottlenecks=bottlenecks,
            proposals=proposals,
            baseline_captured=len(self.steps) > 0,
        )


def main():
    """Entry point for workflow optimizer."""
    import argparse

    parser = argparse.ArgumentParser(description="PGx Pipeline Workflow Optimizer")
    parser.add_argument(
        "--trace",
        type=Path,
        default=Path("results/pipeline_info/execution_trace.txt"),
        help="Path to Nextflow trace file",
    )
    parser.add_argument(
        "--baseline-output",
        type=Path,
        default=Path("tests/baseline_metrics.json"),
        help="Path to save baseline metrics",
    )
    parser.add_argument(
        "--save-baseline",
        action="store_true",
        help="Save current metrics as baseline",
    )
    args = parser.parse_args()

    optimizer = WorkflowOptimizer(trace_file=args.trace)
    analysis = optimizer.analyze()

    print(f"Total Duration: {analysis.total_duration:.1f}s")
    print(f"Baseline Captured: {analysis.baseline_captured}")

    if analysis.bottlenecks:
        print("\nBottlenecks:")
        for bottleneck in analysis.bottlenecks:
            print(f"  - {bottleneck}")

    if analysis.proposals:
        print("\nOptimization Proposals:")
        for i, proposal in enumerate(analysis.proposals, 1):
            print(f"\n  {i}. {proposal.step_name}")
            print(f"     Proposal: {proposal.proposal}")
            print(f"     Expected: {proposal.expected_improvement}")
            print(f"     Risk: {proposal.risk_level}")
            print(f"     Measurement: {proposal.measurement_required}")

    if args.save_baseline:
        optimizer.capture_baseline(args.baseline_output)
        print(f"\nBaseline saved to {args.baseline_output}")


if __name__ == "__main__":
    main()
