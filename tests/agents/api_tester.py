"""
API Tester Agent

Validates all API endpoints and contracts (status codes, schemas, error cases).
Tests backward compatibility expectations where applicable.

Deliverables:
- API contract tests
- Negative tests (invalid inputs, missing auth, etc.)
- Mock server or contract schema enforcement
"""

import json
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import jsonschema


@dataclass
class ContractTestResult:
    """Result of a contract validation test."""

    contract_name: str
    passed: bool
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)


class APITester:
    """
    API and contract tester for PGx Pipeline.

    Validates intermediate and final outputs against their schemas/contracts.
    """

    def __init__(self, schemas_dir: Path):
        self.schemas_dir = schemas_dir
        self.results: list[ContractTestResult] = []

    def load_schema(self, schema_name: str) -> dict | None:
        """Load a JSON schema file."""
        schema_path = self.schemas_dir / schema_name
        if not schema_path.exists():
            return None
        with open(schema_path) as f:
            return json.load(f)

    def validate_against_schema(
        self,
        data: dict,
        schema: dict,
    ) -> list[str]:
        """Validate data against JSON schema, return list of errors."""
        validator = jsonschema.Draft202012Validator(schema)
        errors = []
        for error in validator.iter_errors(data):
            path = ".".join(str(p) for p in error.path) or "root"
            errors.append(f"{path}: {error.message}")
        return errors

    def test_report_schema(self, report_path: Path) -> ContractTestResult:
        """Test final report against schema."""
        schema = self.load_schema("pgx_report_schema.json")
        if schema is None:
            return ContractTestResult(
                contract_name="pgx_report_schema",
                passed=False,
                errors=["Schema file not found"],
            )

        if not report_path.exists():
            return ContractTestResult(
                contract_name="pgx_report_schema",
                passed=False,
                errors=[f"Report file not found: {report_path}"],
            )

        with open(report_path) as f:
            report = json.load(f)

        errors = self.validate_against_schema(report, schema)

        return ContractTestResult(
            contract_name="pgx_report_schema",
            passed=len(errors) == 0,
            errors=errors,
        )

    def test_alleles_contract(self, alleles_path: Path) -> ContractTestResult:
        """Test alleles.json against contract."""
        required_fields = ["sample_id", "reference", "calls"]

        if not alleles_path.exists():
            return ContractTestResult(
                contract_name="alleles_contract",
                passed=False,
                errors=[f"Alleles file not found: {alleles_path}"],
            )

        with open(alleles_path) as f:
            data = json.load(f)

        errors = []
        for field in required_fields:
            if field not in data:
                errors.append(f"Missing required field: {field}")

        # Validate each call entry
        for i, call in enumerate(data.get("calls", [])):
            if "gene" not in call:
                errors.append(f"calls[{i}]: missing 'gene' field")
            if "diplotype" not in call:
                errors.append(f"calls[{i}]: missing 'diplotype' field")
            if "call_quality" not in call:
                errors.append(f"calls[{i}]: missing 'call_quality' field")
            elif "score" not in call.get("call_quality", {}):
                errors.append(f"calls[{i}].call_quality: missing 'score' field")

        return ContractTestResult(
            contract_name="alleles_contract",
            passed=len(errors) == 0,
            errors=errors,
        )

    def test_phenotypes_contract(self, phenotypes_path: Path) -> ContractTestResult:
        """Test phenotypes.json against contract."""
        required_fields = ["sample_id", "phenotypes"]

        if not phenotypes_path.exists():
            return ContractTestResult(
                contract_name="phenotypes_contract",
                passed=False,
                errors=[f"Phenotypes file not found: {phenotypes_path}"],
            )

        with open(phenotypes_path) as f:
            data = json.load(f)

        errors = []
        for field in required_fields:
            if field not in data:
                errors.append(f"Missing required field: {field}")

        # Validate each phenotype entry
        for i, pheno in enumerate(data.get("phenotypes", [])):
            if "gene" not in pheno:
                errors.append(f"phenotypes[{i}]: missing 'gene' field")
            if "diplotype" not in pheno:
                errors.append(f"phenotypes[{i}]: missing 'diplotype' field")
            if "phenotype" not in pheno:
                errors.append(f"phenotypes[{i}]: missing 'phenotype' field")

        return ContractTestResult(
            contract_name="phenotypes_contract",
            passed=len(errors) == 0,
            errors=errors,
        )

    def test_negative_inputs(self) -> list[ContractTestResult]:
        """Test handling of invalid inputs."""
        results = []

        # Test: Empty report should fail schema validation
        empty_report: dict[str, Any] = {}
        schema = self.load_schema("pgx_report_schema.json")
        if schema:
            errors = self.validate_against_schema(empty_report, schema)
            results.append(
                ContractTestResult(
                    contract_name="negative_empty_report",
                    passed=len(errors) > 0,  # Should have errors
                    errors=[] if len(errors) > 0 else ["Empty report should fail validation"],
                )
            )

        # Test: Missing sample_id should fail
        incomplete_report = {"meta": {}, "results": []}
        if schema:
            errors = self.validate_against_schema(incomplete_report, schema)
            results.append(
                ContractTestResult(
                    contract_name="negative_missing_sample_id",
                    passed=len(errors) > 0,
                    errors=[] if len(errors) > 0 else ["Missing sample_id should fail"],
                )
            )

        return results

    def run_all_tests(
        self,
        output_dir: Path,
    ) -> tuple[int, int, list[ContractTestResult]]:
        """
        Run all contract tests.

        Returns:
            Tuple of (passed_count, failed_count, results).
        """
        all_results = []

        # Test report schema
        report_files = list(output_dir.glob("*_pgx_report.json"))
        for report_file in report_files:
            all_results.append(self.test_report_schema(report_file))

        # Test intermediate files
        alleles_files = list(output_dir.glob("**/alleles.json"))
        for alleles_file in alleles_files:
            all_results.append(self.test_alleles_contract(alleles_file))

        phenotypes_files = list(output_dir.glob("**/phenotypes.json"))
        for phenotypes_file in phenotypes_files:
            all_results.append(self.test_phenotypes_contract(phenotypes_file))

        # Negative tests
        all_results.extend(self.test_negative_inputs())

        passed = sum(1 for r in all_results if r.passed)
        failed = len(all_results) - passed

        return passed, failed, all_results


def main():
    """Entry point for API tester."""
    import argparse

    parser = argparse.ArgumentParser(description="PGx Pipeline API Tester")
    parser.add_argument(
        "--schemas",
        type=Path,
        default=Path("schemas"),
        help="Path to schemas directory",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results"),
        help="Path to pipeline output directory",
    )
    args = parser.parse_args()

    tester = APITester(schemas_dir=args.schemas)
    passed, failed, results = tester.run_all_tests(args.output)

    print(f"API Contract Tests: {passed}/{passed + failed} passed")

    for result in results:
        status = "PASS" if result.passed else "FAIL"
        print(f"  [{status}] {result.contract_name}")
        for error in result.errors:
            print(f"         ERROR: {error}")

    sys.exit(0 if failed == 0 else 1)


if __name__ == "__main__":
    main()
