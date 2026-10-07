#!/usr/bin/env python3
"""Rerun the frozen three-case baseline once, or check an already rerun suite."""
from pathlib import Path
import argparse
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from sat_string.config import load_config
from sat_string.io import prepare_output_directory, save_case, write_json
from sat_string.metrics import compute_case_metrics, save_metrics
from sat_string.regression import compare_baseline
from sat_string.simulator import STANDARD_CASES, simulate_case


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reference", type=Path, default=Path("tests/fixtures/regression_24h"))
    p.add_argument("--output", type=Path, required=True, help="new output directory")
    p.add_argument("--existing-suite", type=Path, help="check without repeating simulation")
    args = p.parse_args()
    output = prepare_output_directory(args.output)
    suite = args.existing_suite or output / "suite"
    if args.existing_suite is None:
        config = load_config(args.reference / "config.yaml")
        for case in STANDARD_CASES:
            print(f"Running {case.name}", flush=True)
            result = simulate_case(config, case)
            save_case(result, config, suite / case.name)
            save_metrics(compute_case_metrics(result.arrays(), config, result.switch_events), suite / case.name)
    report = compare_baseline(suite, args.reference)
    write_json(output / "regression.json", report)
    print({"passed": report["passed"], "cases": {name: c["passed"] for name, c in report["cases"].items()}}, flush=True)
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
