#!/usr/bin/env python3
"""Run one fixed-controller observation case into a new directory."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys
from time import perf_counter

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from sat_string.config import config_from_dict, load_config
from sat_string.io import prepare_output_directory, save_case
from sat_string.metrics import compute_case_metrics, save_metrics
from sat_string.simulator import STANDARD_CASES, simulate_case


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--config", type=Path, default=Path("tests/fixtures/regression_24h/config.yaml"))
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--hours", type=float, required=True)
    p.add_argument("--case", choices=[c.name for c in STANDARD_CASES], default="controlled_disturbed")
    p.add_argument("--disturbance-satellite", type=int)
    args = p.parse_args()
    data = load_config(args.config).to_dict()
    data["simulation"]["duration_s"] = args.hours * 3600
    if args.disturbance_satellite is not None:
        n = data["formation"]["satellite_count"]
        if not 1 <= args.disturbance_satellite <= n:
            p.error("disturbance satellite must be within the formation")
        data["disturbance"]["center_mode"] = "manual"
        data["disturbance"]["center_s_m"] = (args.disturbance_satellite - (n + 1) / 2) * data["formation"]["target_spacing_m"]
    config = config_from_dict(data)
    case = next(c for c in STANDARD_CASES if c.name == args.case)
    output = prepare_output_directory(args.output)
    start = perf_counter()
    logs = []
    result = simulate_case(config, case, log_callback=logs.append,
                           progress_callback=lambda t, end: print(f"{case.name}: {t/end:.0%} ({perf_counter()-start:.1f} s)", flush=True))
    save_case(result, config, output)
    save_metrics(compute_case_metrics(result.arrays(), config, result.switch_events), output)
    (output / "run.log").write_text("\n".join(logs) + "\n", encoding="utf-8")
    print(f"Saved {output}; elapsed {perf_counter()-start:.1f} s", flush=True)


if __name__ == "__main__":
    main()
