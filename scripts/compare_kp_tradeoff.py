#!/usr/bin/env python3
"""Join existing Kp summary rows and add percentages/marginal costs, no simulation."""
import argparse
import csv
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from sat_string.io import prepare_output_directory, write_json, load_metadata
from sat_string.config import load_config
from sat_string.tradeoff import relative_to_reference, adjacent_cost, dominated_by


def write_csv(path, rows):
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows({k: json.dumps(v, ensure_ascii=False) if isinstance(v, (dict, list)) else v
                         for k, v in row.items()} for row in rows)


def compact(row):
    row = dict(row)
    # These large lists already exist in the source analyses; keep a pointer.
    row["first_response"] = {k: v for k, v in row["first_response"].items() if k != "first_local_cycles"}
    # Same-satellite LOW completion -> next HIGH command, not LH -> LH period.
    with (Path(row["analysis_directory"]) / "local_cycles.csv").open(encoding="utf-8-sig") as f:
        local = list(csv.DictReader(f))
    waits = sorted(float(x["actual_next_wait_s"]) for x in local if x["actual_next_wait_s"])
    row["LOW_completion_to_next_HIGH_s"] = dict(n=len(waits), min=waits[0] if waits else None,
        median=(waits[(len(waits)-1)//2]+waits[len(waits)//2])/2 if waits else None,
        max=waits[-1] if waits else None)
    return row


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--comparison", type=Path, nargs="+", required=True)
    p.add_argument("--days", type=float, nargs="+", default=[30., 60.])
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    rows = []
    for source in args.comparison:
        for row in json.loads(source.read_text(encoding="utf-8")):
            if row["kp_factor"] in [1., 2., 2.5, 3.] and row["days"] in args.days:
                row = compact(row); row["source_comparison"] = str(source.resolve()); rows.append(row)
    ids = [(r["days"], r["kp_factor"]) for r in rows]
    expected = {(day, factor) for day in args.days for factor in [1., 2., 2.5, 3.]}
    if len(set(ids)) != len(ids) or set(ids) != expected:
        raise ValueError("Exactly one saved comparison per requested day/factor required")
    frozen = load_config("tests/fixtures/regression_24h/config.yaml").to_dict()
    adjacent = []; configuration_checks = []
    for row in rows:
        provenance = json.loads((Path(row["analysis_directory"]) / "provenance.json").read_text(encoding="utf-8"))
        case = load_metadata(Path(provenance["source"]))["case"]
        if not case["control_enabled"] or not case["disturbance_enabled"]:
            raise ValueError("Controlled, disturbed case required")
        cfg = load_config(Path(provenance["source"]) / "normalized_config.yaml").to_dict()
        expected_cfg = json.loads(json.dumps(frozen))
        if row["kp_factor"] != 1:
            expected_cfg["controller"]["kp_s2_inv"] = row["kp_factor"] * frozen["derived"]["controller_kp_s2_inv"]
        expected_cfg["derived"]["controller_kp_s2_inv"] = row["kp_factor"] * frozen["derived"]["controller_kp_s2_inv"]
        expected_cfg["simulation"]["duration_s"] = row["days"] * 86400
        if cfg != expected_cfg:
            raise ValueError("Saved case changed more than Kp/duration")
        configuration_checks.append(dict(days=row["days"], kp_factor=row["kp_factor"], passed=True))
    rows.sort(key=lambda r: (r["days"], r["kp_factor"]))
    for day in args.days:
        selected = [r for r in rows if r["days"] == day]; reference = selected[0]
        for row in selected:
            row["relative_to_1"] = relative_to_reference(row, reference)
            row["max_error_dominated_by"] = dominated_by(row, selected, "e_max_m")
            row["rms_dominated_by"] = dominated_by(row, selected, "rms_m")
        adjacent.extend(dict(days=day, **adjacent_cost(a, b)) for a, b in zip(selected, selected[1:]))
    out = prepare_output_directory(args.output)
    write_json(out / "comparison.json", rows); write_csv(out / "comparison.csv", rows)
    write_json(out / "adjacent_cost.json", adjacent); write_csv(out / "adjacent_cost.csv", adjacent)
    write_json(out / "configuration_checks.json", dict(passed=True, cases=configuration_checks))
    lines = ["|日数|Kp倍率|最大誤差 [m]|全期間RMS [m]|終了時誤差 [m]|遷移数|完了往復|最終完了 [日]|追加損失 [m]|",
             "|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for r in rows:
        lines.append(f"|{r['days']:g}|{r['kp_factor']:g}|{r['e_max_m']:.3f}|{r['rms_m']:.3f}|{r['final_max_error_m']:.3f}|{r['transitions']}|{r['completed_cycles']}|{r['last_switch_day']:.5f}|{r['extra_altitude_loss_m']:.3f}|")
    (out / "table.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines), flush=True)


if __name__ == "__main__":
    main()
