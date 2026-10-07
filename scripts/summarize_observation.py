#!/usr/bin/env python3
"""Create a compact Japanese observation report from NEW postprocessing outputs."""
from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import numpy as np
from sat_string.io import load_case_arrays, load_switch_events, prepare_output_directory, write_json


def summarize(cases, output):
    output = prepare_output_directory(output)
    rows, sensitivity, histories = [], [], {}
    for label, folder in cases:
        m = json.loads((folder / "metrics.json").read_text(encoding="utf-8"))
        p = json.loads((folder / "provenance.json").read_text(encoding="utf-8"))
        f, s, z, a = (m[k] for k in ("formation", "switching", "settling", "altitude"))
        rows.append(dict(case=label, horizon_days=z["observation_end_time_s"] / 86400,
                         max_error_m=f["max_abs_error_m"], rms_error_m=f["rms_error_m"],
                         final_error_m=f["final_max_abs_error_m"],
                         commands=s["n_commands"], starts=s["n_transition_starts"], completions=s["n_transition_completions"],
                         last_switch_s=s["last_switch_time_s"], min_switch_interval_s=s["min_switch_interval_s"],
                         recover_s=z["recovery_time_s"], settle_s=z["settling_time_s"], verification_s=z["verification_horizon_s"],
                         max_lhl_cycles=s["max_lhl_cycles"], secondary=s["secondary_switching_detected"],
                         secondary_started=s["secondary_response_started"],
                         mean_altitude_loss_m=a["mean_altitude_loss_m"], extra_loss_vs_low_m=a["extra_loss_vs_all_low_m"],
                         left_speed_m_s=m["wave"]["left_speed_m_s"], right_speed_m_s=m["wave"]["right_speed_m_s"],
                         source=p["source_directory"]))
        alternate = json.loads((folder / "settling_sensitivity.json").read_text(encoding="utf-8"))
        for item in alternate["evaluations"]:
            sensitivity.append(dict(case=label, **item["settling"]))
        histories[label] = (Path(p["source_directory"]), z["observation_end_time_s"])
    write_json(output / "summary.json", {"rows": rows, "settling_sensitivity": sensitivity})
    for name, table in (("comparison.csv", rows), ("settling_sensitivity.csv", sensitivity)):
        with (output / name).open("w", encoding="utf-8-sig", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(table[0]))
            writer.writeheader()
            writer.writerows(table)
    os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "sat_string_matplotlib"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    fonts = {font.name for font in font_manager.fontManager.ttflist}
    available = next((name for name in ("Yu Gothic", "Meiryo", "MS Gothic", "Noto Sans CJK JP") if name in fonts), None)
    if available is None:
        raise RuntimeError("Japanese font unavailable; install one before creating research figures")
    plt.rcParams.update({"font.family": available, "axes.unicode_minus": False, "font.size": 10})
    # Longest observation for each location; shorter prefixes would draw duplicates.
    chosen = {}
    for row in rows:
        location = row["case"].split("_")[0]
        if location not in chosen or row["horizon_days"] > chosen[location]["horizon_days"]:
            chosen[location] = row
    names = {"central": "中央（11番）", "mid": "中間（6番）", "end": "端部（1番）"}
    # Only load the three longest histories; shorter prefixes need no duplicate arrays.
    for row in chosen.values():
        source, end_s = histories[row["case"]]
        arrays = load_case_arrays(source)
        keep = arrays["time_s"] <= end_s + 1e-9
        histories[row["case"]] = ({key: value[keep] for key, value in arrays.items()}, load_switch_events(source))
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
    for location, row in chosen.items():
        arrays, events = histories[row["case"]]
        t = arrays["time_s"]
        days = t / 86400
        label = names.get(location, location)
        axes[0, 0].plot(days, np.max(np.abs(arrays["gap_error_m"]), axis=1), label=label)
        axes[0, 1].plot(days, np.max(np.abs(arrays["gap_rate_m_s"]), axis=1) * 1000, label=label)
        starts = np.sort([e["start_time"] for e in events if e["start_time"] <= t[-1]])
        axes[1, 0].step(days, np.searchsorted(starts, t, side="right"), where="post", label=label)
        axes[1, 1].plot(days, np.mean(arrays["a_m"][0] - arrays["a_m"], axis=1), label=label)
    axes[0, 0].axhline(50, color="gray", ls="--", label="既存許容値50 m")
    axes[0, 1].axhline(1, color="gray", ls="--", label="既存許容値1 mm/s")
    specs = [("最大絶対間隔誤差", "間隔誤差（m）"), ("最大絶対間隔速度", "間隔速度（mm/s）"),
             ("累積姿勢遷移開始数", "遷移開始数（回）"), ("初期高度からの平均損失", "平均高度損失（m）")]
    for axis, (title, ylabel) in zip(axes.flat, specs, strict=True):
        axis.set_title(title)
        axis.set_ylabel(ylabel)
        axis.set_xlabel("観測時刻（日）")
        axis.grid(alpha=.25)
        axis.legend(fontsize=9)
    fig.suptitle("固定制御条件：外乱位置と観測時間による応答の違い")
    fig.tight_layout()
    fig.savefig(output / "observation_response.png", dpi=180)
    plt.close(fig)
    # Event raster shows completed first and secondary cycles without sampling aliasing.
    fig, axes = plt.subplots(len(chosen), 1, figsize=(12, 3 * len(chosen)), squeeze=False)
    for axis, (location, row) in zip(axes.flat, chosen.items(), strict=True):
        arrays, events = histories[row["case"]]
        for rising, color, marker, label in ((True, "#dc2626", "^", "低抗力→高抗力：開始"), (False, "#2563eb", "v", "高抗力→低抗力：開始")):
            own = [e for e in events if (e["to_mode"] == "H") == rising and e["start_time"] <= arrays["time_s"][-1]]
            axis.scatter([e["start_time"] / 86400 for e in own], [e["satellite_id"] for e in own], color=color, marker=marker, s=18, label=label)
        axis.set_title(names.get(location, location) + "：姿勢遷移イベント")
        axis.set_xlim(0, row["horizon_days"])
        axis.set_ylim(.5, arrays["a_m"].shape[1] + .5)
        axis.set_yticks([1, 6, 11, 16, 21])
        axis.set_xlabel("遷移開始時刻（日）")
        axis.set_ylabel("衛星番号")
        axis.grid(alpha=.25)
        axis.legend(fontsize=9, loc="upper right")
    fig.tight_layout()
    fig.savefig(output / "switching_events.png", dpi=180)
    plt.close(fig)
    for row in rows:
        print({k: row[k] for k in ("case", "max_error_m", "rms_error_m", "starts", "last_switch_s", "settle_s", "secondary")})


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--case", action="append", required=True, help="label=postprocessing-directory")
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    summarize([(item.split("=", 1)[0], Path(item.split("=", 1)[1])) for item in args.case], args.output)


if __name__ == "__main__":
    main()
