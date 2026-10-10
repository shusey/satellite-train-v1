#!/usr/bin/env python3
"""Three all-LOW single-pass experiments; reruns require a fresh output root."""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import copy
import csv
import os

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from sat_string.config import load_config, config_from_dict, dump_normalized_config
from sat_string.io import prepare_output_directory, save_case, write_json, load_case_arrays
from sat_string.relative_motion import simulate_all_low_exposure, evaluate_motion, local_background

NAMES = ("earth_fixed", "tad_730", "comoving")
LABELS = ("地球固定", "TAD伝播（730 m/s）", "初期速度で共移動")


def prepare(output):
    prepare_output_directory(output)
    reference = ROOT / "tests/fixtures/regression_24h/config.yaml"
    source = ROOT / "docs/observed_density_design/data/design_candidates.json"
    event = next(x for x in json.loads(source.read_text(encoding="utf-8"))
                 if x["case"] == "largest_fitted_D100")
    assert event["accepted"] and event["r_squared"] > 0.8
    base = load_config(reference)
    vsat = base.derived.reference_radius_m * base.derived.reference_mean_motion_rad_s
    alignment = base.disturbance.start_time_s  # existing 10800 s, phase only
    baseline = base.to_dict()
    baseline.pop("derived", None)
    shared = copy.deepcopy(baseline)
    shared["simulation"]["duration_s"] = 2 * 86400
    shared["disturbance"].update(type="earth_track_gaussian", center_mode="manual",
        amplitude_fraction=event["amplitude_fraction_local_bg"],
        fwhm_m=event["fwhm_km"]*1000)
    configurations = {}
    for name, speed in zip(NAMES, (0.0, 730.0, vsat)):
        data = copy.deepcopy(shared)
        data["disturbance"]["propagation_speed_m_s"] = speed
        data["disturbance"]["center_s_m"] = (vsat-speed)*alignment
        cfg = config_from_dict(data)
        (output / "configs").mkdir(exist_ok=True)
        dump_normalized_config(cfg, output / "configs" / (name+".yaml"))
        configurations[name] = cfg.to_dict()
    # Mechanically verify physical and numerical settings, and the allowed differences.
    def flatten(value, prefix=""):
        result = {}
        for key, val in value.items():
            path = prefix+key
            if isinstance(val, dict):
                result.update(flatten(val, path+"."))
            else:
                result[path] = val
        return result
    reference_flat = flatten(configurations[NAMES[0]])
    differences = {name: [key for key, val in flatten(configurations[name]).items()
                          if val != reference_flat[key]] for name in NAMES}
    assert all(set(keys) <= {"disturbance.propagation_speed_m_s", "disturbance.center_s_m"}
               for keys in differences.values())
    fixed_sections = ("constants", "spacecraft", "formation", "controller", "atmosphere",
                      "recovery", "visualization")
    assert all(configurations[n][s] == base.to_dict()[s]
               for n in NAMES for s in fixed_sections)
    assert all(configurations[n]["simulation"][s] == base.to_dict()["simulation"][s]
               for n in NAMES for s in ("integration_step_s", "output_step_s", "integrator"))
    write_json(output / "study.json", dict(
        source=str(source.relative_to(ROOT)), source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        selected_event=event, alignment_time_s=alignment, configurations=configurations,
        input_differences_vs_earth_fixed=differences,
        fixed_sections_verified=list(fixed_sections),
        tad_reference="https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2007GL030243",
        tad_event="2003-05-29, CHAMP ~400 km, northern dayside auroral-to-equatorial TAD, 730 m/s",
        coordinates="Nonperiodic unwrapped x=a0*lambda; single-pass Earth-like 1-D proxy, no Earth rotation",
        center_phase_rule="center_s_m=(v_sat_initial-v_dist)*10800; same nominal central encounter epoch",
        time_envelope="none; legacy start/duration fields do not gate earth_track_gaussian",
        comoving="constant initial satellite speed, independent of actual evolving train; ongoing forcing",
        background="unchanged software baseline exponential; observed event background NOT substituted",
    ))
    # Capture existing saved results and fixtures before simulations. Hash all files in
    # these protected roots, excluding only this fresh output directory.
    protected = {}
    for directory in ("outputs", "runs", "tests/fixtures", "docs"):
        for path in (ROOT / directory).rglob("*"):
            if path.is_file() and output.resolve() not in path.resolve().parents:
                protected[str(path.relative_to(ROOT)).replace("\\", "/")] = hashlib.sha256(path.read_bytes()).hexdigest()
    write_json(output / "preserved_sha256.json", protected)
    print("Prepared; protected files:", len(protected), flush=True)


def run(output, name):
    destination = prepare_output_directory(output / "runs" / name)
    config = load_config(output / "configs" / (name+".yaml"))
    result, exposure = simulate_all_low_exposure(config,
        lambda t, end: print(f"{name}: {t/end:.0%}", flush=True))
    save_case(result, config, destination)
    np.savez_compressed(destination / "exposure.npz", **exposure)
    write_json(destination / "motion_metrics.json", evaluate_motion(result.arrays(), exposure, config))


def analyze(output, report):
    prepare_output_directory(report)
    (report / "figures").mkdir()
    (report / "data").mkdir()
    os.environ.setdefault("MPLCONFIGDIR", str(output.resolve() / "matplotlib_cache"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams["font.family"] = "Meiryo"
    plt.rcParams["axes.unicode_minus"] = False
    data, summary = [], {}
    for name in NAMES:
        cfg = load_config(output / "configs" / (name+".yaml"))
        arrays = load_case_arrays(output / "runs" / name)
        with np.load(output / "runs" / name / "exposure.npz") as f:
            exposure = {key: f[key] for key in f.files}
        summary[name] = evaluate_motion(arrays, exposure, cfg)
        data.append((cfg, arrays, exposure))
        write_json(report / "data" / (name+"_metrics.json"), summary[name])
        # Compact per-satellite integrals and signed gap states.
        with (report / "data" / (name+"_satellite_exposure.csv")).open("w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(["satellite_id", "J_rho_kg_s_m3", "J_delta_kg_s_m3", "J_drag_m_s", "J_drag_delta_m_s"])
            writer.writerows(zip(range(1,22), *(exposure[key][-1] for key in ("J_rho","J_delta","J_drag","J_drag_delta"))))
        # Small time-series CSV supports figures without needing large local NPZ.
        stride = max(1, round(600/cfg.simulation.output_step_s))
        with (report / "data" / (name+"_history_600s.csv")).open("w", newline="", encoding="utf-8") as f:
            writer=csv.writer(f)
            writer.writerow(["time_s","max_error_m","rms_error_m","max_gap_rate_m_s",
                             "J_delta_spread_kg_s_m3","max_instant_density_difference_kg_m3"])
            writer.writerows(zip(arrays["time_s"][::stride],np.max(np.abs(arrays["gap_error_m"]),axis=1)[::stride],
                np.sqrt(np.mean(arrays["gap_error_m"]**2,axis=1))[::stride],
                np.max(np.abs(arrays["gap_rate_m_s"]),axis=1)[::stride],
                np.ptp(exposure["J_delta"],axis=1)[::stride],np.ptp(arrays["density_kg_m3"],axis=1)[::stride]))
        with (report/"data"/(name+"_gaps_600s.csv")).open("w",newline="",encoding="utf-8") as f:
            writer=csv.writer(f)
            writer.writerow(["time_s"]+[f"gap_{i}_error_m" for i in range(1,21)]
                            +[f"gap_{i}_rate_m_s" for i in range(1,21)])
            writer.writerows(np.column_stack((arrays["time_s"][::stride],arrays["gap_error_m"][::stride],
                                               arrays["gap_rate_m_s"][::stride])))
        rho_delta=arrays["density_kg_m3"]-local_background(arrays["a_m"],cfg)
        summary[name]["exposure"]["max_instant_increment_difference_kg_m3"]=float(np.max(np.ptp(rho_delta,axis=1)))

    def finish(fig, filename):
        fig.tight_layout()
        fig.savefig(report/"figures"/filename,dpi=160)
        plt.close(fig)

    fig, axes=plt.subplots(1,3,figsize=(14,4))
    for ax,(cfg,a,j),label in zip(axes,data,LABELS):
        t=a["time_s"]
        center=cfg.disturbance.center_s_m+cfg.disturbance.propagation_speed_m_s*t
        offset=(cfg.derived.reference_radius_m*a["lambda_rad"]-center[:,None])/1000
        mask=((t>=10300)&(t<=11300)) if label!=LABELS[2] else np.ones(t.shape,dtype=bool)
        for sat in (0,10,20):
            ax.plot(t[mask]/3600,offset[mask,sat],label=f"衛星{sat+1}")
        ax.axhspan(-cfg.disturbance.fwhm_m/2000,cfg.disturbance.fwhm_m/2000,color="grey",alpha=.2,label="半値幅")
        ax.set(title=label,xlabel="時刻 [時間]",ylabel="外乱中心からの距離 [km]")
        ax.legend()
    finish(fig,"01_relative_motion.png")

    fig, axes=plt.subplots(1,3,figsize=(14,4))
    for ax,(cfg,a,j),label in zip(axes,data,LABELS):
        t=a["time_s"]
        mask=(t>10400)&(t<11200) if label!=LABELS[2] else np.ones(t.shape,dtype=bool)
        for sat in (0,10,20):
            ax.plot(t[mask]/3600,a["density_kg_m3"][mask,sat]*1e12,label=f"衛星{sat+1}")
        ax.set(title=label,xlabel="時刻 [時間]",ylabel="密度 [10⁻¹² kg/m³]")
        ax.legend()
    finish(fig,"02_density.png")

    fig,axes=plt.subplots(1,2,figsize=(11,4))
    for (cfg,a,j),label in zip(data,LABELS):
        axes[0].plot(np.arange(1,22),j["J_delta"][-1]*1e12,"o-",label=label)
        if label != LABELS[2]:
            axes[1].plot(np.arange(1,22),j["J_delta"][-1]*1e12,"o-",label=label)
    axes[0].set(title="各衛星の積算外乱密度（2日）",xlabel="衛星番号",ylabel="積算密度 [10⁻¹² kg·s/m³]")
    pass_values=np.concatenate([j["J_delta"][-1]*1e12 for _,_,j in data[:2]])
    axes[1].set(title="通過ケースの拡大",xlabel="衛星番号",ylabel="積算密度 [10⁻¹² kg·s/m³]",
                ylim=(pass_values.min()*.98,pass_values.max()*1.02))
    for ax in axes: ax.legend()
    finish(fig,"03_exposure.png")

    fig,axes=plt.subplots(1,2,figsize=(11,4))
    for name,(cfg,a,j),label in zip(NAMES,data,LABELS):
        snapshots=summary[name]["post_passage"]
        rate=np.asarray(snapshots["one_day_after"]["gap_rate_m_s"]) if snapshots else a["gap_rate_m_s"][-1]
        axes[0].plot(np.arange(1,21),rate,"o-",label=label+("（通過1日後）" if snapshots else "（2日・作用中）"))
        if snapshots:
            axes[1].plot(np.arange(1,21),rate,"o-",label=label)
    for ax in axes: ax.set(xlabel="隣接間隔の番号",ylabel="相対速度 [m/s]");ax.legend()
    axes[0].set_title("符号付き速度：残留と作用中を区別")
    axes[1].set_title("通過ケースの残留速度を拡大")
    axes[1].set_ylim(-1e-9,1e-9)
    axes[1].text(.5,.75,"細かな符号・大小は数値精度の範囲",transform=axes[1].transAxes,ha="center")
    finish(fig,"04_gap_rates.png")

    fig,axes=plt.subplots(2,1,figsize=(10,7),sharex=True)
    for (cfg,a,j),label in zip(data,LABELS):
        t=a["time_s"]/86400
        axes[0].plot(t,np.maximum(np.max(np.abs(a["gap_error_m"]),axis=1),1e-8),label=label)
        axes[1].plot(t,np.maximum(np.sqrt(np.mean(a["gap_error_m"]**2,axis=1)),1e-8),label=label)
    for ax in axes: ax.set_yscale("log");ax.legend();ax.grid(alpha=.2);ax.set_ylabel("間隔誤差 [m]")
    axes[0].set_title("最大絶対間隔誤差（縦軸は対数）")
    axes[1].set(title="瞬間RMS間隔誤差",xlabel="時刻 [日]")
    finish(fig,"05_errors.png")

    fig,axes=plt.subplots(1,3,figsize=(14,4))
    for ax,(cfg,a,j),label in zip(axes,data,LABELS):
        prediction=-3*np.diff(j["J_drag_delta"][-1])
        actual=a["gap_rate_m_s"][-1]
        ax.scatter(prediction,actual)
        lo,hi=min(prediction.min(),actual.min()),max(prediction.max(),actual.max())
        ax.plot([lo,hi],[lo,hi],"--",color="grey",label="一致線")
        ax.set(title=label,xlabel="−3 × 隣接積算抗力増分差 [m/s]",ylabel="2日時点の隣接相対速度 [m/s]")
        ax.legend()
        if label != LABELS[2]:
            ax.set_xlim(-1e-9,1e-9)
            ax.set_ylim(-1e-9,1e-9)
            ax.text(.5,.8,"ほぼゼロ：相関は評価しない",transform=ax.transAxes,ha="center")
    finish(fig,"06_drag_rate.png")
    write_json(report/"data"/"summary.json",summary)
    if (output/"validation"/"resolution.json").exists():
        (report/"data"/"resolution.json").write_bytes((output/"validation"/"resolution.json").read_bytes())
    if (output/"regression_24h"/"regression.json").exists():
        (report/"data"/"regression.json").write_bytes((output/"regression_24h"/"regression.json").read_bytes())
    if (output/"old_reference.json").exists():
        (report/"data"/"old_reference.json").write_bytes((output/"old_reference.json").read_bytes())
    (report/"study.json").write_bytes((output/"study.json").read_bytes())
    protected=json.loads((output/"preserved_sha256.json").read_text(encoding="utf-8"))
    changed=[p for p,h in protected.items() if not (ROOT/p).exists()
             or hashlib.sha256((ROOT/p).read_bytes()).hexdigest()!=h]
    write_json(report/"data"/"preservation.json",dict(checked_files=len(protected),changed=changed,passed=not changed))
    assert not changed
    print({name:dict(max_error_m=m["formation"]["max_abs_error_m"],
                    J_delta_spread=m["exposure"]["J_delta_spread_kg_s_m3"])
           for name,m in summary.items()},flush=True)


def old_reference(output):
    """Only read previously saved old-pulse runs; no physical recalculation."""
    short=ROOT/"outputs/long_horizon_q_20261008/regression_24h/suite/all_low_drag_baseline"
    long=ROOT/"outputs/long_horizon_q_20261008/runs/low_60d"
    if not (short/"results.npz").exists() or not (long/"results.npz").exists():
        write_json(output/"old_reference.json",dict(available=False))
        return
    metadata=json.loads((short/"metadata.json").read_text(encoding="utf-8"))
    cfg=config_from_dict(metadata["configuration"])
    long_meta=json.loads((long/"metadata.json").read_text(encoding="utf-8"))["configuration"]
    expected=copy.deepcopy(metadata["configuration"])
    expected["simulation"]["duration_s"]=60*86400
    assert expected==long_meta
    arrays=load_case_arrays(short)
    increment=arrays["density_kg_m3"]-local_background(arrays["a_m"],cfg)
    jdelta=np.trapezoid(increment,arrays["time_s"],axis=0)
    with np.load(long/"results.npz") as f:
        maximum=float(np.max(np.abs(f["gap_error_m"])))
        final_rate=f["gap_rate_m_s"][-1]
    write_json(output/"old_reference.json",dict(
        available=True,short_run=str(short.relative_to(ROOT)),long_run=str(long.relative_to(ROOT)),
        max_error_24h_m=float(np.max(np.abs(arrays["gap_error_m"]))),
        max_error_60d_m=maximum,
        max_rate_at_24h_m_s=float(np.max(np.abs(arrays["gap_rate_m_s"][-1]))),
        max_rate_at_60d_m_s=float(np.max(np.abs(final_rate))),
        J_delta_spread_kg_s_m3=float(np.ptp(jdelta)),
        J_delta_kg_s_m3=jdelta.tolist(),
        quadrature="postprocessed trapezoidal at saved 10 s; extra increment effectively zero outside old pulse",
        note="Old width 100 km, amplitude 50%, half-cosine 1 h, rotating comoving pulse. Not matched input.",
    ))


def validate(output):
    """Short numerical-resolution check, not a rerun of the 2-day cases."""
    path=output/"validation"
    prepare_output_directory(path)
    report={}
    for name in NAMES[:2]:
        base=load_config(output/"configs"/(name+".yaml"))
        values=[]
        for step in (4.0,2.0,1.0):
            d=base.to_dict()
            d["simulation"].update(duration_s=900,integration_step_s=step,output_step_s=10)
            vs=base.derived.reference_radius_m*base.derived.reference_mean_motion_rad_s
            # Same track experiment, translated to a short validation interval.
            d["disturbance"]["center_s_m"]=(vs-base.disturbance.propagation_speed_m_s)*450
            cfg=config_from_dict(d)
            result,j=simulate_all_low_exposure(cfg)
            metrics=evaluate_motion(result.arrays(),j,cfg)
            values.append(dict(step_s=step,final_error=result.gap_error_m[-1].tolist(),
                final_rate=result.gap_rate_m_s[-1].tolist(),
                max_error_m=metrics["formation"]["max_abs_error_m"],
                max_rate_at_end_m_s=metrics["formation"]["final"]["max_abs_rate_m_s"],
                J_delta_spread_kg_s_m3=metrics["exposure"]["J_delta_spread_kg_s_m3"]))
        report[name]=dict(
            resolutions=values,
            max_error_difference_2s_vs_1s_m=float(np.max(np.abs(np.array(values[1]["final_error"])-values[2]["final_error"]))),
            max_rate_difference_2s_vs_1s_m_s=float(np.max(np.abs(np.array(values[1]["final_rate"])-values[2]["final_rate"]))),
            note="900 s validation with alignment 450 s; main runs unchanged at 2 s / output 10 s",
        )
    write_json(path/"resolution.json",report)
    print(json.dumps(report,ensure_ascii=False,indent=2))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output",type=Path,required=True)
    sub=p.add_subparsers(dest="action",required=True)
    sub.add_parser("prepare")
    sub.add_parser("validate")
    sub.add_parser("old-reference")
    r=sub.add_parser("run");r.add_argument("case",choices=NAMES)
    a=sub.add_parser("analyze");a.add_argument("--report",type=Path,required=True)
    args=p.parse_args()
    if args.action=="prepare": prepare(args.output)
    elif args.action=="run": run(args.output,args.case)
    elif args.action=="validate": validate(args.output)
    elif args.action=="old-reference": old_reference(args.output)
    else: analyze(args.output,args.report)

if __name__=="__main__":
    main()
