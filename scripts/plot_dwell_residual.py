#!/usr/bin/env python3
"""Seven Japanese figures from saved dwell runs, with no physical recomputation."""
import argparse
import json
import os
from pathlib import Path
import sys
import tempfile
os.environ.setdefault("MPLCONFIGDIR",str(Path(tempfile.gettempdir())/"sat_string_matplotlib"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from sat_string.config import load_config
from sat_string.io import load_case_arrays,load_switch_events,prepare_output_directory
from sat_string.long_horizon import decompose_q


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    manifest=json.loads(args.manifest.read_text(encoding="utf-8"))
    names={f.name for f in font_manager.fontManager.ttflist}
    for name in ["Yu Gothic","Meiryo","Noto Sans CJK JP","IPAexGothic"]:
        if name in names:
            plt.rcParams["font.family"]=name
            break
    else:raise RuntimeError("Japanese font required")
    plt.rcParams.update({"axes.unicode_minus":False,"font.size":10,"axes.grid":True,"grid.alpha":.25})
    output=prepare_output_directory(args.output)
    colors={150:"tab:blue",300:"tab:orange",600:"tab:green"}
    def save(fig,name):
        fig.tight_layout();fig.savefig(output/name,dpi=160,bbox_inches="tight");plt.close(fig)
    data=[]
    with np.load(manifest["all_low"]+"/results.npz") as f:low=f["a_m"]
    for entry in manifest["cases"]:
        arrays=load_case_arrays(entry["run"])
        report=json.loads((Path(entry["analysis"])/"residual_metrics.json").read_text(encoding="utf-8"))
        data.append((entry["dwell_s"],arrays,report,load_switch_events(entry["run"]),load_config(Path(entry["run"])/"normalized_config.yaml")))
    fig,axes=plt.subplots(3,1,figsize=(12,9),sharex=True)
    for dwell,a,r,events,cfg in data:
        stride=max(1,len(a["time_s"])//15000);ix=np.unique(np.r_[np.arange(0,len(a["time_s"]),stride),len(a["time_s"])-1]);day=a["time_s"][ix]/86400
        label=f"dwell = {dwell}秒";color=colors[dwell]
        axes[0].plot(day,np.max(np.abs(a["gap_error_m"][ix]),axis=1),label=label,color=color,lw=1)
        axes[1].plot(day,np.sqrt(np.mean(a["gap_error_m"][ix]**2,axis=1)),label=label,color=color,lw=1)
        axes[2].plot(day,np.mean(low[ix]-a["a_m"][ix],axis=1),label=label,color=color,lw=1)
    axes[0].set(title="中央外乱・dwellだけ変更：保存時系列の比較",ylabel="時刻ごとの最大絶対間隔誤差 [m]")
    axes[1].set_ylabel("時刻ごとの全間隔RMS [m]")
    axes[2].set(ylabel="全LOW比の追加平均高度損失 [m]",xlabel="開始からの時間 [日]")
    for ax in axes:ax.legend()
    save(fig,"histories.png")

    fig,axes=plt.subplots(3,1,figsize=(12,9),sharex=True)
    for ax,(dwell,a,r,events,cfg) in zip(axes,data):
        for target,marker,color,label in [("H","^","tab:red","LOW→HIGH指令"),("L","v","tab:blue","HIGH→LOW指令")]:
            selected=[e for e in events if e["to_mode"]==target and e["command_time"]<=a["time_s"][-1]]
            ax.scatter([e["command_time"]/86400 for e in selected],[e["satellite_id"] for e in selected],s=8,marker=marker,color=color,label=label,alpha=.65)
        ax.set(title=f"dwell = {dwell}秒：指令と遷移開始は同時（遅延0秒）",ylabel="衛星番号",yticks=range(1,22,2));ax.legend(loc="upper right")
    axes[-1].set_xlabel("開始からの時間 [日]");save(fig,"switching.png")

    fig,axes=plt.subplots(3,1,figsize=(12,9),sharex=True)
    for dwell,a,r,events,cfg in data:
        local=[x for x in r["local_cycles"] if x["satellite_id"]==12]
        n=[x["local_cycle_number"] for x in local];label=f"dwell = {dwell}秒";color=colors[dwell]
        axes[0].plot(n,[x["V_res_seeds_m_s"]*1000 for x in local],"o-",ms=3,label=label,color=color)
        axes[1].plot(n,[x["own_q_D_m_s2"]/1e-7 for x in local],"o-",ms=3,label=label,color=color)
        axes[2].plot(n,[x["actual_next_wait_s"]/86400 if x["actual_next_wait_s"] is not None else np.nan for x in local],"o-",ms=3,label=label,color=color)
    axes[0].set(title="衛星12の局所L-H-L終了状態（各条件の周囲状態・時刻は異なる）",ylabel="接続間隔の残留速度最大値 [mm/s]")
    axes[1].set_ylabel("衛星12の残留D項 [10^-7 m/s²]")
    axes[2].set(ylabel="LOW完了から同衛星の次回指令まで [日]",xlabel="衛星12の局所L-H-L番号")
    for ax in axes:ax.legend()
    save(fig,"local_residuals.png")

    fig,axes=plt.subplots(2,2,figsize=(12,8))
    for dwell,a,r,events,cfg in data:
        quiet=r["quiet_cycles"];n=[x["response_group_number"] for x in quiet];label=f"dwell = {dwell}秒";color=colors[dwell]
        axes[0,0].plot(n,[x["V_res_all_m_s"]*1000 for x in quiet],".-",label=label,color=color)
        axes[0,1].plot(n,[x["Q_D_res_all_m_s2"]/1e-7 for x in quiet],".-",label=label,color=color)
        axes[1,0].plot(n,[x["E_n_active_max_error_m"] for x in quiet],".-",label=label,color=color)
        axes[1,1].plot(n,[x["extra_loss_increment_m"] for x in quiet],".-",label=label,color=color)
    axes[0,0].set_ylabel("全間隔の残留速度最大値 [mm/s]")
    axes[0,1].set_ylabel("全衛星の残留D項絶対最大値 [10^-7 m/s²]")
    axes[1,0].set_ylabel("活動区間内の最大間隔誤差 [m]")
    axes[1,1].set_ylabel("活動区間の追加平均高度損失増分 [m]")
    for ax in axes.flat:ax.set_xlabel("全LOWへ戻る活動区間番号");ax.legend()
    fig.suptitle("全LOW境界での残留状態：参加集合が異なるため番号の直接対応はない")
    save(fig,"quiet_residuals.png")

    fig,axes=plt.subplots(1,2,figsize=(12,5))
    for dwell,a,r,events,cfg in data:
        finite=[x for x in r["quiet_cycles"] if not x["right_censored"] and x["prediction"]["predicted_wait_s"] is not None]
        valid=[x for x in finite if x["prediction_valid_on_observed_window"]]
        invalid=[x for x in finite if not x["prediction_valid_on_observed_window"]]
        axes[0].scatter([x["prediction"]["predicted_wait_s"]/86400 for x in valid],[x["actual_next_wait_s"]/86400 for x in valid],s=20,label=f"dwell = {dwell}秒",color=colors[dwell])
        if invalid:axes[0].scatter([x["prediction"]["predicted_wait_s"]/86400 for x in invalid],[x["actual_next_wait_s"]/86400 for x in invalid],marker="x",color="gray",label=f"{dwell}秒：仮定不成立")
        axes[1].scatter([x["completion_time_s"]/86400 for x in valid],[x["relative_prediction_error"]*100 for x in valid],s=15,color=colors[dwell],label=f"dwell = {dwell}秒")
    lim=max(axes[0].get_xlim()[1],axes[0].get_ylim()[1]);axes[0].plot([0,lim],[0,lim],"k--",label="予測＝実測")
    axes[0].set(xlabel="全LOW境界からの予測待ち時間 [日]",ylabel="実測待ち時間 [日]",title="次のLOW→HIGH指令：全LOW境界からの予測")
    axes[1].set(xlabel="全LOW境界の時刻 [日]",ylabel="（実測−予測）／実測 [%]",title="q_D変化10%以内・介在指令なしの区間")
    for ax in axes:ax.legend()
    save(fig,"prediction.png")

    fig,axes=plt.subplots(3,2,figsize=(13,10))
    for row,(dwell,a,r,events,cfg) in enumerate(data):
        local=next(x for x in r["local_cycles"] if x["satellite_id"]==12 and x["local_cycle_number"]==1)
        start=local["high_slew_start_s"];mask=(a["time_s"]>=start-60)&(a["time_s"]<=start+1600)
        p,d=decompose_q(a["gap_error_m"][mask],a["gap_rate_m_s"][mask],cfg.derived.controller_kp_s2_inv,cfg.derived.controller_kd_s_inv)
        time=a["time_s"][mask]-start;ax=axes[row,0]
        for values,label in [(p[:,11],"P項"),(d[:,11],"D項"),(p[:,11]+d[:,11],"判定量q")]:ax.plot(time,values/1e-7,label=label)
        ax.axhline(cfg.controller.threshold_on_m_s2/1e-7,color="black",ls="--",label="HIGH閾値")
        ax.axhline(cfg.controller.threshold_off_m_s2/1e-7,color="gray",ls=":",label="LOW閾値")
        ax.set(title=f"dwell = {dwell}秒：最初の衛星12のHIGHパルス",ylabel="制御判定量 [10^-7 m/s²]")
        ar=axes[row,1];ar.plot(time,a["area_m2"][mask,11],label="衛星12の面積")
        for key,label in [("high_complete_s","HIGH完了"),("low_command_s","LOW指令"),("low_complete_s","LOW完了")]:
            ar.axvline(local[key]-start,ls="--",label=label)
        ar.set(title=f"HIGH等価時間 {local['high_equivalent_time_s']:.0f}秒",ylabel="断面積 [m²]")
        for panel in [ax,ar]:panel.legend(fontsize=8,ncol=2);panel.set_xlabel("最初のHIGH開始から [秒]")
    save(fig,"first_pulse.png")

    fig,axes=plt.subplots(3,1,figsize=(12,9))
    for ax,(dwell,a,r,events,cfg) in zip(axes,data):
        if dwell==300:
            endpoint=r["group_cycles"][0]
        else:
            previous=json.loads((Path(manifest["prefix_analysis_root"])/f"dwell_{dwell}_30d/residual_metrics.json").read_text(encoding="utf-8"))
            source=previous["quiet_cycles"][-1]["completion_time_s"]
            endpoint=next(x for x in r["quiet_cycles"] if x["completion_time_s"]==source)
        j=endpoint["predicted_satellite_id"]-1
        mask=(a["time_s"]>=endpoint["state_sample_time_s"])&(a["time_s"]<=endpoint["actual_next_command_s"])
        e=a["gap_error_m"][mask];v=a["gap_rate_m_s"][mask]
        stride=max(1,len(e)//5000);ix=np.unique(np.r_[np.arange(0,len(e),stride),len(e)-1])
        p,d=decompose_q(e[ix],v[ix],cfg.derived.controller_kp_s2_inv,cfg.derived.controller_kd_s_inv)
        time=a["time_s"][mask][ix]/86400
        for values,label in [(p[:,j],"P項"),(d[:,j],"D項"),(p[:,j]+d[:,j],"判定量q")]:ax.plot(time,values/1e-7,label=label)
        for value,style,label in [(cfg.controller.threshold_on_m_s2,"--","HIGH閾値"),(cfg.controller.threshold_off_m_s2,":","LOW閾値")]:
            ax.axhline(value/1e-7,ls=style,color="gray",label=label)
        ax.set(title=f"dwell = {dwell}秒・衛星{j+1}：全LOW境界から実際の次回指令まで",xlabel="開始からの時間 [日]",ylabel="制御判定量 [10^-7 m/s²]");ax.legend(ncol=5)
    save(fig,"residual_drift.png")


if __name__=="__main__":main()
