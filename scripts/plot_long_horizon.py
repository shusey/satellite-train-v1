#!/usr/bin/env python3
"""Four Japanese scientific figures from saved runs; no simulation."""
import argparse
import json
import os
from pathlib import Path
import sys
import tempfile
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir())/"sat_string_matplotlib"))
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
    p=argparse.ArgumentParser(description=__doc__)
    for name in ["run","all-low","analysis","output"]:
        p.add_argument("--"+name,type=Path,required=True)
    args=p.parse_args()
    names={f.name for f in font_manager.fontManager.ttflist}
    for name in ["Yu Gothic","Meiryo","Noto Sans CJK JP","IPAexGothic"]:
        if name in names:
            plt.rcParams["font.family"]=name
            break
    else:
        raise RuntimeError("Japanese font required; install a Japanese font before plotting")
    plt.rcParams.update({"axes.unicode_minus":False,"font.size":10,"axes.grid":True,"grid.alpha":.25})
    output=prepare_output_directory(args.output)
    a=load_case_arrays(args.run)
    with np.load(args.all_low/"results.npz") as f:
        extra=np.mean(f["a_m"]-a["a_m"],axis=1)
    report=json.loads((args.analysis/"waves.json").read_text(encoding="utf-8"))
    waves=report["waves"];cfg=load_config(args.run/"normalized_config.yaml")
    day=a["time_s"]/86400
    def save(fig,name):
        fig.savefig(output/name,dpi=160,bbox_inches="tight");plt.close(fig)
    fig,axes=plt.subplots(3,1,figsize=(12,9),sharex=True)
    axes[0].plot(day,np.max(np.abs(a["gap_error_m"]),axis=1))
    axes[0].set_ylabel("最大絶対間隔誤差 [m]")
    axes[0].set_title("中央外乱・固定基準制御：長時間応答（各波の活動期間を着色）")
    axes[1].plot(day,np.max(np.abs(a["gap_rate_m_s"]),axis=1)*1000)
    axes[1].set_ylabel("最大絶対間隔速度 [mm/s]")
    axes[2].plot(day,extra);axes[2].set_ylabel("全LOW比の追加平均高度損失 [m]")
    axes[2].set_xlabel("外乱前の開始時刻からの時間 [日]")
    for w in waves:
        for ax in axes:
            ax.axvspan(w["start_time_s"]/86400,(w["end_time_s"] or a["time_s"][-1])/86400,alpha=.14,color="orange")
        axes[0].text(w["start_time_s"]/86400,axes[0].get_ylim()[1]*.9,str(w["wave_number"]))
    fig.tight_layout();save(fig,"long_time_history.png")

    events=load_switch_events(args.run)
    fig,axes=plt.subplots(2,1,figsize=(12,9))
    for target,color,marker,label in [("H","tab:red","^","LOW→HIGH開始"),("L","tab:blue","v","HIGH→LOW開始")]:
        selected=[r for r in events if r["to_mode"]==target and r["start_time"]<=a["time_s"][-1]]
        axes[0].scatter([r["start_time"]/86400 for r in selected],[r["satellite_id"] for r in selected],s=15,color=color,marker=marker,label=label)
    axes[0].set(xlabel="時間 [日]",ylabel="衛星番号",title="姿勢遷移開始イベント（指令と開始は遅延0 sで一致）",yticks=range(1,22,2))
    axes[0].legend()
    for w in waves:
        arrivals=np.asarray([np.nan if x is None else x for x in w["first_l_to_h_times_s"]])
        axes[1].plot(range(1,22),(arrivals-w["start_time_s"])/3600,"o-",label=f"第{w['wave_number']}波")
    axes[1].set(xlabel="衛星番号",ylabel="その波の開始からLOW→HIGH開始まで [h]",title="伝播パターン：中央の遅い切替と両側の発生点",xticks=range(1,22,2))
    axes[1].legend(ncol=4);fig.tight_layout();save(fig,"switching_events.png")

    n=[w["wave_number"] for w in waves]
    fig,axes=plt.subplots(2,2,figsize=(12,8))
    axes[0,0].plot(n,[w["E_n_active_max_error_m"] for w in waves],"o-",label="活動期間の最大誤差（未完了含む）")
    axes[0,0].plot(n,[w["error_at_observed_end_m"] if w["complete"] else np.nan for w in waves],"s--",label="完了した波の終了時誤差")
    axes[0,0].set_ylabel("間隔誤差 [m]");axes[0,0].legend()
    axes[0,1].plot(n[:-1],[w["next_start_interval_s"]/86400 for w in waves[:-1]],"o-")
    axes[0,1].set_ylabel("次波までの開始時刻間隔 [日]")
    axes[1,0].plot(n,[w["extra_loss_active_increment_m"] for w in waves],"o-")
    axes[1,0].set_ylabel("活動期間の追加平均高度損失増分 [m]")
    d=report["cycle_end_differences"]
    axes[1,1].plot([r["to_wave"] for r in d],[r["error_l2_difference_m"] for r in d],"o-")
    axes[1,1].set_ylabel("前波終了状態との誤差差ベクトル長 [m]")
    for ax in axes.flat:ax.set_xlabel("波番号");ax.set_xticks(n)
    fig.suptitle("波ごとの振幅・間隔・コスト・終了状態（有限観測）")
    fig.tight_layout();save(fig,"wave_trends.png")

    with np.load(args.analysis/"detail_7_to_10_days.npz") as f:
        detail={k:f[k] for k in f.files}
    t=detail["time_s"];j=8
    seed=next(r for r in report["triggers"] if r["wave_number"]==2 and r["satellite_id"]==9)
    fig,axes=plt.subplots(3,2,figsize=(14,12))
    for edge in [7,8]:
        axes[0,0].plot(t/86400,detail["gap_error_m"][:,edge],label=f"間隔{edge+1}（衛星{edge+1}–{edge+2}）")
        axes[0,1].plot(t/86400,detail["gap_rate_m_s"][:,edge]*1000,label=f"間隔{edge+1}")
    axes[0,0].set_ylabel("間隔誤差 [m]");axes[0,1].set_ylabel("間隔速度 [mm/s]")
    colors=["tab:orange","tab:green","tab:blue"]
    for ax in [axes[1,0],axes[2,0]]:
        zoom=ax is axes[2,0]
        mask=(t>=seed["command_time_s"]-600)&(t<=seed["return_command_s"]+240) if zoom else np.ones(len(t),dtype=bool)
        x=t[mask]-seed["command_time_s"] if zoom else t[mask]/86400
        for key,label,color in zip(["q_P_m_s2","q_D_m_s2","q_m_s2"],["P項","D項","合計q"],colors):
            ax.plot(x,detail[key][mask,j]/1e-7,label=label,color=color)
        ax.axhline(cfg.controller.threshold_on_m_s2/1e-7,color="red",ls="--",label="HIGH閾値 h_on")
        ax.axhline(cfg.controller.threshold_off_m_s2/1e-7,color="purple",ls=":",label="LOW閾値 h_off")
        ax.set_ylabel("衛星9の制御判定量 [10^-7 m/s²]")
    for sat in [9,13,11]:
        axes[1,1].plot(t/86400,detail["q_m_s2"][:,sat-1]/1e-7,label=f"衛星{sat}")
    axes[1,1].axhline(cfg.controller.threshold_on_m_s2/1e-7,color="red",ls="--",label="h_on")
    axes[1,1].axhline(cfg.controller.threshold_off_m_s2/1e-7,color="purple",ls=":",label="h_off")
    axes[1,1].set_ylabel("制御判定量 q [10^-7 m/s²]")
    mask=(t>=seed["command_time_s"]-600)&(t<=seed["return_command_s"]+240)
    x=t[mask]-seed["command_time_s"]
    axes[2,1].plot(x,detail["area_m2"][mask,j],color="black",label="実効面積")
    axes[2,1].set_ylabel("衛星9の面積 [m²]")
    mode_axis=axes[2,1].twinx()
    mode_axis.step(x,detail["mode"][mask,j],where="post",alpha=.45,color="tab:blue")
    mode_axis.set_yticks([0,1,2,3],labels=["LOW","LOW→HIGH","HIGH","HIGH→LOW"])
    slew=cfg.controller.slew_time_s
    return_delay=slew+cfg.controller.dwell_time_s
    for ax in axes[2,:]:
        ax.axvline(0,color="red",ls="--",label="指令・姿勢遷移開始")
        ax.axvline(slew,color="gray",ls=":",label="HIGH姿勢完了")
        ax.axvspan(slew,return_delay,color="gray",alpha=.12,label="待機（dwell）期間")
        ax.axvline(return_delay,color="purple",ls="--",label="LOWへの指令・開始")
        ax.axvline(return_delay+slew,color="gray",ls="-.",label="LOW姿勢完了")
        ax.set_xlabel("第2波の衛星9の指令からの時間 [s]")
    for ax in axes[:2,:].flat:ax.set_xlabel("時間 [日]")
    for ax in axes.flat:ax.legend(fontsize=8,loc="best")
    fig.suptitle("7–10日の二次応答：P項の蓄積による再横断と、D項・dwellによる戻り")
    fig.tight_layout();save(fig,"q_cause_7_to_10_days.png")

    if a["time_s"][-1]>45*86400:
        late=a["time_s"]>=45*86400
        qp,qd=decompose_q(a["gap_error_m"][late],a["gap_rate_m_s"][late],
            cfg.derived.controller_kp_s2_inv,cfg.derived.controller_kd_s_inv)
        fig,axes=plt.subplots(3,1,figsize=(13,10),sharex=True)
        for ax,sat in zip(axes,[9,11,13]):
            for values,label,color in zip([qp,qd,a["q_m_s2"][late]],["P項","D項","合計q"],colors):
                ax.plot(day[late],values[:,sat-1]/1e-7,label=label,color=color)
            ax.axhline(cfg.controller.threshold_on_m_s2/1e-7,color="red",ls="--",label="h_on")
            ax.axhline(cfg.controller.threshold_off_m_s2/1e-7,color="purple",ls=":",label="h_off")
            ax.set_ylabel(f"衛星{sat}の判定量 [10^-7 m/s²]");ax.legend(ncol=5,fontsize=9)
        axes[0].set_title("45日以降の再応答：P項・D項と左右非同期の変化")
        axes[-1].set_xlabel("時間 [日]");fig.tight_layout();save(fig,"q_late_response.png")


if __name__=="__main__":main()
