#!/usr/bin/env python3
"""Nine Japanese figures, using saved Kd runs and residual endpoints only."""
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
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--manifest",type=Path,required=True);p.add_argument("--output",type=Path,required=True)
    args=p.parse_args();manifest=json.loads(args.manifest.read_text(encoding="utf-8"))
    names={f.name for f in font_manager.fontManager.ttflist}
    for name in ["Yu Gothic","Meiryo","Noto Sans CJK JP","IPAexGothic"]:
        if name in names:plt.rcParams["font.family"]=name;break
    else:raise RuntimeError("Japanese font required")
    plt.rcParams.update({"axes.unicode_minus":False,"font.size":10,"axes.grid":True,"grid.alpha":.25})
    out=prepare_output_directory(args.output);colors={.3:"tab:blue",1.:"tab:orange",3.:"tab:green"}
    def save(fig,name):
        fig.tight_layout();fig.savefig(out/name,dpi=150,bbox_inches="tight");plt.close(fig)
    cases=[]
    with np.load(Path(manifest["all_low"])/"results.npz") as f:low=f["a_m"]
    for entry in manifest["cases"]:
        factor=entry["kd_factor"];directory=Path(entry["analysis"])
        report=json.loads((directory/"residual_metrics.json").read_text(encoding="utf-8"))
        first=json.loads((directory/"first_response.json").read_text(encoding="utf-8"))
        with np.load(directory/"cycle_end_states.npz") as f:states={k:f[k] for k in f.files}
        cases.append((factor,load_case_arrays(entry["run"]),report,load_switch_events(entry["run"]),
            load_config(Path(entry["run"])/"normalized_config.yaml"),states,first))
    fig,axs=plt.subplots(3,1,figsize=(12,9),sharex=True)
    for factor,a,r,events,cfg,states,first in cases:
        # Plot subsampling only; all numerical metrics use every saved 10-second sample.
        ix=np.unique(np.r_[np.arange(0,len(a["time_s"]),max(1,len(a["time_s"])//20000)),len(a["time_s"])-1])
        day=a["time_s"][ix]/86400;label=f"Kd = 基準の{factor:g}倍"
        axs[0].plot(day,np.max(np.abs(a["gap_error_m"][ix]),axis=1),label=label,color=colors[factor],lw=.9)
        axs[1].plot(day,np.sqrt(np.mean(a["gap_error_m"][ix]**2,axis=1)),label=label,color=colors[factor],lw=.9)
        axs[2].plot(day,np.mean(low[ix]-a["a_m"][ix],axis=1),label=label,color=colors[factor],lw=.9)
    axs[0].set(title="中央外乱・Kdのみ変更（dwell 300秒）：保存時系列の比較",ylabel="時刻ごとの最大絶対間隔誤差 [m]")
    axs[1].set_ylabel("時刻ごとの全間隔RMS [m]")
    axs[2].set(ylabel="全LOW比の追加平均高度損失 [m]",xlabel="開始からの時間 [日]")
    for ax in axs:ax.legend()
    zoom=axs[2].inset_axes([.60,.12,.36,.42])
    for factor,a,r,events,cfg,states,first in cases:
        if factor==3:continue
        ix=np.unique(np.r_[np.arange(0,len(a["time_s"]),max(1,len(a["time_s"])//10000)),len(a["time_s"])-1])
        zoom.plot(a["time_s"][ix]/86400,np.mean(low[ix]-a["a_m"][ix],axis=1),color=colors[factor],lw=1)
    zoom.set_title("0.3倍・1倍の拡大",fontsize=9);zoom.set_ylabel("追加損失 [m]",fontsize=8)
    zoom.tick_params(labelsize=8)
    save(fig,"histories.png")

    fig,axs=plt.subplots(len(cases),1,figsize=(12,9),sharex=True,squeeze=False)
    for ax,(factor,a,r,events,cfg,states,first) in zip(axs[:,0],cases):
        for target,marker,color,label in [("H","^","tab:red","LOW→HIGH指令"),("L","v","tab:blue","HIGH→LOW指令")]:
            selected=[e for e in events if e["to_mode"]==target and e["command_time"]<=a["time_s"][-1]]
            ax.scatter([e["command_time"]/86400 for e in selected],[e["satellite_id"] for e in selected],
                s=7,marker=marker,color=color,label=label,alpha=.6)
        ax.set(title=f"Kd = 基準の{factor:g}倍：指令＝遷移開始（遅延0秒）",ylabel="衛星番号",yticks=range(1,22,2));ax.legend()
    axs[-1,0].set_xlabel("開始からの時間 [日]");save(fig,"switching.png")

    fig,axs=plt.subplots(3,1,figsize=(12,9),sharex=True)
    for factor,a,r,events,cfg,states,first in cases:
        local=[x for x in r["local_cycles"] if x["satellite_id"]==12]
        n=[x["local_cycle_number"] for x in local];kw=dict(color=colors[factor],label=f"Kd = {factor:g}倍")
        axs[0].plot(n,[x["V_res_all_m_s"]*1000 for x in local],"o-",ms=3,**kw)
        axs[1].plot(n,[x["V_res_rms_all_m_s"]*1000 for x in local],"o-",ms=3,**kw)
        axs[2].plot(n,[x["actual_next_wait_s"]/86400 if x["actual_next_wait_s"] is not None else np.nan for x in local],"o-",ms=3,**kw)
    axs[0].set(title="衛星12のL-H-L完了時：各条件で周囲状態・時刻は異なる",ylabel="全間隔の残留速度最大値 [mm/s]")
    axs[1].set_ylabel("全間隔の残留速度RMS [mm/s]")
    axs[2].set(ylabel="LOW完了から同衛星の次回指令まで [日]",xlabel="衛星12の局所cycle番号（対数軸）")
    for ax in axs:ax.set_xscale("log");ax.legend()
    save(fig,"local_residuals.png")

    fig,axs=plt.subplots(2,2,figsize=(12,8))
    for factor,a,r,events,cfg,states,first in cases:
        rows=r["quiet_cycles"];n=[x["response_group_number"] for x in rows]
        for ax,key,scale in [(axs[0,0],"V_res_all_m_s",1000),(axs[0,1],"V_res_rms_all_m_s",1000),
                (axs[1,0],"actual_next_wait_s",1/86400),(axs[1,1],"worsening_weight_fraction",1)]:
            ax.plot(n,[x[key]*scale if x[key] is not None else np.nan for x in rows],".-",ms=3,color=colors[factor],label=f"Kd = {factor:g}倍")
    for ax in axs.flat:ax.set_xscale("log");ax.set_xlabel("全LOWへ戻る活動区間番号（対数軸）");ax.legend()
    axs[0,0].set_ylabel("全間隔の残留速度最大値 [mm/s]");axs[0,1].set_ylabel("全間隔の残留速度RMS [mm/s]")
    axs[1,0].set_ylabel("全LOW境界から次回指令まで [日]");axs[1,1].set_ylabel("誤差悪化方向の |e・速度| 加重割合")
    fig.suptitle("全LOW境界の比較：参加集合が異なるため番号に直接対応はない")
    save(fig,"quiet_residuals.png")

    fig,axs=plt.subplots(2,len(cases),figsize=(13,8))
    for column,(factor,a,r,events,cfg,states,first) in enumerate(cases):
        for row,section,title in [(0,"group_cycles","最初の時間群終了"),(1,"quiet_cycles","最初の全LOW境界")]:
            ax=axs[row,column];records=r[section]
            if records:
                endpoint=records[0];k=endpoint["end_state_row"];e=states["gap_error_m"][k];v=states["gap_rate_m_s"][k]*1000
                points=ax.scatter(e,v,c=np.arange(1,len(e)+1),cmap="viridis",s=35)
                ax.set_title(f"Kd {factor:g}倍・{title}\n{endpoint['state_sample_time_s']/86400:.4f}日")
                fig.colorbar(points,ax=ax,label="間隔番号")
            else:ax.set_title(f"Kd {factor:g}倍：該当境界なし")
            ax.axhline(0,color="gray",lw=.8);ax.axvline(0,color="gray",lw=.8)
            ax.set(xlabel="間隔誤差 [m]",ylabel="符号付き間隔速度 [mm/s]")
    fig.suptitle("終了状態の分布：第1・第3象限は誤差が増す方向（時間群は固定の分類間隔）")
    save(fig,"error_velocity.png")

    fig,axs=plt.subplots(1,2,figsize=(12,5));extent=0.
    for factor,a,r,events,cfg,states,first in cases:
        rows=[x for x in r["quiet_cycles"] if not x["right_censored"] and x["prediction"]["predicted_wait_s"] is not None]
        valid=[x for x in rows if x["prediction_valid_on_observed_window"]];invalid=[x for x in rows if not x["prediction_valid_on_observed_window"]]
        px=[x["prediction"]["predicted_wait_s"]/86400 for x in valid];ay=[x["actual_next_wait_s"]/86400 for x in valid]
        extent=max([extent]+px+ay)
        axs[0].scatter(px,ay,s=18,color=colors[factor],label=f"Kd = {factor:g}倍：有効")
        if invalid:axs[0].scatter([x["prediction"]["predicted_wait_s"]/86400 for x in invalid],[x["actual_next_wait_s"]/86400 for x in invalid],marker="x",s=18,color=colors[factor],alpha=.5,label=f"{factor:g}倍：仮定不成立")
        axs[1].scatter([x["completion_time_s"]/86400 for x in valid],[x["relative_prediction_error"]*100 for x in valid],s=16,color=colors[factor],label=f"Kd = {factor:g}倍")
    axs[0].plot([0,extent],[0,extent],"k--",lw=.8,label="予測＝実測")
    axs[0].set(xlabel="全LOW境界からの予測待ち時間 [日]",ylabel="実測待ち時間 [日]",title="一定D項の予測と実測")
    axs[1].set(xlabel="予測に使った全LOW境界 [日]",ylabel="（実測－予測）／実測待ち時間 [%]",title="有効区間のみの相対誤差")
    for ax in axs:ax.legend(fontsize=8)
    save(fig,"predictions.png")

    fig,axs=plt.subplots(3,1,figsize=(12,10))
    for ax,(factor,a,r,events,cfg,states,first) in zip(axs,cases):
        seed=first["first_satellite_ids"][0];local=next((x for x in r["local_cycles"] if x["satellite_id"]==seed and x["local_cycle_number"]==1),None)
        onset=first["first_command_time_s"];finish=local["low_complete_s"] if local else onset+2000
        ix=np.flatnonzero((a["time_s"]>=onset-1000)&(a["time_s"]<=finish+1000))
        pterm,dterm=decompose_q(a["gap_error_m"][ix],a["gap_rate_m_s"][ix],cfg.derived.controller_kp_s2_inv,cfg.derived.controller_kd_s_inv)
        x=(a["time_s"][ix]-onset)/3600;j=seed-1
        for y,label,color in [(pterm[:,j],"位置項 q_P","tab:blue"),(dterm[:,j],"速度項 q_D","tab:orange"),(pterm[:,j]+dterm[:,j],"合計 q","black")]:
            ax.plot(x,y/1e-7,label=label,color=color,lw=1)
        ax.axhline(cfg.controller.threshold_on_m_s2/1e-7,color="tab:red",ls="--",label="HIGH閾値 h_on")
        ax.axhline(cfg.controller.threshold_off_m_s2/1e-7,color="tab:green",ls=":",label="LOW閾値 h_off")
        if local:
            for key,label,color in [("high_complete_s","HIGH完了","gray"),("low_command_s","LOW指令","purple"),("low_complete_s","LOW完了","brown")]:
                ax.axvline((local[key]-onset)/3600,color=color,ls=":",lw=.8,label=label)
        ax.set(title=f"Kd = {factor:g}倍・初回起点の衛星{seed}",ylabel="制御判定量 [10^-7 m/s²]")
        ax.legend(ncol=4,fontsize=8,loc="best")
    axs[-1].set_xlabel("初回LOW→HIGH指令からの時間 [時間]")
    save(fig,"first_pulse_q.png")

    fig,axs=plt.subplots(1,4,figsize=(16,5),sharey=True)
    categories=[("P_dominant","P項が2倍以上","tab:blue"),("D_dominant","D項が2倍以上","tab:orange"),("mixed","混合","tab:green")]
    for ax,key,title in [(axs[0],"command_contribution_counts","全HIGH指令"),
            (axs[1],"quiet_onset_command_contribution_counts","全LOW後の活動開始"),
            (axs[2],"later_in_packet_command_contribution_counts","同じ活動区間の後続指令"),
            (axs[3],"temporal_group_onset_contribution_counts","固定分類間隔で分けた時間群の起点")]:
        for i,(factor,a,r,events,cfg,states,first) in enumerate(cases):
            counts=r["summary"][key];total=sum(counts.values());bottom=0
            for kind,label,color in categories:
                height=counts[kind]/total if total else 0
                ax.bar(i,height,bottom=bottom,color=color,label=label if i==0 else None);bottom+=height
            ax.text(i,1.02,f"{total}指令",ha="center",fontsize=9)
        ax.set(title=title,xlabel="Kdの基準倍率",xticks=range(len(cases)),xticklabels=[f"{x[0]:g}" for x in cases],ylim=(0,1.13));ax.legend(fontsize=8,loc="lower right")
    axs[0].set_ylabel("絶対値比で分類した指令割合")
    fig.suptitle("P/Dの大きさの分類：符号の相殺を含むため因果寄与率ではない")
    save(fig,"pd_contributions.png")
    fig,axs=plt.subplots(1,2,figsize=(12,5))
    for factor,a,r,events,cfg,states,first in cases:
        arrivals=first["first_arrivals_s_by_satellite"]
        participants=first["first_group"]["participant_satellite_ids"]
        stop=first["first_group"]["observed_end_s"]
        ids=sorted(int(k) for k in arrivals if int(k) in participants and arrivals[k]<=stop)
        axs[0].plot(ids,[arrivals[str(j)]/3600 for j in ids],"o-",ms=3,color=colors[factor],label=f"Kd = {factor:g}倍")
        # First local return per satellite; nearby actuation can still be active.
        local=[x for x in first["first_local_cycles"] if x["satellite_id"] in participants and x["high_command_s"]<=stop]
        ids=[x["satellite_id"] for x in local]
        axs[1].scatter(ids,[x["V_res_seeds_m_s"]*1000 for x in local],s=25,color=colors[factor],label=f"Kd = {factor:g}倍")
    axs[0].set(title="最初の時間群に属する各衛星の初回HIGH指令",xlabel="衛星番号",ylabel="開始からの指令時刻 [時間]",xticks=range(1,22,2))
    axs[1].set(title="各衛星の初回LOW復帰：接続間隔の残留速度",xlabel="衛星番号",ylabel="接続間隔の残留速度最大値 [mm/s]",xticks=range(1,22,2))
    for ax in axs:ax.legend()
    save(fig,"first_response.png")
    print(f"Saved nine figures: {out}",flush=True)


if __name__=="__main__":main()
