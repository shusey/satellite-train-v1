#!/usr/bin/env python3
"""Five Japanese figures for intermediate Kp comparison, saved runs only."""
import argparse
import json
import os
from pathlib import Path
import tempfile
import sys
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "sat_string_matplotlib"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from sat_string.io import prepare_output_directory


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--comparison",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args();rows=json.loads(args.comparison.read_text(encoding="utf-8"))
    names={f.name for f in font_manager.fontManager.ttflist}
    for name in ["Yu Gothic","Meiryo","Noto Sans CJK JP","IPAexGothic"]:
        if name in names:plt.rcParams['font.family']=name;break
    else:raise RuntimeError('Japanese font required')
    plt.rcParams.update({'axes.unicode_minus':False,'axes.grid':True,'grid.alpha':.25,'font.size':10})
    out=prepare_output_directory(args.output);colors={1.:'tab:blue',2.:'tab:orange',2.5:'tab:purple',3.:'tab:green'}
    def save(fig,name):
        fig.tight_layout();fig.savefig(out/name,dpi=150,bbox_inches='tight');plt.close(fig)
    longest=max(r['days'] for r in rows);selected=[r for r in rows if r['days']==longest]
    fmax,amax=plt.subplots(figsize=(11,4));frms,arms=plt.subplots(figsize=(11,4))
    for row in selected:
        provenance=json.loads((Path(row['analysis_directory'])/'provenance.json').read_text(encoding='utf-8'))
        with np.load(Path(provenance['source'])/'results.npz') as arrays:
            t=arrays['time_s'];ix=np.unique(np.r_[np.arange(0,len(t),max(1,len(t)//20000)),len(t)-1])
            e=arrays['gap_error_m'][ix];day=t[ix]/86400
        label=f"Kp＝基準の{row['kp_factor']:g}倍";color=colors[row['kp_factor']]
        amax.plot(day,np.max(np.abs(e),axis=1),label=label,color=color,lw=.9)
        arms.plot(day,np.sqrt(np.mean(e**2,axis=1)),label=label,color=color,lw=.9)
    for ax,title,ylabel in [(amax,'最大間隔誤差','時刻ごとの最大絶対間隔誤差 [m]'),
                            (arms,'瞬間RMS','時刻ごとの全間隔RMS [m]')]:
        ax.set(title=f'中央外乱・Kpのみ変更：{longest:g}日間の{title}',xlabel='開始からの時間 [日]',ylabel=ylabel)
        ax.legend(ncol=2)
    save(fmax,'max_error_history.png');save(frms,'rms_history.png')
    fig,axs=plt.subplots(1,2,figsize=(11,4))
    for day in sorted({r['days'] for r in rows}):
        group=[r for r in rows if r['days']==day];x=[r['kp_factor'] for r in group]
        for ax,key in zip(axs,['e_max_m','rms_m']):ax.plot(x,[r[key] for r in group],marker='o',label=f'{day:g}日観測')
    for ax,label in zip(axs,['全期間の最大絶対間隔誤差 [m]','全期間・全間隔RMS [m]']):
        ax.set(xlabel='Kpの基準倍率',ylabel=label,xticks=[1,2,2.5,3]);ax.legend()
    fig.suptitle('Kpと編隊精度：中央外乱・他の条件は固定');save(fig,'accuracy_vs_kp.png')
    fig,axs=plt.subplots(1,2,figsize=(11,4))
    for day in sorted({r['days'] for r in rows}):
        group=[r for r in rows if r['days']==day];x=[r['kp_factor'] for r in group]
        for ax,key in zip(axs,['transitions','extra_altitude_loss_m']):ax.plot(x,[r[key] for r in group],marker='o',label=f'{day:g}日観測')
    for ax,label in zip(axs,['姿勢遷移開始の総数 [回]','全LOW比の追加平均高度損失 [m]']):
        ax.set(xlabel='Kpの基準倍率',ylabel=label,xticks=[1,2,2.5,3]);ax.legend()
    fig.suptitle('Kpと切替・軌道コスト：中央外乱・他の条件は固定');save(fig,'cost_vs_kp.png')
    fig,axs=plt.subplots(1,2,figsize=(11,4))
    x=[r['kp_factor'] for r in selected]
    axs[0].plot(x,[r['first_response']['first_command_time_s']/3600 for r in selected],marker='o')
    axs[0].set(ylabel='最初のHIGH指令時刻 [時間]',title='初回の反応')
    for label,key,style in [('実測','actual_next_wait_s','o-'),('予測','prediction','s--')]:
        values=[(r['first_group'][key]['predicted_wait_s'] if key=='prediction' else r['first_group'][key])/86400 for r in selected]
        axs[1].plot(x,values,style,label=label)
    axs[1].set(ylabel='第1群終了から次回HIGHまで [日]',title='残留速度による再応答');axs[1].legend()
    for ax in axs:ax.set(xlabel='Kpの基準倍率',xticks=[1,2,2.5,3])
    fig.suptitle('Kpと初回反応・次回反応：有限観測内の比較');save(fig,'response_timing.png')
    print(f'Saved five figures: {out}',flush=True)


if __name__=='__main__':main()
