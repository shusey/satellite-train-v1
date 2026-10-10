#!/usr/bin/env python3
"""Five Japanese research figures from the offline comparison cache."""
import argparse
import csv
import json
import os
from pathlib import Path
import sys
import tempfile
os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir())/'sat_string_matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from sat_string.io import prepare_output_directory


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--comparison',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();base=args.comparison;rows=json.loads((base/'comparison.json').read_text(encoding='utf-8'))['cases']
    fonts={f.name for f in font_manager.fontManager.ttflist}
    for name in ['Yu Gothic','Meiryo','Noto Sans CJK JP','IPAexGothic']:
        if name in fonts:plt.rcParams['font.family']=name;break
    else:raise RuntimeError('Japanese font required')
    plt.rcParams.update({'axes.unicode_minus':False,'axes.grid':True,'grid.alpha':.25,'font.size':10})
    out=prepare_output_directory(args.output)
    colors=['black','tab:orange','tab:blue','tab:green']
    def save(fig,name):
        fig.tight_layout();fig.savefig(out/name,dpi=150,bbox_inches='tight');plt.close(fig)
    for metric,ylabel,filename in [('max_error_m','時刻ごとの最大絶対間隔誤差 [m]','max_error.png'),
                                  ('rms_error_m','時刻ごとの全間隔RMS [m]','instant_rms.png')]:
        fig,axs=plt.subplots(1,2,figsize=(12,4))
        for row,color in zip(rows,colors):
            with np.load(base/(row['case']+'_plot.npz')) as data:
                for ax in axs:ax.plot(data['time_s']/86400,data[metric],label=row['label'],color=color,lw=.85)
        for ax in axs:ax.set(xlabel='開始からの時間 [日]',ylabel=ylabel);ax.legend(fontsize=8)
        axs[0].set_title('全LOWを含む全体比較');axs[1].set_ylim(0,300);axs[1].set_title('誤差300 m以下を拡大')
        if metric=='max_error_m':
            sec=axs[0].secondary_yaxis('right',functions=(lambda m:m/1000,lambda percent:percent*1000))
            sec.set_ylabel('目標間隔100 kmに対する割合 [%]')
        fig.suptitle('中央単発外乱・60日：'+ylabel.split(' [')[0]);save(fig,filename)
    fig,ax=plt.subplots(figsize=(8,4.5))
    for row,color in zip(rows,colors):
        ax.scatter(row['extra_altitude_loss_m'],row['max_error_m'],color=color,label=row['label'],s=55)
    ax.set(xlabel='全LOW比の終了時追加平均高度損失 [m]',ylabel='全60日の最大絶対間隔誤差 [m]',
        title='精度と追加軌道コスト：中央単発外乱・60日',yscale='log');ax.legend();save(fig,'accuracy_cost.png')
    with (base/'exceedance.csv').open(encoding='utf-8-sig') as f:stats=list(csv.DictReader(f))
    fig,ax=plt.subplots(figsize=(10,4.5));x=np.arange(5)
    for k,(row,color) in enumerate(zip(rows,colors)):
        vals=[100*float(r['exceedance_fraction']) for r in stats if r['case']==row['case']]
        ax.bar(x+(k-1.5)*.2,vals,width=.2,label=row['label'],color=color)
    ax.set(xticks=x,xticklabels=['50','100','200','500','1000'],xlabel='仮想許容絶対間隔誤差 [m]',
        ylabel='いずれかの間隔が許容値を超える時間割合 [%]',title='要求精度によって変わる制御必要性：60日',ylim=(0,105))
    ax.set_title('要求精度によって変わる制御必要性：60日',pad=58)
    ax.legend(loc='lower center',bbox_to_anchor=(.5,1.01),ncol=2)
    save(fig,'exceedance_fraction.png')
    fig,axs=plt.subplots(2,2,figsize=(12,7),sharex=True)
    for ax,row in zip(axs.flat,rows):
        with np.load(base/(row['case']+'_plot.npz')) as d:
            for j,gap in enumerate([1,9,10,11,12,20]):
                ax.plot(d['time_s']/86400,d['gap_error_m'][:,j],label=f'間隔{gap}（衛星{gap}–{gap+1}）',lw=.8)
        ax.set(title=row['label'],xlabel='開始からの時間 [日]',ylabel='符号付き間隔誤差 [m]');ax.legend(fontsize=7,ncol=2)
    fig.suptitle('中央・隣接・端部の代表間隔：各パネルの縦軸範囲は異なる');save(fig,'representative_gaps.png')
    print('Five Japanese figures saved',flush=True)


if __name__=='__main__':main()
