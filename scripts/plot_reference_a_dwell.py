"""Five Japanese figures for the seven-day Reference A dwell intervention."""
import argparse,csv,json,os
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'outputs/mpl_cache'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager


def rows(path):
    with path.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))


def values(data,key):return np.array([float(x[key]) for x in data if x[key]!=''])


def main():
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,default=ROOT/'outputs/reference_a_dwell_20261017/analysis')
    p.add_argument('--output',type=Path,default=ROOT/'docs/reference_a_dwell_150/figures');args=p.parse_args()
    src=args.input;out=args.output;out.mkdir(parents=True,exist_ok=False)
    names={f.name for f in font_manager.fontManager.ttflist}
    for font in ['Yu Gothic','Meiryo','Noto Sans CJK JP']:
        if font in names:plt.rcParams['font.family']=font;break
    plt.rcParams.update({'font.size':11,'axes.unicode_minus':False,'axes.spines.top':False,'axes.spines.right':False})
    labels=['D300','D150'];colors=['#b87732','#177f9a'];captions=['dwell 300秒（既存）','dwell 150秒（新規）']
    report=json.loads((src/'comparison.json').read_text(encoding='utf-8'));stats=report['cases']
    histories={k:rows(src/f'{k}_history_60s.csv') for k in labels}
    cycles={k:rows(src/f'{k}_cycles.csv') for k in labels}
    edges={k:rows(src/f'{k}_incident_gaps.csv') for k in labels}
    def save(fig,name):fig.savefig(out/f'{name}.png',dpi=180,bbox_inches='tight');plt.close(fig)
    def cdf(ax,x,label,color,scale=1):
        x=np.sort(np.asarray(x))*scale
        if len(x):ax.step(x,np.arange(1,len(x)+1)/len(x),where='post',label=label,color=color)
    fig,axs=plt.subplots(2,2,figsize=(13,8))
    maxlim=max(stats[k]['formation']['max_error_m'] for k in labels)*1.07
    for col,(key,caption,color) in enumerate(zip(labels,captions,colors)):
        arr=np.load(src/f'{key}_plot_series.npz')
        ax=axs[0,col];ax.plot(arr['time_s']/86400,arr['error_m'],lw=.7,alpha=.55)
        ax.set_ylim(-maxlim,maxlim);ax.set_title(caption+'：全20間隔');ax.set_ylabel('符号付き間隔誤差 [m]')
        ax.axhline(0,color='grey',lw=.7);ax.axhline(100,color='#999',ls=':');ax.axhline(-100,color='#999',ls=':')
        x=values(histories[key],'time_s')/86400
        axs[1,0].plot(x,values(histories[key],'max_error_m'),label=caption,color=color)
        axs[1,1].plot(x,values(histories[key],'instantaneous_RMS_m'),label=caption,color=color)
    axs[1,0].set_ylabel('最大絶対間隔誤差 [m]');axs[1,1].set_ylabel('瞬間RMS間隔誤差 [m]')
    axs[1,0].axhline(100,color='grey',ls=':',label='100 m評価閾値')
    for ax in axs.flat:ax.set_xlabel('経過時間 [日]');ax.grid(alpha=.2)
    for ax in axs[1]:ax.legend(fontsize=9)
    fig.suptitle('A-ref：dwell短縮による7日間の間隔誤差の変化')
    fig.tight_layout();save(fig,'01_間隔誤差の時間発展')

    fig,axs=plt.subplots(1,2,figsize=(12,5))
    for key,caption,color in zip(labels,captions,colors):
        cdf(axs[0],abs(values(edges[key],'velocity_after_m_s')),caption,color,1000)
        cdf(axs[1],values(cycles[key],'V_incident_after_m_s'),caption,color,1000)
    axs[0].set_xlabel('LOW復帰時の接続間隔の絶対速度 [mm/s]')
    axs[1].set_xlabel('LOW復帰時の接続間隔の最大絶対速度 [mm/s]')
    for ax in axs:ax.set_ylabel('累積割合');ax.grid(alpha=.2);ax.legend();ax.set_ylim(0,1.02)
    fig.suptitle('残留相対速度の分布：左へ寄るほど小さい')
    fig.text(.5,-.025,'左は「局所サイクル×接続間隔」、右は「局所サイクル」が分母。重複する間隔があり、独立な標本ではない。',ha='center',fontsize=10)
    fig.tight_layout();save(fig,'02_LOW復帰時の残留速度分布')

    fig,axs=plt.subplots(1,2,figsize=(12,5))
    for key,caption,color in zip(labels,captions,colors):
        cdf(axs[0],values(cycles[key],'pulse_g_impulse_m_s'),caption,color,1000)
        x=values(cycles[key],'high_equivalent_s');unique,counts=np.unique(x,return_counts=True)
        axs[1].scatter(unique,counts/len(x),label=caption,color=color,s=40)
    axs[0].set_xlabel('自身の1パルスの速度寄与 [mm/s]');axs[0].set_ylabel('累積割合')
    axs[1].set_xlabel('HIGH換算時間 [s]');axs[1].set_ylabel('該当サイクルの割合');axs[1].set_yscale('log')
    for ax in axs:ax.grid(alpha=.2);ax.legend()
    fig.suptitle('1回の補正量とHIGH換算時間（近傍の相殺前）')
    fig.tight_layout();save(fig,'03_パルス速度寄与と最小入力')

    fig,axs=plt.subplots(1,2,figsize=(12,5))
    for key,caption,color in zip(labels,captions,colors):
        wait=values(cycles[key],'next_wait_s')/60
        cdf(axs[0],wait,caption,color)
        quiet=rows(src/f'{key}_quiet_intervals.csv')
        cdf(axs[1],values(quiet,'actual_next_wait_s')/3600,caption,color)
    axs[0].set_xlabel('同一衛星：LOW完了→次HIGH [分]');axs[0].set_xscale('log')
    axs[1].set_xlabel('全衛星LOW→次HIGH [時間]');axs[1].set_xscale('log')
    for ax in axs:ax.set_ylabel('観測済み待ち時間の累積割合');ax.grid(alpha=.2);ax.legend();ax.set_ylim(0,1.02)
    fig.suptitle('再切替までの時間：局所区間と全LOW区間を区別')
    fig.text(.5,-.025,'次回HIGHを7日以内に観測した区間だけを描画。終了で打ち切られた区間と追跡時間を考慮した割合は表へ保存。',ha='center',fontsize=10)
    fig.tight_layout();save(fig,'04_再切替までの時間分布')

    fig,axs=plt.subplots(1,3,figsize=(13,5))
    specs=[('完了L-H-Lサイクル数','回',[stats[k]['completed_cycles'] for k in labels]),
           ('HIGH換算滞在量','衛星・時間',[stats[k]['formation']['high_equivalent_total_s']/3600 for k in labels]),
           ('全LOW比の追加平均高度損失','m',[stats[k]['formation']['extra_mean_altitude_loss_m'] for k in labels])]
    for ax,(title,unit,x) in zip(axs,specs):
        bars=ax.bar(captions,x,color=colors,width=.58)
        ax.set_title(title);ax.set_ylabel(f'7日間の合計 [{unit}]');ax.grid(axis='y',alpha=.2);ax.set_axisbelow(True)
        for b,value in zip(bars,x):ax.text(b.get_x()+b.get_width()/2,b.get_height(),f'{value:.2f}',ha='center',va='bottom')
        ax.set_ylim(0,max(x)*1.15)
    fig.suptitle('精度と合わせて評価する切替・面積・高度のコスト')
    fig.tight_layout();save(fig,'05_切替負荷と高度コスト')


if __name__=='__main__':main()
