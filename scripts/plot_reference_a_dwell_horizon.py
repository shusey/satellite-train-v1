"""Five Japanese figures comparing reused D300 and checkpoint-extended D150."""
import argparse,json,sys,os
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'outputs/mpl_cache'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from analyze_reference_pulses import read_csv


def col(rows,key):return np.asarray([float(x[key]) for x in rows if x[key]!=''])


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=ROOT/'outputs/reference_a_dwell_60d_20261017')
    ap.add_argument('--output',type=Path,default=ROOT/'docs/reference_a_dwell_60d/figures');args=ap.parse_args()
    src=args.root/'analysis_60d';out=args.output;out.mkdir(parents=True,exist_ok=False)
    names={f.name for f in font_manager.fontManager.ttflist}
    for name in ['Yu Gothic','Meiryo','Noto Sans CJK JP']:
        if name in names:plt.rcParams['font.family']=name;break
    plt.rcParams.update({'font.size':11,'axes.unicode_minus':False,'axes.spines.top':False,'axes.spines.right':False})
    labels=['D300','D150'];captions=['dwell 300秒（既存）','dwell 150秒（延長）'];colors=['#b87732','#177f9a']
    history={k:read_csv(src/f'{k}_history_60s.csv') for k in labels}
    summary=json.loads((src/'comparison.json').read_text(encoding='utf-8'))['cases']
    def save(fig,name):fig.savefig(out/f'{name}.png',dpi=180,bbox_inches='tight');plt.close(fig)
    fig,axs=plt.subplots(2,1,figsize=(12,8),sharex=True)
    for key,caption,color in zip(labels,captions,colors):
        rows=history[key];t=col(rows,'time_s')/86400
        axs[0].plot(t,col(rows,'max_error_m'),label=caption,color=color,lw=1)
        axs[1].plot(t,col(rows,'instantaneous_RMS_m'),label=caption,color=color,lw=1)
    axs[0].axhline(100,color='grey',ls=':',label='100 m評価閾値')
    axs[0].axhline(500,color='grey',ls='--',label='500 m評価閾値')
    axs[0].set_ylabel('最大絶対間隔誤差 [m]');axs[1].set_ylabel('瞬間RMS間隔誤差 [m]');axs[1].set_xlabel('経過時間 [日]')
    for ax in axs:ax.grid(alpha=.2);ax.legend(ncol=2,fontsize=9)
    fig.suptitle('A-ref：60日間の精度（初期半径分散0.5 m、B差なし）')
    fig.tight_layout();save(fig,'01_60日の間隔誤差')

    comparisons={7:json.loads((ROOT/'docs/reference_a_dwell_150/data/comparison.json').read_text())['cases']}
    for day in [30,60]:comparisons[day]=json.loads((args.root/f'analysis_{day}d/comparison.json').read_text())['cases']
    fig,axs=plt.subplots(1,3,figsize=(13,4.7))
    for key,caption,color in zip(labels,captions,colors):
        for ax,metric,ylabel in zip(axs,['max_error_m','rms_error_m','final_max_error_m'],['観測区間内の最大誤差 [m]','観測区間全体のRMS [m]','終了時最大誤差 [m]']):
            ax.plot([7,30,60],[comparisons[d][key]['formation'][metric] for d in [7,30,60]],'o-',label=caption,color=color)
            ax.set_xlabel('観測終了日');ax.set_ylabel(ylabel);ax.set_xticks([7,30,60]);ax.grid(alpha=.2);ax.legend(fontsize=9)
    fig.suptitle('7日で見えた違いは30日・60日でも維持されるか')
    fig.tight_layout();save(fig,'02_観測期間による性能差')

    fig,axs=plt.subplots(2,1,figsize=(12,8),sharex=True)
    tt=np.arange(0,60*86400+1,3600,dtype=float)
    for key,caption,color in zip(labels,captions,colors):
        events=read_csv(src/f'{key}_events.csv');starts=np.sort(col(events,'start_time'))
        axs[0].step(tt/86400,np.searchsorted(starts,tt,side='right'),where='post',label=caption,color=color)
        exposure=np.zeros_like(tt)
        for event in events:
            start=float(event['start_time']);duration=float(event['complete_time'])-start
            x=np.maximum(0,tt-start);z=np.minimum(x/duration,1.)
            integral=duration*(z**3-.5*z**4)+np.maximum(0,x-duration)
            exposure+=(1 if event['to_mode']=='H' else -1)*integral
        assert np.isclose(exposure[-1],summary[key]['formation']['high_equivalent_total_s'],rtol=1e-9,atol=1e-4)
        axs[1].plot(tt/86400,exposure/3600,label=caption,color=color)
    axs[0].set_ylabel('累積姿勢遷移 [回]');axs[1].set_ylabel('累積HIGH換算量 [衛星・時間]');axs[1].set_xlabel('経過時間 [日]')
    for ax in axs:ax.grid(alpha=.2);ax.legend()
    fig.suptitle('切替とHIGH面積曝露の長期負荷')
    fig.tight_layout();save(fig,'03_累積制御負荷')

    fig,axs=plt.subplots(1,2,figsize=(12,5))
    for key,caption,color in zip(labels,captions,colors):
        rows=history[key];t=col(rows,'time_s')/86400
        axs[0].plot(t,col(rows,'extra_mean_altitude_loss_m'),label=caption,color=color)
        segments=summary[key]['formation']['segments']
        axs[1].plot(range(3),[x['extra_mean_altitude_loss_per_day_m'] for x in segments],'o-',label=caption,color=color)
    axs[0].set_xlabel('経過時間 [日]');axs[0].set_ylabel('全LOW比の追加平均高度損失 [m]')
    axs[1].set_xticks([0,1,2],['0–7日','7–30日','30–60日']);axs[1].set_ylabel('区間平均の追加高度損失 [m/日]');axs[1].set_xlabel('評価区間')
    for ax in axs:ax.grid(alpha=.2);ax.legend()
    fig.suptitle('同じprofileの全LOWを基準にした高度コスト')
    fig.tight_layout();save(fig,'04_高度コストと区間別増分')

    fig,axs=plt.subplots(1,2,figsize=(12,5))
    for key,caption,color in zip(labels,captions,colors):
        rows=read_csv(src/f'{key}_cycles.csv')
        for ax,field,scale in [(axs[0],'V_incident_after_m_s',1000),(axs[1],'next_wait_s',1/3600)]:
            x=np.sort(col(rows,field))*scale
            ax.step(x,np.arange(1,len(x)+1)/len(x),where='post',label=caption,color=color)
    axs[0].set_xlabel('LOW復帰時の接続最大速度 [mm/s]');axs[1].set_xlabel('同機のLOW完了→次HIGH [時間]');axs[1].set_xscale('log')
    for ax in axs:ax.set_ylabel('累積割合');ax.grid(alpha=.2);ax.legend()
    fig.suptitle('60日間の残留速度と再切替待ち時間')
    fig.text(.5,-.02,'待ち時間は次回HIGHを観測した区間のみ。末尾打切りの件数と観測時間は集計へ保存。',ha='center',fontsize=10)
    fig.tight_layout();save(fig,'05_残留速度と再切替')


if __name__=='__main__':main()
