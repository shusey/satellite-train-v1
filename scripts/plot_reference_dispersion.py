"""Six Japanese figures from saved postprocessing caches; no physical runs."""
import argparse,json,os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR',str(Path(__file__).resolve().parents[1]/'outputs/reference_dispersion_20261012/matplotlib'))
import numpy as np
import matplotlib.pyplot as plt
import yaml
from run_reference_dispersion import ROOT,SPEC,DOC,DATA
from plot_nrlmsis_validation import japanese_style


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--comparison-only',action='store_true')
    args=parser.parse_args()
    spec=yaml.safe_load(SPEC.read_text());out=ROOT/spec['output'];days=60
    dest=DOC/'figures';dest.mkdir(exist_ok=args.comparison_only);japanese_style()
    cases=['A','B','AB'];labels={'A':'A：初期半長軸分散のみ','B':'B：持続的な抗力個体差のみ','AB':'AB：両方'}
    arrays={c:dict(np.load(out/f'{c}_plot_{days}d.npz')) for c in cases}
    def save(fig,name):
        if (dest/name).exists():raise FileExistsError(dest/name)
        fig.tight_layout();fig.savefig(dest/name,dpi=160,bbox_inches='tight');plt.close(fig)
    fig,axes=plt.subplots(3,1,figsize=(10,10),sharex=True)
    for ax,c in zip(axes,cases):
        a=arrays[c];t=a['time_s']/86400
        off=json.loads((DATA/f'60d/{c}_OFF_metrics.json').read_text())
        bad=off['first_label_order_reversal_s'];use=np.ones(len(t),bool) if bad is None else t<bad/86400
        ax.plot(t[use],np.abs(a['off_error'][use]).max(axis=1),'--',color='0.4',label='全LOW：順序逆転前')
        if bad is not None:ax.axvline(bad/86400,color='0.4',ls=':',label=f'全LOW順序逆転 {bad/86400:.2f}日')
        ax.plot(t,np.abs(a['on_error']).max(axis=1),label='現行制御')
        ax.axhline(500,color='0.6',ls=':',label='Reference評価閾値500 m')
        ax.set(title=labels[c],ylabel='最大絶対間隔誤差 [m]',yscale='symlog')
        ax.set_yscale('symlog',linthresh=500);ax.legend(fontsize=8);ax.grid(alpha=.2)
        ax.set_ylim(bottom=0)
        ax.yaxis.set_label_coords(-.065,.5)
    axes[-1].set_xlabel('開始からの時間 [日]');save(fig,'01_error_comparison.png')
    if args.comparison_only:return
    for number,c in enumerate(cases,2):
        a=arrays[c];t=a['time_s']/86400;worst=int(np.argmax(np.abs(a['on_error']).max(axis=0)))
        fig,axes=plt.subplots(2,1,figsize=(10,7),sharex=True)
        axes[0].plot(t,a['on_error'],color='0.7',lw=.4)
        axes[0].plot(t,a['on_error'][:,worst],label=f'最大誤差の間隔{worst+1}')
        for sign in (-1,1):axes[0].axhline(sign*500,color='0.4',ls=':')
        axes[0].set(title=labels[c]+'：制御あり・全20間隔',ylabel='符号付き間隔誤差 [m]');axes[0].legend()
        axes[1].plot(t,a['on_rate']*1000,color='0.7',lw=.4)
        axes[1].plot(t,a['on_rate'][:,worst]*1000,label=f'間隔{worst+1}')
        axes[1].set(xlabel='開始からの時間 [日]',ylabel='相対間隔速度 [mm/s]');axes[1].legend()
        save(fig,f'{number:02}_{c}_gaps.png')
    checkpoint=json.loads((out/'batch_60d/checkpoint.json').read_text())
    fig,axes=plt.subplots(2,1,figsize=(10,7),sharex=True)
    for j,c in enumerate(cases):
        ev=checkpoint['events'][j];times=sorted(x['start_time']/86400 for x in ev if x['start_time']<=60*86400)
        axes[0].step([0]+times+[60],[0]+list(range(1,len(times)+1))+[len(times)],where='post',label=labels[c])
        a=arrays[c];axes[1].step(a['activity_time_s']/86400,a['high_count'],where='post',alpha=.7,label=labels[c])
    axes[0].set(title='Reference条件の累積制御負荷',ylabel='累積姿勢遷移 [回]');axes[0].legend()
    axes[1].set(xlabel='開始からの時間 [日]',ylabel='同時に安定HIGHの衛星 [機]');axes[1].legend()
    save(fig,'05_switching_activity.png')
    fig,ax=plt.subplots(figsize=(10,5))
    for c in cases:
        a=arrays[c];ax.plot(a['time_s']/86400,np.mean(a['off_a_m']-a['a_m'],axis=1),label=labels[c])
    ax.set(title='同一profileの全LOWに対する制御コスト',xlabel='開始からの時間 [日]',ylabel='追加平均高度損失 [m]')
    ax.legend();ax.grid(alpha=.2);save(fig,'06_extra_altitude_loss.png')


if __name__=='__main__':main()
