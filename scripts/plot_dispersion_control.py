"""Six Japanese research figures; saved arrays only, no simulation."""
import argparse,json
import numpy as np
import matplotlib.pyplot as plt
from run_dispersion_control import OUT,DATA,DOC,configuration
from plot_nrlmsis_validation import japanese_style


def first_trigger_figure(cfg,save,days):
    labels={'A':'A：初期半長軸残差','B':'B：持続的な抗力個体差'}
    fig,axs=plt.subplots(2,2,figsize=(12,7))
    for col,c in enumerate(['A','B']):
        a=dict(np.load(OUT/f'{c}_first_trigger_{days}d.npz'));mins=(a['time_s']-a['command_time_s'])/60
        ax=axs[0,col]
        for key,label in [('qP','位置誤差由来のP項'),('qD','速度差由来のD項'),('q','合計q')]:ax.plot(mins,a[key]/1e-7,label=label)
        ax.axhline(cfg.controller.threshold_on_m_s2/1e-7,ls='--',color='k',label='HIGH起動閾値')
        ax.axhline(cfg.controller.threshold_off_m_s2/1e-7,ls=':',color='k',label='LOW帰還閾値')
        ax.axvline(0,color='0.6',ls=':');ax.set(ylabel='判定量 [10^-7 m/s²]',title=f'{labels[c]}：最初のHIGH指令（衛星{int(a["satellite_id"])}）');ax.legend(fontsize=8)
        j=int(a['satellite_id'])-1;incident=[g for g in (j-1,j) if 0<=g<20]
        for g in incident:axs[1,col].plot(mins,a['gap_rate_m_s'][:,g]*1000,label=f'間隔{g+1}')
        axs[1,col].set(xlabel='初回HIGH指令から [分]',ylabel='関連する間隔速度 [mm/s]');axs[1,col].legend()
    save(fig,'06_first_trigger_PD.png')


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--days',type=int,choices=[30,60],required=True)
    parser.add_argument('--first-trigger-only',action='store_true',help='Repair only figure6 in this study output')
    args=parser.parse_args();_,cfg=configuration(args.days);japanese_style()
    figures=DOC/'figures';figures.mkdir(exist_ok=args.first_trigger_only)
    cases=['A','B','AB'];labels={'A':'A：初期半長軸残差','B':'B：持続的な抗力個体差','AB':'AB：両方'}
    def save(fig,name):
        fig.tight_layout();fig.savefig(figures/name,dpi=160,bbox_inches='tight');plt.close(fig)
    if args.first_trigger_only:
        first_trigger_figure(cfg,save,args.days);return
    arrays={c:dict(np.load(OUT/f'{c}_plot_{args.days}d.npz')) for c in cases}
    fig,axes=plt.subplots(3,1,figsize=(10,9),sharex=True)
    for ax,c in zip(axes,cases):
        a=arrays[c];t=a['time_s']/86400
        off=json.loads((DATA/f'{args.days}d/{c}_OFF_metrics.json').read_text())
        bad=off['first_label_order_reversal_s'];before=np.ones(t.shape,dtype=bool) if bad is None else t<bad/86400
        envelope=np.max(np.abs(a['off_error']),axis=1)
        ax.plot(t[before],envelope[before],'--',color='0.4',label='全LOW（制御なし）')
        if bad is not None:
            ax.plot(t[~before],envelope[~before],':',color='0.4',label='順序逆転後：固定番号の診断値')
            ax.axvline(bad/86400,color='0.5',ls=':')
        ax.plot(t,np.max(np.abs(a['on_error']),axis=1),label='現行制御（再調整なし）')
        ax.axhline(100,color='0.6',ls=':',label='仮想許容値100 m')
        ax.set(ylabel='最大絶対間隔誤差 [m]',title=labels[c],yscale='symlog',ylim=(0,None))
        ax.set_yscale('symlog',linthresh=100);ax.legend(fontsize=8);ax.grid(alpha=.2)
    axes[-1].set_xlabel('開始からの時間 [日]');save(fig,'01_precision_comparison.png')
    for c,name in [('A','02_A_gaps.png'),('B','03_B_gaps.png')]:
        a=arrays[c];t=a['time_s']/86400
        fig,axs=plt.subplots(2,1,figsize=(10,7),sharex=True)
        worst=int(np.argmax(np.max(np.abs(a['on_error']),axis=0)))
        axs[0].plot(t,a['on_error'],color='0.7',lw=.5)
        axs[0].plot(t,a['on_error'][:,worst],label=f'最大誤差の間隔{worst+1}')
        axs[0].axhline(100,ls=':',color='0.4');axs[0].axhline(-100,ls=':',color='0.4')
        axs[0].set(ylabel='符号付き間隔誤差 [m]',title=labels[c]+'：制御あり・全20間隔')
        axs[0].legend()
        axs[1].plot(t,a['on_rate']*1000,color='0.7',lw=.5)
        axs[1].plot(t,a['on_rate'][:,worst]*1000,label=f'間隔{worst+1}')
        axs[1].set(xlabel='開始からの時間 [日]',ylabel='相対間隔速度 [mm/s]');axs[1].legend()
        save(fig,name)
    fig,axs=plt.subplots(2,1,figsize=(10,7),sharex=True)
    b=arrays['B'];t=b['activity_time_s']/86400
    axs[0].step(t,b['high_count'],where='post',label='安定HIGH')
    axs[0].step(t,b['nonLOW_count'],where='post',alpha=.65,label='姿勢変更中を含む非LOW')
    axs[0].set(ylabel='衛星数 [機]',title='持続的な抗力個体差への制御活動（全10秒保存の状態変化）');axs[0].legend()
    for c in cases:
        events=json.loads((OUT/('batch_60d' if args.days==60 else 'batch_30d')/'checkpoint.json').read_text())
        own=events['events'][cases.index(c)];times=sorted(z['start_time']/86400 for z in own)
        axs[1].step([0]+times+[args.days],[0]+list(range(1,len(times)+1))+[len(times)],where='post',label=labels[c])
    axs[1].set(xlabel='開始からの時間 [日]',ylabel='累積姿勢遷移数 [回]');axs[1].legend();save(fig,'04_control_activity.png')
    fig,ax=plt.subplots(figsize=(10,5))
    for c in cases:
        a=arrays[c];ax.plot(a['time_s']/86400,np.mean(a['off_a_m']-a['a_m'],axis=1),label=labels[c])
    ax.set(xlabel='開始からの時間 [日]',ylabel='同じ配置の全LOW比：追加平均高度損失 [m]',title='精度改善のために支払った追加高度損失')
    ax.legend();ax.grid(alpha=.2);save(fig,'05_extra_altitude_loss.png')
    first_trigger_figure(cfg,save,args.days)


if __name__=='__main__':main()
