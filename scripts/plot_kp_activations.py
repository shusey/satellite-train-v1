"""Physical position errors at HIGH triggers, using saved activation CSV only."""
import csv
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from sat_string.long_horizon import decompose_q


def plot_activations(manifest,out):
    fig,axs=plt.subplots(2,3,figsize=(13,8),sharex=True)
    fields=[('command_max_abs_error_m','全列の最大誤差'),
        ('command_incident_max_abs_error_m','接続間隔の最大誤差'),
        ('command_rms_error_m','全列の瞬間RMS')]
    for i,entry in enumerate(manifest['cases']):
        with (Path(entry['analysis'])/'threshold_events.csv').open(encoding='utf-8-sig') as f:rows=list(csv.DictReader(f))
        for row,predicate,title in [(0,lambda x:True,'全HIGH指令'),
                (1,lambda x:x['quiet_packet_onset']=='True','全LOW後の活動開始')]:
            selected=[x for x in rows if predicate(x)]
            for col,(key,label) in enumerate(fields):
                values=[float(x[key]) for x in selected];ax=axs[row,col]
                if values:ax.boxplot([values],positions=[i],widths=.45,showfliers=False,
                    medianprops=dict(color='black'),manage_ticks=False)
                ax.set(title=title+'\n'+label,ylabel='間隔誤差 [m]')
                ax.text(i,.98,f'{len(values)}指令',transform=ax.get_xaxis_transform(),ha='center',va='top',fontsize=8,
                    bbox=dict(facecolor='white',alpha=.8,edgecolor='none',pad=1))
    for ax in axs.flat:
        ax.set_xticks(range(len(manifest['cases'])),[f"{x['kp_factor']:g}" for x in manifest['cases']])
        ax.set_xlabel('Kpの基準倍率');ax.grid(alpha=.2)
    fig.suptitle('HIGH起動時の位置誤差：中央値・四分位範囲・ひげ（外れ値は省略）')
    fig.tight_layout();fig.savefig(out/'activation_errors.png',dpi=150,bbox_inches='tight');plt.close(fig)


def plot_quiet_growth(cases,out):
    """First full response end -> next response; avoid nearby propagation gaps."""
    fig,axs=plt.subplots(len(cases),1,figsize=(12,9))
    for ax,(factor,a,r,events,cfg,states,first) in zip(np.atleast_1d(axs),cases):
        candidates=[x for x in r['group_cycles'][:1] if x['actual_next_wait_s'] is not None
            and x['all_low_at_sample'] and x['predicted_satellite_id'] is not None]
        if not candidates:
            ax.set_title(f'Kp {factor:g}倍：第1群終了後の観測済み静穏区間なし');continue
        row=candidates[0];j=row['predicted_satellite_id']-1
        ix=np.flatnonzero((a['time_s']>=row['state_sample_time_s'])&(a['time_s']<=row['actual_next_command_s']))
        stride=max(1,len(ix)//5000);ix=ix[::stride]
        p,d=decompose_q(a['gap_error_m'][ix],a['gap_rate_m_s'][ix],cfg.derived.controller_kp_s2_inv,cfg.derived.controller_kd_s_inv)
        x=(a['time_s'][ix]-row['state_sample_time_s'])/86400
        for y,label,color in [(p[:,j],'位置項 q_P','tab:blue'),(d[:,j],'速度項 q_D','tab:orange'),
                (p[:,j]+d[:,j],'合計 q','black')]:ax.plot(x,y/1e-7,label=label,color=color)
        ax.axhline(cfg.controller.threshold_on_m_s2/1e-7,color='red',ls='--',label='HIGH閾値 h_on')
        ax.axhline(cfg.controller.threshold_off_m_s2/1e-7,color='green',ls=':',label='LOW閾値 h_off')
        ax.set(title=f"Kp {factor:g}倍・予測衛星{j+1}：全LOWの{row['state_sample_time_s']/86400:.3f}日から次回指令まで",
            xlabel='全LOW境界からの経過時間 [日]',ylabel='制御判定量 [10^-7 m/s²]')
        ax.legend(ncol=5,fontsize=8)
    fig.suptitle('位置誤差の蓄積：第1応答群終了から次回切替（条件ごとに異なる終了状態）')
    fig.tight_layout();fig.savefig(out/'quiet_q_growth.png',dpi=150,bbox_inches='tight');plt.close(fig)
