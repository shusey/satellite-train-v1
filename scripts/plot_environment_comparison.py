"""Five Japanese figures at most; compact saved post-processing only."""
import json
from pathlib import Path
import numpy as np
from plot_nrlmsis_validation import japanese_style
import matplotlib.pyplot as plt
from analyze_environment_comparison import ROOT,OUT,DOC,csv_save

LABELS=dict(quiet='静穏側（1月1日開始）',active='活動期（3月1日開始）',strong='強嵐期（5月1日開始）')


def save(fig,name):
    fig.tight_layout();fig.savefig(DOC/'figures'/name,dpi=150,bbox_inches='tight');plt.close(fig)


def main(refresh=False):
    japanese_style();(DOC/'figures').mkdir(exist_ok=True)
    if (DOC/'figures'/'01_weather.png').exists() and not refresh:raise FileExistsError('Figures already exist')
    protocol=json.loads((OUT/'protocol.json').read_text())
    comparison=json.loads((DOC/'data'/'comparison.json').read_text())['results']
    omni=np.loadtxt(OUT/'raw'/'omni2_2024.dat')
    utc=np.array([np.datetime64(f'{int(r[0])}-01-01')+np.timedelta64(int(r[1])-1,'D')+np.timedelta64(int(r[2]),'h') for r in omni])
    rows=[]
    fig,axs=plt.subplots(3,1,figsize=(10,8),sharex=True)
    for group,selection in protocol['selected'].items():
        start=np.datetime64(selection['start']);end=np.datetime64(selection['end_exclusive']);use=(utc>=start)&(utc<end)
        day=(utc[use]-start)/np.timedelta64(1,'D');o=omni[use]
        weather=np.genfromtxt(OUT/f'weather/{group}/drivers.csv',delimiter=',',names=True,dtype=None,encoding='utf-8')
        f107=weather['f107_previous_day'][(day*8).astype(int)]
        axs[0].plot(day,o[:,38]/10,label=LABELS[group],lw=.8)
        axs[1].plot(day,o[:,40],lw=.8)
        axs[2].plot(day,f107,lw=.8)
        for i,d in enumerate(day):rows.append(dict(group=group,utc=str(utc[use][i]),kp=float(o[i,38]/10),dst_nT=float(o[i,40]),f107_previous_day_sfu=float(f107[i])))
    axs[0].set(ylabel='Kp指数',title='間隔誤差を見る前に選定した宇宙天気期間');axs[0].legend(fontsize=9)
    axs[1].set(ylabel='Dst [nT]');axs[1].axhline(-250,color='k',ls='--',lw=.7,label='本研究の強嵐基準');axs[1].legend()
    axs[2].set(ylabel='前日F10.7 [sfu]',xlabel='各期間の開始から [日]')
    save(fig,'01_weather.png');csv_save(DOC/'data'/'space_weather_hourly.csv',rows)
    fig,axs=plt.subplots(3,1,figsize=(10,9),sharex=True)
    for group in LABELS:
        with np.load(OUT/f'{group}_hourly.npz') as d:
            day=d['time_s']/86400
            axs[0].plot(day,d['max_error_m'],label=LABELS[group])
            axs[1].plot(day,d['rms_error_m'])
            axs[2].plot(day,np.max(np.abs(d['current_arc_error_m']),axis=1))
    axs[0].set(ylabel='最大絶対間隔誤差 [m]',title='全LOW・無制御：主要定義は固定基準半径')
    axs[0].axhline(100,color='k',ls='--',lw=.7,label='仮想許容100 m')
    axs[0].legend()
    axs[1].set(ylabel='瞬間RMS間隔誤差 [m]')
    axs[2].set(ylabel='現在半径の弧長誤差 [m]',xlabel='各期間の開始から [日]',title='補助定義：共通の軌道収縮も含む')
    save(fig,'02_spacing.png')
    with np.load(OUT/'strong_mechanism.npz') as d:
        hours=(d['time_s']-d['time_s'][0])/3600
        fig,axs=plt.subplots(4,1,figsize=(10,9),sharex=True)
        for gap in (0,9,19):
            label=f'間隔{gap+1}（衛星{gap+1}–{gap+2}）'
            axs[0].plot(hours,d['adjacent_density'][:,gap]*1e12,label=label)
            axs[1].plot(hours,d['adjacent_Jrho'][:,gap]*1e11)
            axs[2].plot(hours,d['gap_rate_m_s'][:,gap]*1e6)
            axs[3].plot(hours,d['gap_error_m'][:,gap])
        axs[0].set(ylabel='隣接密度差 [×10^-12 kg/m³]',title='強嵐期：最大3時間apの前後（10秒保存データ）');axs[0].legend(fontsize=9)
        axs[1].set(ylabel='累積密度差 [×10^-11 kg s/m³]')
        axs[2].set(ylabel='間隔の変化速度 [µm/s]')
        axs[3].set(ylabel='間隔誤差 [m]',xlabel='図の先頭時刻から [時間]')
        save(fig,'03_mechanism.png')
    fig,ax=plt.subplots(figsize=(9,5));width=.24;x=np.arange(3)
    for j,m in enumerate(comparison):
        ax.bar(x+(j-1)*width,[r['exceedance_fraction']*100 for r in m['virtual_tolerances']],width,label=LABELS[m['group']])
    ax.set(xticks=x,xticklabels=['±100 m','±500 m','±1 km'],ylabel='60日中の超過時間 [%]',xlabel='仮想的な許容間隔誤差',title='固定基準半径の間隔誤差：10秒標本での超過割合')
    ax.legend();save(fig,'04_tolerances.png')
    observation=DOC/'data'/'grace_integrated_summary.json'
    if observation.exists():
        values=json.loads(observation.read_text())['results'];x=np.arange(len(values))
        labels=['静穏側：10月25日','回復期擾乱：11月1日','強嵐：10月29日'];fig,axs=plt.subplots(1,2,figsize=(11,5))
        obs=np.array([r['simultaneous_obs']['integral_observed_support'] for r in values])*1e7
        model=np.array([r['simultaneous_model']['integral_observed_support'] for r in values])*1e7
        envelope=np.array([r['simultaneous_obs']['fully_correlated_error_envelope'] for r in values])*1e7
        axs[0].errorbar(x-.06,obs,yerr=envelope,fmt='o',capsize=5,label='観測由来・誤差包絡付き')
        axs[0].plot(x+.06,model,'s',label='同一地点のNRLMSIS 2.0')
        axs[0].set(xticks=x,xticklabels=labels,ylabel='A−B積算密度差 [×10^-7 kg s/m³]',title='欠測を除いた観測支持区間のみ');axs[0].legend(fontsize=9)
        axs[1].bar(x,[r['simultaneous_obs']['coverage_fraction']*100 for r in values])
        axs[1].set(xticks=x,xticklabels=labels,ylabel='有効区間の時間被覆率 [%]',title='1日全部の積算履歴ではない',ylim=(0,100))
        fig.suptitle('2003年のGRACE約220 km対・3度平均：100 kmの直接証拠ではありません')
        save(fig,'05_grace_support.png')


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--refresh',action='store_true',help='Revise only this study\'s published figures; no run data changed')
    main(parser.parse_args().refresh)
