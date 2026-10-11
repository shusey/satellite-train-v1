"""Five Japanese screening figures from saved physics only."""
import json
import numpy as np
import matplotlib.pyplot as plt
from plot_nrlmsis_validation import japanese_style
from run_dispersion_screening import ROOT,OUT,DATA,DOC,configuration
from sat_string.controller import gap_kinematics
from sat_string.passive_dispersion import theory_initial_gap_rate

def main():
    japanese_style();figures=DOC/'figures';figures.mkdir(exist_ok=True)
    if (figures/'01_initial_a.png').exists():raise FileExistsError('Published figures already exist')
    groups=json.loads((DATA/'group_summary.json').read_text())['groups']
    metrics=json.loads((DATA/'metrics.json').read_text())['results'];spec,cfg=configuration(7)
    def save(fig,name):
        fig.tight_layout();fig.savefig(figures/name,dpi=150,bbox_inches='tight');plt.close(fig)
    for factor,name,title,label,transform in [('a','01_initial_a.png','初期半長軸分散：全LOW・乱数5通り','指定した半長軸標準偏差 [m]',1),
            ('B','02_B_dispersion.png','抗力面積・質量比の分散：全LOW・乱数5通り','指定した C_D A/m の相対標準偏差 [%]',100)]:
        selected=[g for g in groups if g['factor']==factor];x=np.arange(len(selected));fig,ax=plt.subplots(figsize=(9,5))
        for days,style,color in [('7','o','C0'),('30','s','C1')]:
            valid=[(i,g['horizons'][days]) for i,g in enumerate(selected) if days in g['horizons']]
            if not valid:continue
            xx=np.array([i for i,h in valid]);yy=np.array([h['median_max_error_m'] for i,h in valid])
            lo=np.array([h['min_max_error_m'] for i,h in valid]);hi=np.array([h['max_max_error_m'] for i,h in valid])
            ax.errorbar(xx,yy,yerr=[yy-lo,hi-yy],fmt=style,color=color,capsize=5,
                label=f'{days}日：中央値・範囲（30日は延長した乱数のみ）')
        stress=next(r for r in metrics if r['factor']==factor and r['kind']=='stress' and r['days']==7)
        i=next(i for i,g in enumerate(selected) if g['sigma']==stress['nominal_sigma'])
        ax.plot(i+.08,stress['max_error_m'],'D',color='k',label='交互配置：負荷確認（分布・上限ではない）')
        baseline=next(r for r in metrics if r['factor']=='zero' and r['days']==7)
        ax.axhline(baseline['max_error_m'],ls=':',color='0.4',label='分散なし：7日')
        ax.set(xticks=x,xticklabels=[f"{g['sigma']*transform:g}" for g in selected],xlabel=label,ylabel='最大絶対間隔誤差 [m]',title=title,yscale='log')
        ax.legend(fontsize=8);save(fig,name)
    directory=OUT/'initial_7d';meta=json.loads((directory/'metadata.json').read_text());profiles=meta['profiles']
    t=np.load(directory/'time_s.npy');ix=np.arange(0,len(t),360);days=t[ix]/86400
    a=np.load(directory/'a_m.npy',mmap_mode='r');u=np.load(directory/'lambda_rad.npy',mmap_mode='r')
    selected_curves={}
    for factor,sigma in [('a',.1),('B',.003)]:
        g=next(g for g in groups if g['factor']==factor and g['sigma']==sigma)
        ident=g['horizons']['7']['representative_case_id'];k=next(i for i,p in enumerate(profiles) if p['id']==ident)
        e,v,_=gap_kinematics(a[ix,k],u[ix,k],reference_radius_m=cfg.derived.reference_radius_m,
            target_spacing_rad=cfg.derived.target_spacing_rad,mu_m3_s2=cfg.constants.mu_m3_s2)
        selected_curves[factor]=(e,v,profiles[k])
    e,v,p=selected_curves['a'];j=int(np.argmax(np.max(np.abs(e),axis=0)));fig,axs=plt.subplots(2,1,figsize=(9,7),sharex=True)
    axs[0].plot(days,e,color='0.75',lw=.5);axs[0].plot(days,e[:,j],label=f'最大誤差を示した間隔{j+1}')
    rate=theory_initial_gap_rate(cfg,p['delta_a_m']);axs[0].plot(days,rate[j]*t[ix],'k--',label='初期角速度差だけの一次予測')
    axs[0].set(ylabel='符号付き間隔誤差 [m]',title=f"初期半長軸0.1 m：代表乱数番号 {p['seed']}（全LOW）");axs[0].legend()
    axs[1].plot(days,v*1000,color='0.75',lw=.5);axs[1].plot(days,v[:,j]*1000,label=f'間隔{j+1}')
    axs[1].set(xlabel='開始から [日]',ylabel='間隔の変化速度 [mm/s]');axs[1].legend();save(fig,'03_initial_drift.png')
    hourly=np.genfromtxt(DATA/'representative_B_hourly.csv',delimiter=',',names=True)
    e,v,p=selected_curves['B'];j=int(np.argmax(np.max(np.abs(e),axis=0)));fig,axs=plt.subplots(3,1,figsize=(9,8),sharex=True)
    for gap in dict.fromkeys((0,9,j)):
        selected=hourly[hourly['gap']==gap+1];day=selected['time_s']/86400
        axs[0].plot(day,selected['excess_radius_gap_m'],label=f'間隔{gap+1}：ΔB/B0={np.diff(p["epsilon"])[gap]*100:.2f}%')
        axs[1].plot(day,selected['gap_rate_m_s']*1000)
        axs[2].plot(day,selected['gap_error_m'])
    axs[0].set(ylabel='分散なしとの差：隣接半長軸差 [m]',title=f"C_D A/m 分散0.3%：代表乱数番号 {p['seed']}（全LOW）");axs[0].legend(fontsize=8)
    axs[1].set(ylabel='間隔の変化速度 [mm/s]');axs[2].set(ylabel='間隔誤差 [m]',xlabel='開始から [日]')
    save(fig,'04_B_mechanism.png')
    fig,ax=plt.subplots(figsize=(9,5))
    with np.load(ROOT/spec['baseline_run']/'results.npz') as d:
        bt=d['time_s'];use=(bt>0)&(bt<=7*86400)&(bt%3600==0)
        be=np.max(np.abs(d['gap_error_m'][use]),axis=1);ax.plot(bt[use]/86400,be,label='分散なし・NRLMSISのみ')
    for f,label in [('a','初期半長軸0.1 m・代表乱数'),('B','C_D A/m 0.3%・代表乱数')]:
        ax.plot(days[1:],np.max(np.abs(selected_curves[f][0][1:]),axis=1),label=label)
    for threshold in (100,500,1000):ax.axhline(threshold,color='0.7',ls='--',lw=.6);ax.text(6.9,threshold,f'{threshold} m',ha='right',va='bottom',fontsize=8)
    ax.set(xlabel='開始から [日]',ylabel='最大絶対間隔誤差 [m]',yscale='log',title='同じ7日・同じ宇宙天気：小さい指定分散との比較')
    ax.legend(fontsize=9);save(fig,'05_baseline_comparison.png')

if __name__=='__main__':main()
