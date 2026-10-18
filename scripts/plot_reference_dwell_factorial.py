"""Six Japanese-labelled scientific figures from saved factorial analyses."""
import os
from pathlib import Path
from run_reference_dwell_factorial import OUTPUT
from sat_string.dwell_factorial import CASE_ORDER
os.environ.setdefault('MPLCONFIGDIR',str(OUTPUT/'matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import json
import csv


def main():
    dest=OUTPUT/'figures';dest.mkdir(exist_ok=True) # Regenerate this experiment's figures only; never touch raw runs.
    plt.rcParams.update({'font.family':'Yu Gothic','axes.unicode_minus':False,'font.size':11})
    colors=['#4c566a','#0072b2','#d55e00','#009e73']
    report=json.loads((OUTPUT/'analysis/comparison.json').read_text())
    stats=[report['cases'][k] for k in CASE_ORDER]
    histories=[]
    for label in CASE_ORDER:
        with (OUTPUT/'analysis'/label/f'{label}_history_60s.csv').open(encoding='utf-8-sig') as f:
            histories.append({k:np.array(v) for k,v in _columns(csv.DictReader(f)).items()})
    def save(fig,name):
        fig.tight_layout();fig.savefig(dest/name,dpi=150);plt.close(fig)
    for key,title,ylabel,name in [
        ('max_error_m','最大間隔誤差：HIGH側・LOW側待機時間の分離','最大絶対間隔誤差 [m]','01_max_error.png'),
        ('instantaneous_RMS_m','瞬間RMS：全20間隔の誤差','瞬間RMS間隔誤差 [m]','02_RMS.png')]:
        fig,ax=plt.subplots(figsize=(10,4.6))
        for label,color,h in zip(CASE_ORDER,colors,histories):ax.plot(h['time_s']/86400,h[key],label=label,color=color,lw=1)
        ax.axvline(7,color='grey',ls=':',lw=1);ax.set(xlabel='経過時間 [日]',ylabel=ylabel,title=title,xlim=(0,30))
        ax.grid(alpha=.25);ax.legend(ncol=1);save(fig,name)
    fig,axs=plt.subplots(1,2,figsize=(11,4.8))
    for ax,key,title in [(axs[0],'residual_abs_incident_gap_velocity_m_s','接続間隔を個別に集計'),
                         (axs[1],'residual_max_incident_velocity_per_cycle_m_s','各cycleの接続間隔最大値')]:
        for i,(s,color) in enumerate(zip(stats,colors)):
            d=s[key];ax.vlines(i,d['p25']*1000,d['p75']*1000,color=color,lw=9)
            ax.scatter(i,d['median']*1000,c=color,edgecolor='black',zorder=3)
            ax.scatter(i,d['p90']*1000,c=color,marker='_',s=130)
        ax.set_xticks(range(4),CASE_ORDER,rotation=20,ha='right')
        ax.set(title=title,ylabel='LOW復帰時の相対速度絶対値 [mm/s]');ax.grid(axis='y',alpha=.25)
    fig.suptitle('残留速度の分布：点＝中央値、太線＝25–75%、短線＝90%')
    save(fig,'03_residual_velocity.png')
    fig,axs=plt.subplots(1,2,figsize=(11,4.8))
    for ax,values,ylabel in [(axs[0],[s['completed_cycles'] for s in stats],'完了L-H-L回数'),
                             (axs[1],[s['formation']['high_equivalent_total_s']/3600 for s in stats],'HIGH換算量 [衛星・時間]')]:
        bars=ax.bar(range(4),values,color=colors)
        ax.bar_label(bars,fmt='%.1f',fontsize=10,padding=3)
        ax.set_xticks(range(4),CASE_ORDER,rotation=20,ha='right');ax.set(ylabel=ylabel);ax.margins(y=.16)
    fig.suptitle('30日間の操作回数と抗力投入量');save(fig,'04_activity.png')
    fig,ax=plt.subplots(figsize=(10,4.6))
    for label,color,h in zip(CASE_ORDER,colors,histories):
        ax.plot(h['time_s']/86400,h['extra_mean_altitude_loss_m'],color=color,label=label)
    ax.set(title='同じA-ref全LOWからの追加平均高度損失',xlabel='経過時間 [日]',ylabel='追加平均高度損失 [m]',xlim=(0,30))
    ax.legend(ncol=1);ax.grid(alpha=.25);save(fig,'05_altitude_cost.png')
    names=['最大誤差','RMS','完了cycle','HIGH換算量','追加高度損失','残留速度中央値']
    keys=['max_error_m','RMS_m','cycles','HIGH_equivalent_satellite_hours','extra_altitude_loss_m','residual_cycle_max_velocity_median_mm_s']
    baseline=[stats[0]['formation']['max_error_m'],stats[0]['formation']['rms_error_m'],stats[0]['completed_cycles'],
              stats[0]['formation']['high_equivalent_total_s']/3600,stats[0]['formation']['extra_mean_altitude_loss_m'],
              stats[0]['residual_max_incident_velocity_per_cycle_m_s']['median']*1000]
    lookup={r['metric']:r for r in report['contrasts']}
    fig,ax=plt.subplots(figsize=(10,4.8));y=np.arange(6)
    for offset,key,label,color in [(-.24,'HIGH','HIGH側の短縮','#0072b2'),(0,'LOW','LOW側の短縮','#d55e00'),(.24,'interaction','相互作用','#009e73')]:
        ax.barh(y+offset,[100*lookup[k][key]/b for k,b in zip(keys,baseline)],height=.22,label=label,color=color)
    ax.set_yticks(y,names);ax.invert_yaxis();ax.axvline(0,color='black',lw=.8)
    ax.set(title='記述的2×2対比：統計的な有意差検定ではない',xlabel='H300-L300の値で規格化した対比 [%]')
    ax.legend();ax.grid(axis='x',alpha=.25);save(fig,'06_contrasts.png')


def _columns(rows):
    result={}
    for row in rows:
        for key,value in row.items():result.setdefault(key,[]).append(float(value))
    return result


if __name__=='__main__':main()
