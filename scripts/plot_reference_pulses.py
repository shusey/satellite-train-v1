"""Five Japanese research figures from offline pulse attribution tables."""
from pathlib import Path
import argparse,csv,json,os
os.environ.setdefault('MPLCONFIGDIR',str(Path(__file__).resolve().parents[1]/'outputs/mpl_cache'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import LogNorm

ROOT=Path(__file__).resolve().parents[1]


def read(path):
    with path.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))


def col(rows,key):return np.array([float(x[key]) for x in rows])


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',default='outputs/reference_pulse_mechanism_20261017_v3')
    ap.add_argument('--output',default='docs/reference_pulse_mechanism/figures_final');args=ap.parse_args()
    src=ROOT/args.input;dst=ROOT/args.output;dst.mkdir(parents=True,exist_ok=False)
    names={f.name for f in font_manager.fontManager.ttflist}
    for f in ['Yu Gothic','Meiryo','Noto Sans CJK JP']:
        if f in names:plt.rcParams['font.family']=f;break
    plt.rcParams.update({'font.size':11,'axes.unicode_minus':False,'axes.spines.top':False,'axes.spines.right':False})
    colors=['#167694','#ce6f23','#7853a5'];cases=['A','B','AB']
    labels=['A：初期半径差','B：抗力個体差','AB：両方']
    def save(fig,name):
        fig.savefig(dst/f'{name}.png',dpi=180,bbox_inches='tight');plt.close(fig)
    local={c:read(src/f'{c}_pulse_cycles.csv') for c in cases}
    fig,axs=plt.subplots(1,3,figsize=(14,5),sharey=True)
    for ax,c,label,color in zip(axs,cases,labels,colors):
        vals=[1000*col(local[c],key) for key in ['own_excess_g_impulse_m_s','V_incident_before_m_s','V_incident_after_m_s']]
        ax.boxplot(vals,showfliers=False,tick_labels=['自身のパルス\n速度寄与','操作前の\n隣接最大速度','操作後の\n隣接最大速度'],patch_artist=True,
                   boxprops={'facecolor':color,'alpha':.35},medianprops={'color':'black'})
        ax.set_title(label);ax.grid(axis='y',alpha=.2)
    axs[0].set_ylabel('速度量 [mm/s]')
    fig.suptitle('1回のHIGH操作による速度寄与と、操作前後の残留速度')
    fig.text(.5,-.02,'箱：25–75百分位、中央線：中央値、ひげ：1.5×四分位範囲。自身の寄与は近傍操作との相殺前。全完了サイクルを集計。',ha='center',fontsize=10)
    fig.tight_layout();save(fig,'01_パルス量と残留速度')

    rows=read(src/'A_first_pulse_trace.csv');t=col(rows,'time_s');qp=col(rows,'qP');qd=col(rows,'qD')
    fig,axs=plt.subplots(3,1,figsize=(12,8),sharex=True)
    axs[0].plot(t,qp/1e-7,label='位置由来 P');axs[0].plot(t,qd/1e-7,label='速度由来 D');axs[0].plot(t,(qp+qd)/1e-7,'k',label='判定量 q')
    axs[0].axhline(1.5,color='#c55',ls='--',label='HIGH開始閾値');axs[0].axhline(1.4,color='#886',ls=':',label='LOW帰還閾値')
    axs[0].set_ylabel(r'制御判定量 [$10^{-7}$ m/s²]');axs[0].legend(ncol=3,fontsize=9)
    axs[1].plot(t,col(rows,'area_m2'),color='#d78b24');axs[1].set_ylabel('衛星3の面積 [m²]')
    for key in rows[0]:
        if key.endswith('_v_m_s'):axs[2].plot(t,col(rows,key)*1000,label=f'間隔{key.split("_")[0][3:]}')
    axs[2].axhline(0,color='grey',lw=.7);axs[2].set_ylabel('相対速度 [mm/s]');axs[2].legend()
    for ax in axs:
        ax.axvspan(120,420,color='#f4ce76',alpha=.2);ax.axvline(540,color='grey',ls='--');ax.axvline(110,color='#c55',ls=':',alpha=.6);ax.grid(alpha=.2)
    axs[-1].set_xlabel('シミュレーション開始からの時間 [s]')
    fig.suptitle('A-ref：最初の衛星3の応答（全体の制御開始時）')
    fig.text(.5,.01,'110 s：初めてLOW閾値以下 ／ 着色120–420 s：HIGH完了後の保持 ／ 540 s：LOW完了 ／ 930 s：同機が再起動',ha='center',fontsize=10)
    fig.tight_layout(rect=(0,.03,1,.96));save(fig,'02_最初のパルスと相対速度')

    examples=json.loads((src/'examples.json').read_text(encoding='utf-8'))
    fig,axs=plt.subplots(2,2,figsize=(13,7))
    for column,c in enumerate(['A','B']):
        ex=next(x for x in examples if x['case']==c and x['selection']=='first_quiet')
        rows=read(src/f'{c}_first_quiet_trace.csv');tt=col(rows,'time_s');sel=(tt>=ex['start_s'])&(tt<=ex['end_s'])
        tau=tt[sel]-ex['start_s'];qp=col(rows,'qP')[sel];qd=col(rows,'qD')[sel];pred=qp[0]+qd[0]+qd[0]/43200*tau
        ax=axs[0,column];ax.plot(tau/3600,(qp+qd)/1e-7,label='保存結果の q');ax.plot(tau/3600,pred/1e-7,'--',label='D一定の予測')
        ax.axhline(1.5,color='#b55',ls=':',label='HIGH開始閾値');ax.set_title(f'{c}-ref：最初の全LOW区間・次に起動する衛星{ex["trigger_satellite"]}');ax.legend(fontsize=9)
        ax=axs[1,column];ax.plot(tau/3600,(qp-qp[0])/1e-7,label='位置由来 P の変化');ax.plot(tau/3600,(qd-qd[0])/1e-7,label='速度由来 D の変化');ax.legend(fontsize=9);ax.set_xlabel('全衛星LOWへ戻ってからの時間 [h]')
    for ax in axs.flat:ax.set_ylabel(r'制御判定量 [$10^{-7}$ m/s²]');ax.grid(alpha=.2)
    fig.suptitle('次の切替：Aでは残留速度、Bでは速度の再成長も寄与')
    fig.tight_layout();save(fig,'03_全LOW期間と再切替')

    fig,axs=plt.subplots(1,3,figsize=(14,5),sharex=True,sharey=True);hexbins=[]
    for ax,c,label in zip(axs,cases,labels):
        rows=read(src/f'{c}_incident_gaps.csv');x=col(rows,'v_before_m_s')*1000;y=col(rows,'v_after_m_s')*1000
        im=ax.hexbin(x,y,gridsize=45,mincnt=1,bins='log',extent=(-7,7,-7,7),cmap='viridis')
        hexbins.append(im)
        ax.plot([-7,7],[-7,7],ls='--',color='gray');ax.axvline(0,color='gray',lw=.6);ax.axhline(0,color='gray',lw=.6)
        ax.set_title(label);ax.set_xlabel('パルス開始時の相対速度 [mm/s]');ax.set_aspect('equal')
    norm=LogNorm(vmin=1,vmax=max(float(x.get_array().max()) for x in hexbins))
    for im in hexbins:im.set_norm(norm)
    axs[0].set_ylabel('LOW完了時の相対速度 [mm/s]');fig.colorbar(im,ax=axs.tolist(),label='記録件数（対数色尺度）',shrink=.7)
    fig.suptitle('速度の符号反転は多いが、位置誤差の補正に必要な反転も含む')
    fig.text(.5,.02,'各局所サイクルの隣接間隔を集計。同じ間隔が重複登場するため、独立な試行の確率とは解釈しない。',ha='center',fontsize=10)
    save(fig,'04_符号付き速度の前後比較')

    seg=read(src/'interval_summary.csv');fig,axs=plt.subplots(1,2,figsize=(12,5))
    for c,label,color in zip(cases,labels,colors):
        rows=[x for x in seg if x['case']==c];xx=np.arange(3)
        axs[0].plot(xx,1000*col(rows,'median_own_excess_g_impulse_m_s'),'o-',label=label,color=color)
        axs[1].plot(xx,col(rows,'cycles')/(col(rows,'end_day')-col(rows,'start_day')),'o-',label=label,color=color)
    for ax in axs:ax.set_xticks([0,1,2],['0–7日','7–30日','30–60日']);ax.set_xlabel('観測区間');ax.grid(alpha=.2);ax.legend()
    axs[0].set_ylabel('1パルスの速度寄与・中央値 [mm/s]');axs[1].set_ylabel('完了サイクルの開始数 / 日')
    fig.suptitle('1回の補正量と切替頻度は後半に増加（全完了サイクル）')
    fig.tight_layout();save(fig,'05_期間別パルス量と頻度')


if __name__=='__main__':main()
