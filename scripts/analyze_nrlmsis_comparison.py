#!/usr/bin/env python3
"""Matched 60-day NRLMSIS A/B metrics, mechanism checks and five Japanese figures."""
import argparse
import csv
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from sat_string.config import load_config
from sat_string.io import load_case_arrays,load_switch_events,load_metadata,prepare_output_directory,write_json
from sat_string.evaluation import compute_evaluation_metrics
from sat_string.necessity import exceedance_metrics,improvement
from sat_string.nrlmsis_analysis import density_statistics,density_to_gap_acceleration,pd_components,current_radius_arc_error
from sat_string.kd_analysis import threshold_events
from sat_string.residual import residual_analysis
from sat_string.control_cost import clean_json
from plot_nrlmsis_validation import japanese_style
import matplotlib.pyplot as plt


def csv_rows(path,rows):
    if not rows:return
    with path.open('w',encoding='utf-8-sig',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def integrate(y,t):
    return np.r_[np.zeros((1,y.shape[1])),np.cumsum((y[:-1]+y[1:])*.5*np.diff(t)[:,None],axis=0)]


def daily(y,t,reduction=np.mean):
    # Exclude the endpoint, which alone is not an extra day of observations.
    y=np.asarray(y);size=int(round(86400/(t[1]-t[0])))
    blocks=y[:-1].reshape(-1,size,*y.shape[1:])
    return (np.arange(len(blocks))+.5),reduction(blocks,axis=1)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--all-low',type=Path,required=True);p.add_argument('--controlled',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    out=prepare_output_directory(args.output);japanese_style()
    config=load_config(args.all_low/'normalized_config.yaml')
    controlled_config=load_config(args.controlled/'normalized_config.yaml')
    if config.to_dict()!=controlled_config.to_dict():raise ValueError('Physical/controller settings differ')
    frozen=load_config('tests/fixtures/regression_24h/config.yaml').to_dict()
    frozen['atmosphere']['density_model']='nrlmsis2'
    frozen['disturbance']['enabled']=False
    frozen['simulation']['duration_s']=5184000
    if config.to_dict()!=frozen:raise ValueError('Physical/control/numerical settings differ from frozen baseline')
    environments=[json.loads((x/'density_environment.json').read_text(encoding='utf-8')) for x in [args.all_low,args.controlled]]
    for key in ['environment','cache_seconds','weather_sha256','parameter_version','pymsis_version']:
        if environments[0][key]!=environments[1][key]:raise ValueError('Environmental input differs: '+key)
    metas=[load_metadata(x) for x in [args.all_low,args.controlled]]
    if [m['case']['control_enabled'] for m in metas]!=[False,True] or any(m['case']['disturbance_enabled'] for m in metas):
        raise ValueError('Expected artificial-pulse-free A/B case flags')
    low=load_case_arrays(args.all_low);controlled=load_case_arrays(args.controlled)
    t=low['time_s']
    if t[-1]!=5184000 or not np.array_equal(t,controlled['time_s']):raise ValueError('Matched 60-day time grids required')
    for key in ['a_m','lambda_rad','mode','area_m2']:
        np.testing.assert_array_equal(low[key][0],controlled[key][0])
    checks=dict(configuration_equal=True,environment_equal=True,weather_equal=True,initial_states_equal=True,
        time_grid_equal=True,only_case_flag_difference='control_enabled',artificial_pulse_enabled=False,
        physical_control_numerical_settings_equal_frozen_baseline=True,
        allowed_changes_vs_frozen_24h=['atmosphere.density_model','disturbance.enabled','simulation.duration_s'])
    write_json(out/'input_checks.json',checks)
    items=[('all_low','全LOW・制御なし',low,[]),('baseline','基準制御',controlled,load_switch_events(args.controlled))]
    comparison={};series=[];gaps=[];exceed=[];rows=[];hourly=[]
    for key,label,arrays,events in items:
        e=arrays['gap_error_m'];v=arrays['gap_rate_m_s'];rho=arrays['density_kg_m3'];a=arrays['a_m']
        maximum=np.abs(e).max(axis=1);rms=np.sqrt(np.mean(e**2,axis=1))
        metrics=compute_evaluation_metrics(arrays,config,events,low_baseline=low,evaluation_start_s=0)
        metrics['wave']['applicable_as_central_disturbance_front']=False
        metrics['wave']['nrlmsis_interpretation']='Legacy center-based wave fields retained; distributed ongoing density has no prescribed central disturbance origin'
        metrics['density']=density_statistics(rho)
        metrics['formation']['max_error_percent_of_spacing']=100*maximum.max()/config.formation.target_spacing_m
        metrics['formation']['rms_error_percent_of_spacing']=100*metrics['formation']['rms_error_m']/config.formation.target_spacing_m
        metrics['formation']['final_error_percent_of_spacing']=100*maximum[-1]/config.formation.target_spacing_m
        peak_t,peak_gap=np.unravel_index(np.abs(e).argmax(),e.shape)
        metrics['formation'].update(max_abs_error_time_s=float(t[peak_t]),max_abs_error_gap_id=int(peak_gap+1),
            max_abs_error_signed_value_m=float(e[peak_t,peak_gap]))
        metrics['formation']['prefix_horizons']={}
        for days in [1,7,14,30,60]:
            stop=int(np.searchsorted(t,days*86400,side='right'));window=e[:stop]
            metrics['formation']['prefix_horizons'][str(days)]=dict(max_abs_error_m=float(np.abs(window).max()),
                rms_error_m=float(np.sqrt(np.mean(window**2))),final_max_abs_error_m=float(np.abs(window[-1]).max()))
        arc=current_radius_arc_error(a,arrays['lambda_rad'],config.formation.target_spacing_m)
        arc_max=np.abs(arc).max(axis=1)
        metrics['formation']['current_radius_arc_diagnostic']=dict(max_abs_error_m=float(arc_max.max()),
            rms_error_m=float(np.sqrt(np.mean(arc**2))),final_max_abs_error_m=float(arc_max[-1]),
            definition='adjacent mean current radius * phase difference - fixed 100 km; near-circular arc proxy, not exact instantaneous real-space distance',
            hypothetical_tolerances=[exceedance_metrics(t,arc_max,z) for z in [100,500,1000]])
        del arc,arc_max
        metrics['hypothetical_tolerances']=[exceedance_metrics(t,maximum,z) for z in [100,500,1000]]
        ppart,dpart=pd_components(e,v,config.derived.controller_kp_s2_inv,config.derived.controller_kd_s_inv)
        q=arrays['q_m_s2']
        metrics['q']=dict(maximum_positive_m_s2=float(q.max()),maximum_absolute_m_s2=float(np.abs(q).max()),
            threshold_on_m_s2=config.controller.threshold_on_m_s2,
            maximum_P_absolute_m_s2=float(np.abs(ppart).max()),maximum_D_absolute_m_s2=float(np.abs(dpart).max()),
            gap_rate_max_m_s=float(np.abs(v).max()),gap_rate_final_max_m_s=float(np.abs(v[-1]).max()))
        acc=density_to_gap_acceleration(a,rho,arrays['area_m2'],config)
        vi=integrate(acc,t)+v[0];ei=integrate(v,t)+e[0]
        metrics['mechanism']=dict(max_gap_rate_reconstruction_difference_m_s=float(np.max(np.abs(vi-v))),
            max_gap_error_reconstruction_difference_m=float(np.max(np.abs(ei-e))),
            method='trapezoidal integration of exact saved gap acceleration and rate; sampled 10 s diagnostic, not another propagator',
            forcing='ongoing geographic/time-varying density; constant-D quiet forecast not assumed valid')
        metrics['threshold_events']=threshold_events(arrays,config,events)
        report,states=residual_analysis(arrays,config,events,groups=[],low_arrays=low)
        write_json(out/(key+'_cycles.json'),report)
        if states:np.savez_compressed(out/(key+'_cycle_states.npz'),**states)
        write_json(out/(key+'_metrics.json'),clean_json(metrics));comparison[key]=metrics
        for j in range(e.shape[1]):
            gaps.append(dict(case=key,gap_id=j+1,max_abs_error_m=float(np.abs(e[:,j]).max()),rms_error_m=float(np.sqrt(np.mean(e[:,j]**2))),
                final_error_m=float(e[-1,j]),final_rate_m_s=float(v[-1,j])))
        exceed.extend(dict(case=key,**record) for record in metrics['hypothetical_tolerances'])
        series.append((key,label,arrays,maximum,rms))
        rows.append(dict(case=key,max_abs_error_m=metrics['formation']['max_abs_error_m'],
            percent_of_100km=metrics['formation']['max_error_percent_of_spacing'],rms_error_m=metrics['formation']['rms_error_m'],
            final_max_abs_error_m=float(maximum[-1]),transitions=metrics['switching']['n_transition_starts'],
            completed_cycles=sum(metrics['switching']['lhl_cycles_per_satellite']),
            high_stable_total_s=metrics['switching']['high_stable_total_s'],
            high_equivalent_total_s=metrics['switching']['high_equivalent_total_s'],
            mean_initial_altitude_loss_m=metrics['altitude']['mean_altitude_loss_m'],
            extra_mean_loss_vs_nrl_all_low_m=metrics['altitude']['extra_loss_vs_all_low_m']))
        # Hourly summaries retain physical data needed for lightweight reproduction.
        for k in range(0,len(t),360):
            rec=dict(case=key,time_s=float(t[k]),max_error_m=float(maximum[k]),rms_error_m=float(rms[k]),
                mean_density_kg_m3=float(rho[k].mean()),fleet_density_spread_kg_m3=float(np.ptp(rho[k])),
                maximum_q_m_s2=float(q[k].max()),mean_altitude_m=float(arrays['altitude_m'][k].mean()))
            for j in [0,9,19]:rec.update({f'gap_{j+1}_error_m':float(e[k,j]),f'gap_{j+1}_rate_m_s':float(v[k,j])})
            hourly.append(rec)
        del ppart,dpart,acc,vi,ei
    write_json(out/'comparison.json',dict(inputs=checks,max_error=improvement(rows[1]['max_abs_error_m'],rows[0]['max_abs_error_m']),
        rms=improvement(rows[1]['rms_error_m'],rows[0]['rms_error_m']),cases=rows))
    csv_rows(out/'comparison.csv',rows);csv_rows(out/'gap_metrics.csv',gaps);csv_rows(out/'tolerance_exceedance.csv',exceed);csv_rows(out/'hourly_summary.csv',hourly)
    # Figure 2: daily envelopes preserve orbital extremes without plotting aliasing.
    rho=low['density_kg_m3'];mean=rho.mean(axis=1);spread=np.ptp(rho,axis=1)
    days,dm=daily(mean,t);_,lo=daily(mean,t,np.min);_,hi=daily(mean,t,np.max)
    fig,axs=plt.subplots(2,1,figsize=(10,7),sharex=True)
    axs[0].fill_between(days,lo*1e12,hi*1e12,alpha=.2,label='1日内の最小〜最大')
    axs[0].plot(days,dm*1e12,label='全衛星平均密度の日平均');axs[0].set(ylabel='質量密度 [×10^-12 kg/m³]');axs[0].legend()
    _,avg=daily(100*spread/mean,t);_,high=daily(100*spread/mean,t,np.max)
    axs[1].plot(days,avg,label='列全体の密度差比：日平均');axs[1].plot(days,high,label='列全体の密度差比：1日内最大')
    axs[1].set(xlabel='開始からの時間 [日]',ylabel='(最大密度−最小密度) / 平均 [%]');axs[1].legend()
    fig.suptitle('NRLMSIS 2.0：共通密度と差動成分（全LOW・列長2000 km）');fig.tight_layout();fig.savefig(out/'density_common_differential.png',dpi=150);plt.close(fig)
    # Figure 3: include signed gaps and actual q threshold alongside max error.
    fig,axs=plt.subplots(3,1,figsize=(10,9),sharex=True)
    for key,label,arr,mx,rms in series:
        style='--' if key=='baseline' else '-'
        axs[0].plot(t[::30]/86400,mx[::30],label=label,lw=1,ls=style)
        for j in [0,9,19]:axs[1].plot(t[::120]/86400,arr['gap_error_m'][::120,j],label=f'{label}：間隔{j+1}',lw=.8,ls=style)
        # Daily maximum avoids hiding rare q peaks when thinning lines.
        days,peak=daily(arr['q_m_s2'].max(axis=1),t,np.max)
        axs[2].plot(days,peak*1e7,label=f'{label}：日内最大判定量',ls=style)
    axs[0].set(ylabel='最大絶対間隔誤差 [m]');axs[0].legend()
    axs[1].set(ylabel='符号付き間隔誤差 [m]');axs[1].legend(ncol=2,fontsize=8)
    axs[2].axhline(config.controller.threshold_on_m_s2*1e7,color='k',ls='--',label='HIGH起動閾値')
    axs[2].set(ylabel='判定量 [×10^-7 m/s²]',xlabel='開始からの時間 [日]');axs[2].legend()
    fig.suptitle('全LOWと基準制御：60日間の間隔誤差と切替判定');fig.tight_layout();fig.savefig(out/'error_and_q_history.png',dpi=150);plt.close(fig)
    fig,ax=plt.subplots(figsize=(10,4.5))
    for key,label,_,mx,rms in series:ax.plot(t[::30]/86400,rms[::30],label=label,ls='--' if key=='baseline' else '-')
    ax.set(xlabel='開始からの時間 [日]',ylabel='全20間隔の瞬間RMS [m]',title='全LOWと基準制御：各時刻の間隔誤差RMS');ax.legend();fig.tight_layout();fig.savefig(out/'instantaneous_rms.png',dpi=150);plt.close(fig)
    # Figure 5: the same representative gap, raw short response / daily mean long response.
    gap=9;drho=np.diff(low['density_kg_m3'],axis=1)[:,gap]
    drag=-.5*low['density_kg_m3']*config.spacecraft.drag_coefficient*low['area_m2']/config.spacecraft.mass_kg*config.constants.mu_m3_s2/low['a_m']
    dg=np.diff(drag,axis=1)[:,gap]
    fig,axs=plt.subplots(4,2,figsize=(12,9));short=t<=16670
    signals=[drho*1e12,dg*1e8,low['gap_rate_m_s'][:,gap]*1000,low['gap_error_m'][:,gap]]
    labels=['隣接密度差 [×10^-12 kg/m³]','接線抗力差 [×10^-8 m/s²]','間隔速度 [mm/s]','間隔誤差 [m]']
    for axrow,y,label in zip(axs,signals,labels):
        axrow[0].plot(t[short]/3600,y[short],lw=.8);axrow[0].set(ylabel=label)
        days,dy=daily(y,t);axrow[1].plot(days,dy,lw=1);axrow[1].set(ylabel=label)
    axs[0,0].set_title('最初の約3軌道：10秒保存値');axs[0,1].set_title('60日：日平均値')
    axs[-1,0].set_xlabel('開始からの時間 [時間]');axs[-1,1].set_xlabel('開始からの時間 [日]')
    fig.suptitle('密度差 → 抗力差 → 間隔速度 → 間隔誤差：全LOW・衛星10と11\n抗力差と間隔加速度は同じ量ではありません（平均軌道モデル）')
    fig.tight_layout();fig.savefig(out/'density_to_spacing.png',dpi=150);plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,4.5));x=np.arange(3)
    for k,(key,label,*_) in enumerate(series):
        values=[100*r['exceedance_fraction'] for r in comparison[key]['hypothetical_tolerances']]
        ax.bar(x+(k-.5)*.32,values,width=.32,label=label)
        for xx,yy in zip(x+(k-.5)*.32,values):ax.text(xx,yy+.4,f'{yy:.2f}%',ha='center',fontsize=9)
    ax.set_xticks(x,['±100 m','±500 m','±1 km']);ax.set(ylabel='60日中の超過時間割合 [%]',xlabel='仮想的な許容間隔誤差',title='全LOWと基準制御：実ミッション要求ではない仮想許容値',ylim=(0,max(5,ax.get_ylim()[1]+3)))
    ax.legend();fig.tight_layout();fig.savefig(out/'tolerance_exceedance.png',dpi=150);plt.close(fig)
    print(rows,flush=True)


if __name__=='__main__':main()
