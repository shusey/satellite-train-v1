"""Postprocess all full-resolution pairs and publish six Japanese figures."""
import argparse
import csv
import json
import os
from pathlib import Path
import sys
import numpy as np
import yaml

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src')); sys.path.insert(0,str(ROOT/'scripts'))
from run_local_disturbance_n_scaling import OUTPUT, SPEC, old_artifacts
from sat_string.config import load_config
from sat_string.controller import gap_kinematics
from sat_string.n_scaling import gap_graph_coordinates, response_summary
from sat_string.control_cost import event_activity, clean_json
from sat_string.io import write_json


def csv_write(path,rows):
    rows=list(rows)
    if not rows:return
    with path.open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def analyze(n,spec):
    src=OUTPUT/f'N{n:03d}'; out=src/'analysis'
    if (out/'metrics.json').exists():
        m=json.loads((out/'metrics.json').read_text(encoding='utf-8'))
        if 'endpoint_first_1m_s' not in m:
            with (out/'gaps.csv').open(encoding='utf-8-sig') as f:gaps=list(csv.DictReader(f))
            first=[float(row['first_1m_s']) for row in gaps if row['first_1m_s']]
            m['latest_first_1m_s']=max(first) if first else None
            m['endpoint_first_1m_s']=[float(row['first_1m_s']) if row['first_1m_s'] else None
                                     for row in [gaps[0],gaps[-1]]]
            write_json(out/'metrics.json',m)
        return m
    out.mkdir(exist_ok=False)
    cfg=load_config(src/'normalized_config.yaml')
    t=np.load(src/'time_s.npy'); a=np.load(src/'a_m.npy',mmap_mode='r')
    u=np.load(src/'lambda_rad.npy',mmap_mode='r')
    error=np.lib.format.open_memmap(out/'raw_gap_error_m.npy',mode='w+',dtype='float64',shape=(len(t),2,n-1))
    delta=np.lib.format.open_memmap(out/'delta_gap_error_m.npy',mode='w+',dtype='float64',shape=(len(t),n-1))
    maxrate=0.; min_gap=float('inf')
    for start in range(0,len(t),4096):
        sl=slice(start,start+4096)
        e,v,gap=gap_kinematics(a[sl],u[sl],reference_radius_m=cfg.derived.reference_radius_m,
                target_spacing_rad=cfg.derived.target_spacing_rad,mu_m3_s2=cfg.constants.mu_m3_s2)
        error[sl]=e
        delta[sl]=e[:,1]-e[:,0]
        maxrate=max(maxrate,float(np.max(abs(v[:,1]-v[:,0]))))
        min_gap=min(min_gap,float(np.min(gap)))
    error.flush();delta.flush()
    end=spec['input']['start_time_s']+spec['input']['duration_s']
    m=response_summary(t,delta,amplitude=spec['input']['amplitude_m_s2'],
        terminal_window_s=spec['terminal_window_days']*86400.,pulse_end_s=end,
        recovery_band_m=spec['recovery']['error_band_m'],hold_s=spec['recovery']['hold_time_s'])
    r,side=gap_graph_coordinates(n)
    peak=np.array(m.pop('peak_by_gap_m')); tp=np.array(m.pop('peak_time_by_gap_s'))
    gaprows=[]
    for i in range(n-1):
        first=np.flatnonzero(abs(delta[:,i])>=1.)
        gaprows.append(dict(N=n,gap_id=i+1,side=side[i],graph_distance=int(r[i]),
            peak_m=peak[i],peak_time_s=tp[i],peak_at_end=bool(tp[i]==t[-1]),
            first_1m_s=float(t[first[0]]) if len(first) else None,
            nominal_raw_max_m=float(np.max(abs(error[:,0,i]))),
            disturbed_raw_max_m=float(np.max(abs(error[:,1,i])))))
    csv_write(out/'gaps.csv',gaprows)
    spatial=[]
    for distance in range(n//2):
        ix=np.flatnonzero(r==distance); l,ri=ix
        spatial.append(dict(N=n,r=distance,peak_left_m=peak[l],peak_right_m=peak[ri],
            peak_mean_m=float(np.mean(peak[ix])),peak_max_m=float(np.max(peak[ix])),
            peak_left_time_s=tp[l],peak_right_time_s=tp[ri]))
    csv_write(out/'spatial.csv',spatial)
    events=json.loads((src/'checkpoint.json').read_text(encoding='utf-8'))['events']
    nominal=event_activity(events[0],n,float(t[-1])); disturbed=event_activity(events[1],n,float(t[-1]))
    csv_write(out/'nominal_events.csv',events[0]);csv_write(out/'disturbed_events.csv',events[1])
    envelope=np.max(abs(delta),axis=1)
    history=[]
    for k in range(0,len(t),6):
        history.append(dict(time_s=t[k],max_delta_m=envelope[k],instant_rms_delta_m=np.sqrt(np.mean(delta[k]**2)),
            max_nominal_raw_m=np.max(abs(error[k,0])),max_disturbed_raw_m=np.max(abs(error[k,1]))))
    csv_write(out/'history_60s.csv',history)
    p0=max(peak[r==0]); remote=max(peak[r>0])
    prefixes={str(day):float(np.max(envelope[t<=day*86400])) for day in [1,3,7,14]}
    blocks=[]
    for left,right in [(0,1),(1,3),(3,7),(7,11),(11,14)]:
        select=(t>=left*86400)&(t<=right*86400)
        blocks.append(dict(start_day=left,end_day=right,max_response_m=float(np.max(envelope[select])),
            rms_response_m=float(np.sqrt(np.mean(np.asarray(delta[select])**2)))))
    touched=np.flatnonzero(peak>=1.)
    first_times=[row['first_1m_s'] for row in gaprows if row['first_1m_s'] is not None]
    edge_times=[gaprows[0]['first_1m_s'],gaprows[-1]['first_1m_s']]
    m.update(N=n,central_satellite_id=n//2+1,source_gap_ids=[n//2,n//2+1],
        simulation_days=float(t[-1]/86400),terminal_window_days=spec['terminal_window_days'],
        source_peak_m=float(p0),remote_peak_m=float(remote),remote_to_source_ratio=float(remote/p0),
        nominal=nominal,disturbed=disturbed,
        transitions_per_satellite=disturbed['switch_count']/n,
        high_stable_satellite_hours=disturbed['high_stable_total_s']/3600,
        high_stable_hours_per_satellite=disturbed['high_stable_total_s']/3600/n,
        high_equivalent_satellite_hours=disturbed['high_equivalent_total_s']/3600,
        high_equivalent_hours_per_satellite=disturbed['high_equivalent_total_s']/3600/n,
        minimum_signed_gap_m=min_gap,max_delta_gap_rate_m_s=maxrate,
        nominal_raw_max_m=float(np.max(abs(error[:,0]))),
        disturbed_raw_max_m=float(np.max(abs(error[:,1]))),
        maximum_reached_graph_distance_1m=int(max(r[touched])) if len(touched) else None,
        latest_first_1m_s=max(first_times) if first_times else None,
        endpoint_first_1m_s=edge_times,
        end_gap_1m_reached=bool(peak[0]>=1. and peak[-1]>=1.),
        max_mirror_peak_difference_m=float(np.max(abs(peak[:n//2][::-1]-peak[n//2:]))),
        max_response_percent_of_target=m['P_max_m']/cfg.formation.target_spacing_m*100.,
        prefix_peak_m=prefixes,time_blocks=blocks)
    write_json(out/'metrics.json',clean_json(m))
    return clean_json(m)


def publish(metrics,spec,destination):
    destination.mkdir(parents=True,exist_ok=True)
    data=destination/'data';data.mkdir();fig=destination/'figures';fig.mkdir()
    import shutil
    for n in spec['satellite_counts']:
        src=OUTPUT/f'N{n:03d}'
        for source,name in [(src/'normalized_config.yaml',f'N{n:03d}_config.yaml'),
                            (src/'input_checks.json',f'N{n:03d}_input_checks.json'),
                            (src/'completion.json',f'N{n:03d}_completion.json')]:shutil.copyfile(source,data/name)
        for source in (src/'analysis').glob('*.csv'):
            if source.name=='history_60s.csv':
                with source.open(encoding='utf-8-sig') as stream:history=list(csv.DictReader(stream))
                csv_write(data/f'N{n:03d}_history_600s.csv',history[::10])
            else:shutil.copyfile(source,data/f'N{n:03d}_{source.name}')
        shutil.copyfile(src/'analysis/metrics.json',data/f'N{n:03d}_metrics.json')
    rows=[]
    for m in metrics:
        rows.append({k:m[k] for k in ['N','P_max_m','gain_s2','max_response_gap_id','peak_time_s',
            'peak_at_observation_end','E_ult_m','rms_response_m','source_peak_m','remote_to_source_ratio',
            'recovery_candidate_s','transitions_per_satellite','high_stable_satellite_hours',
            'high_stable_hours_per_satellite','high_equivalent_satellite_hours',
            'high_equivalent_hours_per_satellite','maximum_reached_graph_distance_1m','end_gap_1m_reached']}
            | {'transitions':m['disturbed']['switch_count']})
    csv_write(data/'comparison.csv',rows)
    write_json(data/'comparison.json',dict(conditions=spec,results=metrics))
    ns=np.array([m['N'] for m in metrics],float)
    ps=np.array([m['P_max_m'] for m in metrics])
    slope,intercept=np.polyfit(np.log(ns),np.log(ps),1)
    residual=np.log(ps)-(slope*np.log(ns)+intercept)
    total=np.sum((np.log(ps)-np.mean(np.log(ps)))**2)
    r2=float(1-np.sum(residual**2)/total) if total>0 else None
    supported=r2 is not None and r2>=.95 and slope>0 and all(np.diff(ps)>0) and max(ps)/min(ps)>=1.01
    write_json(data/'scaling_assessment.json',dict(
        diagnostic_log_log_R_squared=r2,fit_acceptance='R2>=0.95, positive slope, strictly increasing six peaks, >=1% peak range',
        approximate_power_law_supported=bool(supported),
        exploratory_exponent=float(slope) if supported else None,
        reason='Only six finite-horizon samples; no theoretical scaling law',
        common_horizon_days=spec['duration_days']))
    plot_figures(metrics,spec,fig)


def plot_figures(metrics,spec,fig):
    os.environ['MPLCONFIGDIR']=str(OUTPUT/'matplotlib_cache')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'Yu Gothic','axes.unicode_minus':False,'font.size':10})
    colors=plt.cm.viridis(np.linspace(.05,.9,6))
    vmax=max(m['P_max_m'] for m in metrics)
    f,axes=plt.subplots(3,2,figsize=(13,11),sharex=True,constrained_layout=True)
    for ax,m in zip(axes.ravel(),metrics):
        src=OUTPUT/f"N{m['N']:03d}";t=np.load(src/'time_s.npy')[::6]/86400
        de=np.load(src/'analysis/delta_gap_error_m.npy',mmap_mode='r')[::6]
        im=ax.imshow(de.T,origin='lower',aspect='auto',extent=[0,t[-1],.5,m['N']-.5],
                     cmap='RdBu_r',vmin=-vmax,vmax=vmax,interpolation='nearest')
        ax.set(title=f"衛星数 {m['N']}機",ylabel='間隔番号')
        ax.axvline(4/24,color='gray',ls=':',lw=.7)
    for ax in axes[-1]:ax.set_xlabel('経過時間［日］')
    f.colorbar(im,ax=axes.ravel().tolist(),label='外乱による追加間隔誤差［m］',shrink=.8)
    f.suptitle('中央局所入力の応答：全条件で共通の色範囲')
    f.savefig(fig/'01_heatmaps.png',dpi=170);plt.close(f)
    f,ax=plt.subplots(figsize=(10,6),constrained_layout=True)
    for m,color in zip(metrics,colors):
        with (OUTPUT/f"N{m['N']:03d}/analysis/spatial.csv").open(encoding='utf-8-sig') as stream:
            spatial=list(csv.DictReader(stream))
        rr=[int(x['r']) for x in spatial]
        ax.plot(rr,[float(x['peak_max_m']) for x in spatial],label=f"{m['N']}機：左右最大",color=color)
        ax.plot(rr,[float(x['peak_mean_m']) for x in spatial],color=color,ls=':',lw=.9)
    ax.set(xlabel='外乱源の隣接間隔からの距離［リンク数］',ylabel='最大追加間隔誤差［m］',title='空間的な応答分布（点線は左右平均）')
    ax.set_ylim(0,1.08*vmax)
    ax.legend(ncol=2);ax.grid(alpha=.25);f.savefig(fig/'02_spatial_peaks.png',dpi=170);plt.close(f)
    f,axes=plt.subplots(1,2,figsize=(12,5),constrained_layout=True)
    ns=[m['N'] for m in metrics];ps=[m['P_max_m'] for m in metrics]
    for ax in axes:
        ax.plot(ns,ps,'o-',label='全列最大応答')
        ax.plot(ns,[m['E_ult_m'] for m in metrics],'s--',label='最後の3日間の最大残留影響')
        ax.set(xlabel='衛星数［機］',ylabel='追加間隔誤差［m］');ax.grid(alpha=.25);ax.legend()
    axes[0].set_title('衛星数と最大応答：線形表示')
    axes[0].set_ylim(0,1.12*max(ps))
    axes[1].set(xscale='log',yscale='log',title='衛星数と最大応答：両対数表示')
    axes[1].set_ylim(10,1.3*max(ps))
    f.savefig(fig/'03_n_scaling.png',dpi=170);plt.close(f)
    f,axes=plt.subplots(1,2,figsize=(13,5),constrained_layout=True)
    for m,color in zip(metrics,colors):
        with (OUTPUT/f"N{m['N']:03d}/analysis/gaps.csv").open(encoding='utf-8-sig') as stream:gg=list(csv.DictReader(stream))
        for side,ls in [('left','-'),('right','--')]:
            g=sorted([x for x in gg if x['side']==side],key=lambda x:int(x['graph_distance']))
            axes[0].plot([int(x['graph_distance']) for x in g],[float(x['peak_time_s'])/86400 for x in g],
                ls,color=color,label=f"{m['N']}機" if side=='left' else None,lw=1)
            g=[x for x in g if x['first_1m_s']]
            axes[1].plot([int(x['graph_distance']) for x in g],[float(x['first_1m_s'])/86400 for x in g],ls,color=color,lw=1)
    for ax in axes:ax.set(xlabel='外乱源からの距離［リンク数］',ylabel='経過時間［日］');ax.grid(alpha=.25)
    axes[0].set_title('最大応答の時刻（実線：左／破線：右）');axes[0].legend(ncol=2)
    axes[1].set_title('追加誤差1 mへの初回到達：伝播の補助指標')
    f.savefig(fig/'04_arrival_times.png',dpi=170);plt.close(f)
    f,ax=plt.subplots(figsize=(10,5),constrained_layout=True)
    for m,color in zip(metrics,colors):
        with (OUTPUT/f"N{m['N']:03d}/analysis/history_60s.csv").open(encoding='utf-8-sig') as stream:h=list(csv.DictReader(stream))
        ax.plot([float(x['time_s'])/86400 for x in h],[float(x['max_delta_m']) for x in h],label=f"{m['N']}機",color=color)
    ax.set(xlabel='経過時間［日］',ylabel='全列の最大追加間隔誤差［m］',title='応答包絡線：終盤にも増加しているか')
    ax.legend(ncol=3);ax.grid(alpha=.25);f.savefig(fig/'05_response_envelopes.png',dpi=170);plt.close(f)
    f,axes=plt.subplots(1,2,figsize=(12,5),constrained_layout=True)
    axes[0].plot(ns,[m['disturbed']['switch_count'] for m in metrics],'o-',label='総姿勢遷移回数')
    axes[0].set(xlabel='衛星数［機］',ylabel='総遷移回数',title='制御負荷：総数')
    a1=axes[1];a2=a1.twinx()
    a1.plot(ns,[m['transitions_per_satellite'] for m in metrics],'o-',label='1機あたり遷移回数')
    a2.plot(ns,[m['high_equivalent_hours_per_satellite'] for m in metrics],'s--',color='tab:orange',label='1機あたりHIGH換算時間')
    a1.set(xlabel='衛星数［機］',ylabel='遷移回数［回／機］',title='制御負荷：衛星数で規格化')
    a2.set_ylabel('HIGH換算時間［時間／機］')
    lines=a1.lines+a2.lines;a1.legend(lines,[x.get_label() for x in lines])
    axes[0].legend()
    for ax in axes:ax.grid(alpha=.25)
    f.savefig(fig/'06_control_burden.png',dpi=170);plt.close(f)


def main():
    p=argparse.ArgumentParser();p.add_argument('--publish',type=Path)
    args=p.parse_args();spec=yaml.safe_load(SPEC.read_text(encoding='utf-8'))
    metrics=[analyze(n,spec) for n in spec['satellite_counts']]
    before=json.loads((OUTPUT/'old_artifacts_before.json').read_text(encoding='utf-8'))
    now=old_artifacts()
    changed=[k for k,v in before['files'].items() if now['files'].get(k)!=v]
    assert not changed,changed
    assert before['fixture_sha256']==now['fixture_sha256']
    write_json(OUTPUT/'non_destruction.json',dict(passed=True,old_files_checked=len(before['files']),
        fixture_hashes_unchanged=True))
    write_json(OUTPUT/'comparison.json',dict(conditions=spec,results=metrics))
    if args.publish:publish(metrics,spec,args.publish)
    print(json.dumps([{k:m[k] for k in ['N','P_max_m','E_ult_m','remote_to_source_ratio','peak_time_s','end_gap_1m_reached']} for m in metrics],indent=2))


if __name__=='__main__':main()
