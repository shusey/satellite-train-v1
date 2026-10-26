"""Join unchanged 0-14d prefixes with new 14-30d suffixes, then compare."""
import argparse
import csv
import json
import os
from pathlib import Path
import shutil
import sys
import numpy as np
import yaml

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'scripts'))
from extend_local_disturbance_n_scaling import OUTPUT,SOURCE,old_artifacts
from run_local_disturbance_n_scaling import SPEC
from analyze_local_disturbance_n_scaling import analyze,csv_write
from sat_string.control_cost import event_activity,clean_json
from sat_string.io import write_json
from sat_string.regression import compare_baseline

DEST=ROOT/'docs/local_disturbance_n_scaling_30d'


def symmetry_diagnostics(n,m):
    """Record tiny mirror differences and threshold-adjacent command branching."""
    from sat_string.config import load_config
    from sat_string.controller import gap_kinematics,switching_function
    folder=OUTPUT/f'N{n:03d}'
    cp=json.loads((folder/'checkpoint.json').read_text(encoding='utf-8'))
    events=cp['events'][1];mismatches=[];unequal=[]
    for j in range(1,n//2+1):
        mirror=n+1-j
        left=[e for e in events if e['satellite_id']==j]
        right=[e for e in events if e['satellite_id']==mirror]
        if len(left)!=len(right):unequal.append([j,mirror,len(left),len(right)])
        for k,(a,b) in enumerate(zip(left,right)):
            if a['start_time']!=b['start_time']:
                mismatches.append(dict(satellite_ids=[j,mirror],command_index=k,
                    first_start_s=min(a['start_time'],b['start_time']),
                    start_difference_s=a['start_time']-b['start_time']))
    mismatches.sort(key=lambda row:row['first_start_s'])
    first=None
    if mismatches:
        first=dict(mismatches[0]);cfg=load_config(folder/'normalized_config.yaml')
        part=folder if first['first_start_s']>=14*86400 else SOURCE/f'N{n:03d}'
        t=np.load(part/'time_s.npy');k=np.searchsorted(t,first['first_start_s'])
        a=np.load(part/'a_m.npy',mmap_mode='r')[k,1];u=np.load(part/'lambda_rad.npy',mmap_mode='r')[k,1]
        e,v,_=gap_kinematics(a,u,reference_radius_m=cfg.derived.reference_radius_m,
            target_spacing_rad=cfg.derived.target_spacing_rad,mu_m3_s2=cfg.constants.mu_m3_s2)
        q=switching_function(e,v,cfg.derived.controller_kp_s2_inv,cfg.derived.controller_kd_s_inv)
        first['q_minus_h_on_m_s2']=[float(q[j-1]-cfg.controller.threshold_on_m_s2) for j in first['satellite_ids']]
    return dict(N=n,max_mirror_peak_difference_m=m['max_mirror_peak_difference_m'],
        unequal_mirror_command_counts=unequal,mirror_command_start_mismatches=len(mismatches),
        maximum_start_difference_s=max(abs(row['start_difference_s']) for row in mismatches) if mismatches else 0,
        first_mismatch=first)


def analyze_pair(n,spec):
    folder=OUTPUT/f'N{n:03d}';out=folder/'analysis'
    m=analyze(n,spec,output=OUTPUT,parts=[SOURCE/f'N{n:03d}',folder])
    prior=json.loads((SOURCE/f'N{n:03d}/analysis/metrics.json').read_text(encoding='utf-8'))
    assert m['prefix_peak_m']['14']==prior['P_max_m']
    cp=json.loads((folder/'checkpoint.json').read_text(encoding='utf-8'))
    original=json.loads((SOURCE/f'N{n:03d}/checkpoint.json').read_text(encoding='utf-8'))
    for before,after in zip(original['events'],cp['events']):assert after[:len(before)]==before
    t=np.load(out/'time_s.npy');delta=np.load(out/'delta_gap_error_m.npy',mmap_mode='r')
    np.testing.assert_array_equal(delta[:len(np.load(SOURCE/f'N{n:03d}/time_s.npy'))],
        np.load(SOURCE/f'N{n:03d}/analysis/delta_gap_error_m.npy',mmap_mode='r'))
    old_events=cp['events'][1]
    windows=[]
    for left,right in [(0,14),(14,21),(21,27),(27,30),(14,30)]:
        sl=(t>=left*86400)&(t<=right*86400)
        d=np.asarray(delta[sl]);a=event_activity(old_events,n,right*86400)
        b=event_activity(old_events,n,left*86400)
        # Starts at the left boundary belong to this interval; final end is included.
        count=sum(left*86400<=e['start_time']<right*86400 for e in old_events)
        cycles=sum(left*86400<=e['complete_time']<right*86400 and e['to_mode']=='L' for e in old_events)
        if right==30:
            count+=sum(e['start_time']==right*86400 for e in old_events)
            cycles+=sum(e['complete_time']==right*86400 and e['to_mode']=='L' for e in old_events)
        windows.append(dict(N=n,start_day=left,end_day=right,max_delta_m=float(np.max(abs(d))),
            rms_delta_m=float(np.sqrt(np.mean(d*d))),transitions=count,cycles=cycles,
            transitions_per_day=count/(right-left),cycles_per_day=cycles/(right-left),
            high_equivalent_satellite_hours=(a['high_equivalent_total_s']-b['high_equivalent_total_s'])/3600))
    csv_write(out/'windows.csv',windows)
    series=[]
    for end in [14,21,30]:
        sl=t<=end*86400
        act=event_activity(old_events,n,end*86400)
        series.append(dict(N=n,observed_days=end,P_max_m=float(np.max(abs(delta[sl]))),
            transitions=act['switch_count'],high_equivalent_satellite_hours=act['high_equivalent_total_s']/3600))
    csv_write(out/'horizon_comparison.csv',series)
    m['extension_windows']=windows
    m['horizon_comparison']=series
    m['new_peak_gain_vs_14d_percent']=100*(m['P_max_m']/prior['P_max_m']-1)
    m['events_after_14d']=len(cp['events'][1])-len(original['events'][1])
    m['completed_cycles']=sum(e['to_mode']=='L' and e['complete_time']<=30*86400 for e in old_events)
    m['LOW_to_HIGH_starts']=sum(e['to_mode']=='H' and e['start_time']<=30*86400 for e in old_events)
    m['HIGH_to_LOW_starts']=sum(e['to_mode']=='L' and e['start_time']<=30*86400 for e in old_events)
    m['last_transition_complete_s']=max(e['complete_time'] for e in old_events if e['complete_time']<=30*86400)
    m['last_transition_start_s']=max(e['start_time'] for e in old_events if e['start_time']<=30*86400)
    m['final_modes']=cp['modes'][1]
    write_json(out/'extended_metrics.json',clean_json(m))
    return clean_json(m)


def figures(metrics):
    os.environ['MPLCONFIGDIR']=str(OUTPUT/'matplotlib_cache')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'Yu Gothic','axes.unicode_minus':False,'font.size':10})
    fig=DEST/'figures';fig.mkdir(exist_ok=True)
    colors=['tab:blue','tab:orange'];vmax=max(m['P_max_m'] for m in metrics)
    f,axes=plt.subplots(2,1,figsize=(12,8),sharex=True,constrained_layout=True)
    for m,col in zip(metrics,colors):
        folder=OUTPUT/f"N{m['N']:03d}/analysis"
        with (folder/'history_60s.csv').open(encoding='utf-8-sig') as stream:h=list(csv.DictReader(stream))
        t=[float(r['time_s'])/86400 for r in h]
        style='--' if m['N']==301 else '-'
        axes[0].plot(t,[float(r['max_delta_m']) for r in h],style,color=col,label=f"{m['N']}機")
        axes[1].plot(t,[float(r['instant_rms_delta_m']) for r in h],style,color=col,label=f"{m['N']}機")
    for ax in axes:
        ax.axvline(14,color='gray',ls='--',label='既存計算との接続：14日')
        ax.grid(alpha=.25);ax.legend()
    axes[0].set(ylabel='最大追加間隔誤差［m］',title='中央局所入力への応答：14日から30日へ延長')
    axes[1].set(xlabel='経過時間［日］',ylabel='瞬間RMS［m］')
    f.savefig(fig/'01_error_history.png',dpi=170);plt.close(f)
    f,axes=plt.subplots(2,1,figsize=(12,8),sharex=True,constrained_layout=True)
    for ax,m in zip(axes,metrics):
        folder=OUTPUT/f"N{m['N']:03d}/analysis";t=np.load(folder/'time_s.npy')[::6]/86400
        d=np.load(folder/'delta_gap_error_m.npy',mmap_mode='r')[::6]
        im=ax.imshow(d.T,origin='lower',aspect='auto',extent=[0,30,.5,m['N']-.5],cmap='RdBu_r',vmin=-vmax,vmax=vmax,interpolation='nearest')
        ax.axvline(14,color='gray',ls='--',lw=.8);ax.set(title=f"衛星数 {m['N']}機",ylabel='間隔番号')
    axes[-1].set_xlabel('経過時間［日］');f.colorbar(im,ax=list(axes),label='追加間隔誤差［m］',shrink=.8)
    f.savefig(fig/'02_heatmaps.png',dpi=170);plt.close(f)
    f,axes=plt.subplots(1,2,figsize=(13,5),constrained_layout=True)
    for m,col in zip(metrics,colors):
        n=m['N']
        for day,root,style in [(14,SOURCE,'--'),(30,OUTPUT,'-')]:
            with (root/f'N{n:03d}/analysis/spatial.csv').open(encoding='utf-8-sig') as stream:rows=list(csv.DictReader(stream))
            axes[0].plot([int(r['r']) for r in rows],[float(r['peak_max_m']) for r in rows],style,color=col,label=f'{n}機・{day}日')
        with (OUTPUT/f'N{n:03d}/analysis/gaps.csv').open(encoding='utf-8-sig') as stream:rows=list(csv.DictReader(stream))
        for side,style in [('left','-'),('right','--')]:
            rr=sorted([r for r in rows if r['side']==side],key=lambda r:int(r['graph_distance']))
            axes[1].plot([int(r['graph_distance']) for r in rr],[float(r['peak_time_s'])/86400 for r in rr],style,color=col,label=f'{n}機・'+('左' if side=='left' else '右'))
    axes[0].set(xlabel='外乱源からの距離［リンク数］',ylabel='最大追加間隔誤差［m］',title='空間最大応答：14日と30日の比較')
    axes[0].set_ylim(0,1.08*vmax)
    axes[1].set(xlabel='外乱源からの距離［リンク数］',ylabel='最大応答時刻［日］',title='30日間の最大応答時刻：初回到達とは区別')
    for ax in axes:ax.legend();ax.grid(alpha=.25)
    f.savefig(fig/'03_spatial_response.png',dpi=170);plt.close(f)
    f,axes=plt.subplots(1,2,figsize=(13,5),constrained_layout=True)
    for m,col in zip(metrics,colors):
        cp=json.loads((OUTPUT/f"N{m['N']:03d}/checkpoint.json").read_text(encoding='utf-8'))
        ev=cp['events'][1];times=np.sort([e['start_time']/86400 for e in ev])
        axes[0].step(np.r_[0,times,30],np.r_[0,np.arange(1,len(times)+1),len(times)],where='post',color=col,label=f"{m['N']}機・総数")
        axes[1].step(np.r_[0,times,30],np.r_[0,np.arange(1,len(times)+1),len(times)]/m['N'],where='post',color=col,label=f"{m['N']}機・1機あたり")
    for ax in axes:ax.axvline(14,color='gray',ls='--');ax.set_xlabel('経過時間［日］');ax.legend();ax.grid(alpha=.25)
    axes[0].set(ylabel='累積遷移回数［回］',title='制御活動は停止したか：総数')
    axes[1].set(ylabel='累積遷移回数［回／機］',title='衛星数で規格化した負荷')
    f.savefig(fig/'04_switching_activity.png',dpi=170);plt.close(f)


def main():
    p=argparse.ArgumentParser();p.add_argument('--analyze-only',action='store_true');args=p.parse_args()
    spec=yaml.safe_load(SPEC.read_text(encoding='utf-8'));spec['satellite_counts']=[21,301];spec['duration_days']=30
    metrics=[analyze_pair(n,spec) for n in [21,301]]
    write_json(OUTPUT/'comparison.json',dict(conditions=spec,results=metrics))
    if args.analyze_only:
        print(json.dumps(metrics,indent=2));return
    DEST.mkdir(parents=True,exist_ok=True);data=DEST/'data';data.mkdir(exist_ok=False)
    for m in metrics:
        n=m['N'];folder=OUTPUT/f'N{n:03d}'
        for name in ['normalized_config.yaml','completion.json','input_checks.json']:
            shutil.copyfile(folder/name,data/f'N{n:03d}_{name}')
        for file in (folder/'analysis').glob('*.csv'):
            if file.name=='history_60s.csv':
                with file.open(encoding='utf-8-sig') as stream:rows=list(csv.DictReader(stream))
                csv_write(data/f'N{n:03d}_history_600s.csv',rows[::10])
            else:shutil.copyfile(file,data/f'N{n:03d}_{file.name}')
        for name in ['extended_metrics.json','join_checks.json']:shutil.copyfile(folder/'analysis'/name,data/f'N{n:03d}_{name}')
    write_json(data/'comparison.json',dict(conditions=spec,results=metrics))
    write_json(data/'symmetry_diagnostics.json',dict(results=[symmetry_diagnostics(m['N'],m) for m in metrics]))
    csv_write(data/'comparison.csv',[dict(N=m['N'],P_max_m=m['P_max_m'],gain_s2=m['gain_s2'],
        E_ult_m=m['E_ult_m'],final_max_m=m['final_max_response_m'],rms_m=m['rms_response_m'],
        peak_day=m['peak_time_s']/86400,peak_gap=m['max_response_gap_id'],
        remote_to_source_ratio=m['remote_to_source_ratio'],transitions=m['disturbed']['switch_count'],
        transitions_per_satellite=m['transitions_per_satellite'],high_eq_sat_h=m['high_equivalent_satellite_hours'],
        high_stable_sat_h=m['high_stable_satellite_hours'],recovery_candidate_s=m['recovery_candidate_s']) for m in metrics])
    figures(metrics)
    before=json.loads((OUTPUT/'old_artifacts_before.json').read_text(encoding='utf-8'));now=old_artifacts(OUTPUT)
    assert all(now['files'].get(k)==v for k,v in before['files'].items())
    assert now['fixture_sha256']==before['fixture_sha256']
    reg=compare_baseline(SOURCE/'regression_24h/suite',ROOT/'tests/fixtures/regression_24h');assert reg['passed']
    write_json(data/'regression_24h.json',dict(**reg,execution='Rechecked existing freshly-generated 14-day-study baseline suite; physical core unchanged'))
    log=(OUTPUT/'tests.txt').read_text(encoding='utf-8');assert '345 passed' in log and 'failed' not in log
    shutil.copyfile(OUTPUT/'tests.txt',data/'tests.txt')
    write_json(data/'verification.json',dict(passed=True,tests_passed=345,new_tests=10,
        unchanged_old_files=len(before['files']),frozen_fixture_hashes_unchanged=True,
        exact_checkpoint_boundaries=True,old_14d_response_exact=True,old_event_prefix_exact=True,
        pulse_replayed=False,fresh_physics_days=[14,30],nominal_transitions=[m['nominal']['switch_count'] for m in metrics],
        regression_24h_passed=True,physical_core_changed=False))
    print(json.dumps([{k:m[k] for k in ['N','P_max_m','E_ult_m','rms_response_m','remote_to_source_ratio','peak_time_s']} for m in metrics],indent=2))


if __name__=='__main__':main()
