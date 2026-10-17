"""Validate offline attribution, compare saved 24h reference, and publish small tables."""
from pathlib import Path
import argparse,json,sys,shutil,hashlib,subprocess
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from sat_string.regression import compare_baseline
from sat_string.pulse_diagnostics import acceleration_terms, cumulative_trapezoid
from analyze_reference_pulses import read_csv,write_json,write_csv,inputs_manifest,SOURCE,OLD,FROZEN


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',default='outputs/reference_pulse_mechanism_20261017_v3')
    ap.add_argument('--output',default='docs/reference_pulse_mechanism/data_final');args=ap.parse_args()
    src=ROOT/args.input;out=ROOT/args.output;out.mkdir(parents=True,exist_ok=False)
    before=json.loads((src/'input_manifest.json').read_text(encoding='utf-8'))
    assert inputs_manifest()==before
    report=compare_baseline(SOURCE/'zero_regression_24h/suite',ROOT/'tests/fixtures/regression_24h')
    assert report['passed'];report['execution']='saved previously regenerated 24h suite compared again; no new simulation'
    write_json(out/'regression_24h.json',report)
    summary=json.loads((src/'summary.json').read_text(encoding='utf-8'));detail={};sampling=[]
    profiles=json.loads(FROZEN.read_text(encoding='utf-8'))
    for ic,case in enumerate(['A','B','AB']):
        rows=read_csv(src/f'{case}_pulse_cycles.csv');q=read_csv(OLD/f'{case}_quiet_cycles.csv');old=read_csv(OLD/f'{case}_local_cycles.csv')
        matched=[abs(float(x['validated_prediction_error_s'])) for x in q if x['prediction_valid_on_observed_window']=='True' and x['prediction_target_matches_actual']=='True' and x['validated_prediction_error_s']]
        detail[case]=dict(HIGH_completion_already_below_off_fraction=float(np.mean([float(x['q_at_HIGH_complete_m_s2'])<=1.4e-7 for x in old])),
            first_off_during_up_slew_fraction=float(np.mean([float(x['first_q_below_off_s'])<float(x['high_complete_s']) for x in old if x['first_q_below_off_s']])),
            validated_target_matched_forecasts=len(matched),matched_prediction_median_abs_error_s=float(np.median(matched)))
        assert len(rows)==summary['cases'][case]['cycles']
        for part,limit in [('batch_30d',30*86400),('batch_60d',60*86400)]:
            times=np.load(SOURCE/part/'time_s.npy')
            arrays={k:np.load(SOURCE/part/f'{k}.npy',mmap_mode='r')[:,ic] for k in ['a_m','density_kg_m3','area_m2']}
            eligible=[x for x in rows if float(x['start_s'])>=times[0] and float(x['end_s'])<=times[-1]]
            # Deterministic stratified sample, not selected by result.
            ix=np.unique(np.linspace(0,len(eligible)-1,50,dtype=int))
            for z in ix:
                row=eligible[z];t0=float(row['start_s']);t1=float(row['end_s']);k0=int(np.searchsorted(times,t0));k1=int(np.searchsorted(times,t1))
                a=arrays['a_m'][k0:k1+1];rho=arrays['density_kg_m3'][k0:k1+1];area=arrays['area_m2'][k0:k1+1];t=times[k0:k1+1]
                terms=acceleration_terms(a,rho,area,np.array(profiles['profiles'][ic]['epsilon']),r_ref=6778137.,mu=398600441800000.,cd=2.2,mass=4.,area_low=.01)
                y=sum(terms.values());fine=np.diff(cumulative_trapezoid(y,t)[-1]);sel=np.unique(np.r_[np.arange(0,len(t),2),len(t)-1])
                coarse=np.diff(cumulative_trapezoid(y[sel],t[sel])[-1]);v=6778137.*np.diff(np.sqrt(398600441800000./a**3),axis=1);truth=v[-1]-v[0]
                sampling.append(dict(case=case,start_s=t0,end_s=t1,
                    fine_10s_max_closure_m_s=float(np.max(abs(fine-truth))),
                    coarse_20s_max_closure_m_s=float(np.max(abs(coarse-truth))),
                    fine_coarse_max_difference_m_s=float(np.max(abs(fine-coarse)))))
        # Independent closure acceptance is deliberately much smaller than ~1 mm/s.
        assert summary['cases'][case]['pulse_closure_max_m_s']<1e-5
        # Compact representative navigation table, totals remain in summary.
        chosen=sorted(set(range(min(3,len(rows))))|set(np.linspace(0,len(rows)-1,30,dtype=int)))
        write_csv(out/f'{case}_pulse_samples.csv',[rows[k] for k in chosen])
    assert inputs_manifest()==before
    write_json(out/'additional_diagnostics.json',detail);write_csv(out/'quadrature_check.csv',sampling)
    for name in ['summary.json','interval_summary.csv','examples.json']:
        shutil.copyfile(src/name,out/name)
    core_files=subprocess.check_output(['git','ls-files','src/sat_string','configs','tests/fixtures/regression_24h'],cwd=ROOT,text=True).splitlines()
    unchanged=subprocess.check_output(['git','diff','--name-only','HEAD','--','src/sat_string','configs','tests/fixtures/regression_24h'],cwd=ROOT,text=True).strip()==''
    assert unchanged
    verify=dict(passed=True,input_files_unchanged=len(before),existing_tracked_core_config_fixture_unchanged=True,
        new_physics_runs=0,regression_comparison_passed=True,new_unit_tests=5,
        pulse_closure_tolerance_m_s=1e-5,pulse_max_observed_closure_m_s=max(x['pulse_closure_max_m_s'] for x in summary['cases'].values()),
        quadrature_cases=len(sampling),quadrature_fine_median_m_s=float(np.median([x['fine_10s_max_closure_m_s'] for x in sampling])),
        quadrature_coarse_median_m_s=float(np.median([x['coarse_20s_max_closure_m_s'] for x in sampling])),
        preexisting_core_config_fixture_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in core_files})
    write_json(out/'verification.json',verify)
    print(json.dumps(detail,ensure_ascii=False));print('Verification passed',flush=True)


if __name__=='__main__':main()
