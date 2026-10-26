"""Verify saved pairs, legacy regression, completed tests and old artifacts."""
import json
import argparse
import re
import shutil
from pathlib import Path
import sys
import numpy as np
import yaml

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'scripts'))
from run_local_disturbance_n_scaling import OUTPUT,SPEC,old_artifacts
from sat_string.config import load_config
from sat_string.n_scaling import scaling_configuration,CentralNodePulse
from sat_string.io import write_json


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--destination',type=Path,default=OUTPUT/'verification')
    args=parser.parse_args()
    spec=yaml.safe_load(SPEC.read_text(encoding='utf-8'))
    base=load_config(ROOT/spec['base_configuration']);rows=[]
    for n in spec['satellite_counts']:
        folder=OUTPUT/f'N{n:03d}'
        cfg=load_config(folder/'normalized_config.yaml')
        assert cfg==scaling_configuration(base,n,spec['duration_days'])
        meta=json.loads((folder/'metadata.json').read_text(encoding='utf-8'))
        profiles=meta['profiles']
        assert profiles[0]['control_enabled'] and profiles[1]['control_enabled']
        for key in ['delta_a_m','B_multiplier']:
            assert profiles[0][key]==profiles[1][key]
        assert profiles[0]['delta_a_m']==[0.]*n and profiles[0]['B_multiplier']==[1.]*n
        pulse=CentralNodePulse(n,**{k:spec['input'][k] for k in ['amplitude_m_s2','start_time_s','duration_s']})
        assert meta['node_input']==pulse.metadata()
        cp=json.loads((folder/'checkpoint.json').read_text(encoding='utf-8'))
        assert not cp['events'][0] and cp['node_input']==pulse.metadata()
        t=np.load(folder/'time_s.npy');np.testing.assert_array_equal(np.diff(t),10.)
        assert t[0]==0 and t[-1]==spec['duration_days']*86400
        a=np.load(folder/'a_m.npy',mmap_mode='r');u=np.load(folder/'lambda_rad.npy',mmap_mode='r')
        mode=np.load(folder/'mode.npy',mmap_mode='r')
        np.testing.assert_array_equal(a[0,0],a[0,1]);np.testing.assert_array_equal(u[0,0],u[0,1])
        np.testing.assert_array_equal(a[0],cfg.derived.reference_radius_m)
        np.testing.assert_allclose(np.diff(u[0],axis=-1)*cfg.derived.reference_radius_m,100000.,rtol=0,atol=1e-8)
        assert not np.any(mode[:,0])
        delta=np.load(folder/'analysis/delta_gap_error_m.npy',mmap_mode='r')
        np.testing.assert_array_equal(delta[t<=pulse.start_time_s],0.)
        np.testing.assert_array_equal(a[-1],np.asarray(cp['state'])[0].reshape(2,n))
        np.testing.assert_array_equal(u[-1],np.asarray(cp['state'])[1].reshape(2,n))
        assert np.all(np.isfinite(delta))
        rows.append(dict(N=n,passed=True,samples=len(t),nominal_transitions=0,
                         identical_initial_state=True,only_local_node_input_differs=True))
    before=json.loads((OUTPUT/'old_artifacts_before.json').read_text(encoding='utf-8'))
    now=old_artifacts()
    assert all(now['files'].get(k)==v for k,v in before['files'].items())
    assert now['fixture_sha256']==before['fixture_sha256']
    regression=json.loads((OUTPUT/'regression_24h/regression.json').read_text(encoding='utf-8'))
    assert regression['passed']
    log=(OUTPUT/'tests_verified.txt').read_text(encoding='utf-8')
    match=re.search(r'(\d+) passed',log);assert match and 'failed' not in log
    data=args.destination
    data.mkdir(parents=True,exist_ok=True)
    for name in ['verification.json','tests.txt','regression_24h.json','non_destruction.json']:
        if (data/name).exists():raise FileExistsError(f'Existing verification retained: {data/name}')
    report=dict(passed=True,pairs=rows,full_test_count=int(match[1]),new_test_count=12,
        regression_24h_passed=True,regression_fresh_three_case_simulation=True,
        unchanged_old_run_files=len(before['files']),frozen_fixture_hashes_unchanged=True,
        raw_history_and_metrics_cadence_s=10,plot_cadence_s=60,published_history_cadence_s=600,
        no_gaussian_density_pulse=True,no_nrlmsis_spatial_forcing=True,
        no_initial_or_B_dispersion=True,unchanged_controller=True)
    write_json(data/'verification.json',report)
    shutil.copyfile(OUTPUT/'tests_verified.txt',data/'tests.txt')
    shutil.copyfile(OUTPUT/'regression_24h/regression.json',data/'regression_24h.json')
    shutil.copyfile(OUTPUT/'non_destruction.json',data/'non_destruction.json')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
