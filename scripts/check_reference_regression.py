"""Fresh zero-dispersion 24h adapter suite versus frozen software reference."""
import csv,json
import numpy as np
from run_reference_dispersion import ROOT,DATA
from sat_string.config import load_config,dump_normalized_config
from sat_string.dispersion_control import integrate_control_profiles,derived_histories
from sat_string.io import write_json
from sat_string.regression import compare_baseline


def main():
    cfg=load_config(ROOT/'tests/fixtures/regression_24h/config.yaml')
    out=ROOT/'outputs/reference_dispersion_20261012/zero_regression_24h'
    out.mkdir(parents=True,exist_ok=False)
    for disturbed,names in [(True,['controlled_disturbed','all_low_drag_baseline']),
                             (False,['controlled_undisturbed'])]:
        profiles=[dict(delta_a_m=[0.]*21,B_multiplier=[1.]*21,
            control_enabled=name.startswith('controlled')) for name in names]
        path=out/('disturbed' if disturbed else 'undisturbed')
        _,events=integrate_control_profiles(cfg,profiles,None,path,disturbance_enabled=disturbed)
        t=np.load(path/'time_s.npy')
        for j,name in enumerate(names):
            folder=out/'suite'/name;folder.mkdir(parents=True)
            a,u,area,mode,rho=[np.load(path/f'{key}.npy',mmap_mode='r')[:,j] for key in
                ('a_m','lambda_rad','area_m2','mode','density_kg_m3')]
            arrays=derived_histories(t,a,u,area,mode,rho,cfg)
            np.savez_compressed(folder/'results.npz',**arrays)
            dump_normalized_config(cfg,folder/'normalized_config.yaml')
            write_json(folder/'metadata.json',dict(case=dict(name=name,
                control_enabled=name.startswith('controlled'),disturbance_enabled=disturbed)))
            with (folder/'switch_events.csv').open('w',newline='',encoding='utf-8') as f:
                writer=csv.DictWriter(f,fieldnames=['satellite_id','command_time','start_time','complete_time','from_mode','to_mode'])
                writer.writeheader();writer.writerows(events[j])
    report=compare_baseline(out/'suite',ROOT/'tests/fixtures/regression_24h')
    report['zero_dispersion_adapter']=True
    write_json(DATA/'regression_24h.json',report)
    print('24h zero-dispersion reference:',report['passed'],flush=True)
    if not report['passed']:raise SystemExit(1)


if __name__=='__main__':main()
