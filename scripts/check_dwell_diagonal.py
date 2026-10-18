"""Six-hour explicit split-dwell trajectories versus saved scalar-dwell runs."""
import json
from run_reference_dwell_factorial import ROOT, OUTPUT, SOURCE, PROFILE, SPEC
import numpy as np
import yaml
from run_nrlmsis_case import make_environment
from sat_string.config import load_config
from sat_string.reference_condition import load_frozen
from sat_string.dwell_factorial import factorial_config
from sat_string.dispersion_control import integrate_control_profiles
from sat_string.io import write_json


def main():
    report={}; spec=yaml.safe_load(SPEC.read_text()); profile=load_frozen(PROFILE)['profiles'][0]
    base=load_config(SOURCE/'normalized_config.yaml')
    for dwell,source,member in [(300,SOURCE,0),(150,ROOT/'outputs/reference_a_dwell_20261017/D150_7d',0)]:
        cfg=factorial_config(base,dwell,dwell,days=.25)
        _,_,provider=make_environment(ROOT/spec['environment'],ROOT/spec['weather'],spec['cache_seconds'])
        out=OUTPUT/f'diagonal_{dwell}_6h'
        _,events=integrate_control_profiles(cfg,[profile],provider,out)
        time=np.load(out/'time_s.npy'); checks={}
        for key in ['a_m','lambda_rad','area_m2','mode','density_kg_m3']:
            actual=np.load(out/f'{key}.npy')[:,0]
            original=np.load(source/f'{key}.npy',mmap_mode='r')[:len(time),member]
            checks[key]=dict(exact=bool(np.array_equal(actual,original)),max_abs_difference=float(np.max(abs(actual.astype(float)-original.astype(float)))))
        expected=[e for e in json.loads((source/'checkpoint.json').read_text())['events'][member] if e['command_time']<=time[-1]]
        checks['events_exact']=events[0]==expected
        checks['passed']=all(x['exact'] for x in checks.values() if isinstance(x,dict)) and checks['events_exact']
        report[str(dwell)]=checks
        print(dwell,checks,flush=True)
        assert checks['passed']
    write_json(OUTPUT/'diagonal_regression.json',report)


if __name__=='__main__':main()
