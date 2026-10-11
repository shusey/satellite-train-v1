"""Publish compact input, test and preservation evidence; no propagation."""
import hashlib,json,re,shutil
from pathlib import Path
import numpy as np
import pymsis
from run_dispersion_screening import ROOT,OUT,DATA,configuration


def main():
    spec,cfg=configuration(7)
    reference=json.loads((OUT/'initial_7d'/'metadata.json').read_text())['configuration']
    old=__import__('sat_string.config',fromlist=['load_config']).load_config(ROOT/spec['baseline_run']/'normalized_config.yaml').to_dict()
    stages=[]
    for name in ('initial_7d','extended_30d','combination_7d'):
        path=OUT/name
        if not (path/'completion.json').exists():continue
        m=json.loads((path/'metadata.json').read_text())
        checked=m['configuration'].copy();checked['simulation']=dict(checked['simulation'])
        duration=checked['simulation']['duration_s'];checked['simulation']['duration_s']=old['simulation']['duration_s']
        if checked!=old or not m['all_low'] or m['control_enabled']:raise ValueError('Unexpected input changes')
        env=json.loads((path/'environment.json').read_text())
        previous=json.loads((ROOT/spec['baseline_run']/'density_environment.json').read_text())
        if env['weather_sha256']!=previous['weather_sha256'] or env['spec']!=previous['environment']:raise ValueError('Weather/geometry mismatch')
        profiles=m['profiles']
        for p in profiles:
            if p['factor']=='a' and np.any(np.asarray(p['B_multiplier'])!=1):raise ValueError('a-only profile changes B')
            if p['factor']=='B' and np.any(np.asarray(p['delta_a_m'])!=0):raise ValueError('B-only profile changes initial a')
        stages.append(dict(name=name,profiles=len(profiles),duration_s=duration,start_s=m['start_s'],
            configuration_changes=['simulation.duration_s'],independent_profile_changes=['initial_a_m','B_multiplier'],
            maximum_mean_delta_a_m=max(abs(np.mean(p['delta_a_m'])) for p in profiles),
            maximum_mean_B_fraction_error=max(abs(np.mean(p['B_multiplier'])-1) for p in profiles),
            minimum_B_multiplier=min(min(p['B_multiplier']) for p in profiles),all_low=True,
            weather_sha256=env['weather_sha256']))
    tests=(OUT/'tests_final.log').read_text()
    match=re.search(r'(\d+) passed in ([\d.]+)s',tests)
    if match is None or 'failed' in tests:raise ValueError('Full test suite not passed')
    final=OUT/'preservation_final'/'preservation.json'
    preservation=final if final.exists() else OUT/'preservation_check'/'preservation.json'
    for original,dest in [(OUT/'regression_24h'/'regression.json','regression_24h.json'),
                          (preservation,'preservation.json')]:
        data=json.loads(original.read_text())
        if not data['passed']:raise ValueError('Regression/preservation failed')
        shutil.copyfile(original,DATA/dest)
    parameter=Path(pymsis.__file__).parent/'msis2.0.parm'
    proof=dict(stages=stages,reference_config=old,packages=dict(pymsis=pymsis.__version__,numpy=np.__version__),
        msis_version='2.0',msis20_parameter_sha256=hashlib.sha256(parameter.read_bytes()).hexdigest(),
        tests=dict(passed=int(match[1]),new_tests=11,elapsed_s=float(match[2]),log=tests),
        baseline_reused=True,core_modifications=False,
        preservation_report='preservation.json',regression_report='regression_24h.json')
    (DATA/'verification.json').write_text(json.dumps(proof,indent=2,ensure_ascii=False),encoding='utf-8')
    print('Published input checks and',match[1],'passing tests',flush=True)


if __name__=='__main__':main()
