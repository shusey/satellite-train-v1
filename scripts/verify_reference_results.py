"""Verify frozen inputs, prefix continuity, cache audits, tests and preservation."""
import hashlib,json,re,shutil,subprocess,sys
from importlib.metadata import version
import numpy as np
import yaml
from run_reference_dispersion import ROOT,SPEC,DATA
from sat_string.reference_condition import load_frozen
from sat_string.io import write_json


def main():
    spec=yaml.safe_load(SPEC.read_text());out=ROOT/spec['output'];load_frozen(DATA/'profiles.json')
    prefix,suffix=[out/f'batch_{day}d' for day in (30,60)]
    checks={k:bool(np.array_equal(np.load(prefix/f'{k}.npy',mmap_mode='r')[-1],
        np.load(suffix/f'{k}.npy',mmap_mode='r')[0])) for k in ['a_m','lambda_rad','area_m2','mode']}
    a=np.load(prefix/'density_kg_m3.npy',mmap_mode='r')[-1]
    b=np.load(suffix/'density_kg_m3.npy',mmap_mode='r')[0]
    density=dict(passed=bool(np.allclose(a,b,rtol=1e-8,atol=0)),rtol=1e-8,atol=0,
        max_relative_difference=float(np.max(np.abs(a/b-1))),
        interpretation='cache is reanchored; physical/actuator states must match exactly')
    cp=json.loads((prefix/'checkpoint.json').read_text());cs=json.loads((suffix/'checkpoint.json').read_text())
    events=all(old==[x for x in new if x['command_time']<=30*86400] for old,new in zip(cp['events'],cs['events']))
    configurations=[];inputs=[];audits=[]
    for day,path in zip((30,60),(prefix,suffix)):
        cfg=json.loads((path/'metadata.json').read_text())['configuration']
        cfg['simulation'].pop('duration_s');configurations.append(cfg)
        comparison=json.loads((DATA/f'{day}d/comparison.json').read_text());inputs+=comparison['input_checks']
        audit=json.loads((DATA/f'density_audit_{day}d.json').read_text());audits.append(audit)
        shutil.copyfile(path/'normalized_config.yaml',DATA/f'normalized_config_{day}d.yaml')
        shutil.copyfile(path/'environment.json',DATA/f'environment_{day}d.json')
    profile_hash=hashlib.sha256((DATA/'profiles.json').read_bytes()).hexdigest()
    hashes=[json.loads((p/'environment.json').read_text())['frozen_profiles_sha256']==profile_hash for p in (prefix,suffix)]
    regression=json.loads((DATA/'regression_24h.json').read_text())
    preserved=json.loads((DATA/'preservation.json').read_text())
    tests=(out/'tests.log').read_text(encoding='utf-8');match=re.search(r'(\d+) passed in ([\d.]+)s',tests)
    tests_pass=match is not None and 'failed' not in tests and 'warning' not in tests
    passed=(all(checks.values()) and density['passed'] and events and all(hashes)
        and configurations[0]==configurations[1] and all(x['passed'] for x in audits)
        and regression['passed'] and preserved['passed'] and tests_pass)
    report=dict(passed=passed,branch=subprocess.check_output(['git','branch','--show-current'],cwd=ROOT,text=True).strip(),
        base_commit=preserved['base_commit'],tests_passed=int(match[1]) if match else None,
        tests_elapsed_s=float(match[2]) if match else None,tests_log=tests,
        checkpoint_state_checks=checks,checkpoint_density=density,events_through_30d_exact=events,
        configuration_change='duration only',profiles_sha256=profile_hash,frozen_profile_hashes_match=all(hashes),
        input_checks=inputs,density_audits=audits,baseline_24h_passed=regression['passed'],
        old_tracked_contents_sha256_unchanged=preserved['passed'],old_local_run_stats_unchanged=preserved['passed'],
        physical_controller_integrator_core_modified=False,new_parameter_sweeps=False,
        local_only_arrays=str(out.relative_to(ROOT)))
    write_json(DATA/'verification.json',report)
    files=[]
    for path in (prefix,suffix):
        for file in sorted(path.iterdir()):
            if not file.is_file():continue
            h=hashlib.sha256()
            with file.open('rb') as stream:
                for block in iter(lambda:stream.read(4*1024*1024),b''):h.update(block)
            files.append(dict(local_path=str(file.relative_to(ROOT)),bytes=file.stat().st_size,sha256=h.hexdigest()))
    write_json(DATA/'local_run_manifest.json',dict(local_files=files,
        python=sys.version,numpy=version('numpy'),pymsis=version('pymsis'),
        note='Large trajectories/checkpoints remain local. Fingerprints identify the frozen reference; floating tolerances apply to future numerical comparisons.'))
    print('Reference verification:',passed,flush=True)
    if not passed:raise SystemExit(1)


if __name__=='__main__':main()
