"""Continue ONLY N=21/301 from saved day14 to day30, never repeat day0-14."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'scripts'))
from run_local_disturbance_n_scaling import OUTPUT as SOURCE,old_artifacts
from sat_string.n_scaling import thirty_day_extension
from sat_string.config import dump_normalized_config,load_config
from sat_string.dispersion_control import integrate_control_profiles
from sat_string.io import write_json

OUTPUT=ROOT/'outputs/local_disturbance_n_scaling_30d_20261026'


def main():
    p=argparse.ArgumentParser();p.add_argument('--count',type=int,choices=[21,301],required=True)
    args=p.parse_args();n=args.count
    OUTPUT.mkdir(parents=True,exist_ok=True)
    manifest=OUTPUT/'old_artifacts_before.json'
    # Written once before launching the two independent N processes.
    if not manifest.exists():write_json(manifest,old_artifacts(OUTPUT))
    source=SOURCE/f'N{n:03d}';out=OUTPUT/f'N{n:03d}'
    completion=json.loads((source/'completion.json').read_text(encoding='utf-8'))
    assert completion['complete'] and completion['end_s']==14*86400
    cp=json.loads((source/'checkpoint.json').read_text(encoding='utf-8'))
    cfg,pulse=thirty_day_extension(cp)
    assert cp['configuration']==load_config(source/'normalized_config.yaml').to_dict()
    for key,part in [('a_m',0),('lambda_rad',1)]:
        np.testing.assert_array_equal(np.load(source/f'{key}.npy',mmap_mode='r')[-1].ravel(),np.asarray(cp['state'])[part])
    if out.exists():
        saved=json.loads((out/'metadata.json').read_text(encoding='utf-8'))
        if not (out/'completion.json').exists() or saved['configuration']!=cfg.to_dict() or saved['node_input']!=pulse.metadata():
            raise FileExistsError('Existing mismatched or incomplete output retained')
        print(f'Reuse completed N={n} day14-30',flush=True);return
    integrate_control_profiles(cfg,cp['profiles'],None,out,checkpoint=cp,
        node_acceleration=pulse,node_input_metadata=pulse.metadata(),
        progress=lambda t,end,wall,counts,gap:print(
            f'N={n} {t/86400:g}/{end/86400:g}d; {wall:.1f}s; events={counts}; min gap={gap:.4f}m',flush=True))
    dump_normalized_config(cfg,out/'normalized_config.yaml')
    write_json(out/'input_checks.json',dict(passed=True,source=str(source.relative_to(ROOT)),
        source_checkpoint_sha256=hashlib.sha256((source/'checkpoint.json').read_bytes()).hexdigest(),
        only_change='simulation.duration_s: 14d -> 30d',fresh_physics_interval_days=[14,30],
        unchanged_controller=True,unchanged_pulse=True,pulse_replayed=False,initial_state_from_checkpoint=True))


if __name__=='__main__':main()
