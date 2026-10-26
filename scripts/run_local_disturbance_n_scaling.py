"""Run six nominal/disturbed pairs; no existing run may be overwritten."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import yaml
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from sat_string.config import load_config, dump_normalized_config
from sat_string.dispersion_control import integrate_control_profiles
from sat_string.n_scaling import CentralNodePulse, scaling_configuration
from sat_string.io import write_json

OUTPUT = ROOT/'outputs/local_disturbance_n_scaling_20261026'
SPEC = ROOT/'configs/local_disturbance_n_scaling.yaml'


def old_artifacts():
    """Size/mtime of all old local artifacts; hash all frozen fixture files."""
    rows = {}
    for base in ['runs', 'outputs']:
        for p in (ROOT/base).rglob('*'):
            if p.is_file() and OUTPUT not in p.parents:
                s = p.stat()
                rows[str(p.relative_to(ROOT))] = [s.st_size, s.st_mtime_ns]
    hashes = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in (ROOT/'tests/fixtures/regression_24h').rglob('*') if p.is_file()}
    return dict(files=rows, fixture_sha256=hashes)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--counts', type=int, nargs='+')
    args = parser.parse_args()
    spec = yaml.safe_load(SPEC.read_text(encoding='utf-8'))
    OUTPUT.mkdir(parents=True, exist_ok=True)
    manifest = OUTPUT/'old_artifacts_before.json'
    if not manifest.exists(): write_json(manifest, old_artifacts())
    base = load_config(ROOT/spec['base_configuration'])
    counts = spec['satellite_counts'] if args.counts is None else args.counts
    for n in counts:
        if n not in spec['satellite_counts']: raise ValueError('N outside preregistered set')
        cfg = scaling_configuration(base, n, spec['duration_days'])
        pulse = CentralNodePulse(n, **{k:v for k,v in spec['input'].items()
                                      if k in ['amplitude_m_s2','start_time_s','duration_s']})
        out = OUTPUT/f'N{n:03d}'
        if out.exists():
            if not (out/'completion.json').exists():
                raise FileExistsError(f'Incomplete run retained: {out}')
            saved=json.loads((out/'metadata.json').read_text(encoding='utf-8'))
            if saved['configuration']!=cfg.to_dict() or saved['node_input']!=pulse.metadata():
                raise ValueError(f'Completed run has different inputs; choose a new output directory: {out}')
            print(f'Reuse completed N={n}', flush=True)
            continue
        expected = base.to_dict(); actual = cfg.to_dict()
        for obj in [expected, actual]:
            obj['formation'].pop('satellite_count')
            obj['simulation'].pop('duration_s')
            obj['disturbance'].pop('enabled')
        assert expected == actual
        profiles = [dict(name=name, control_enabled=True,
                         delta_a_m=[0.]*n, B_multiplier=[1.]*n)
                    for name in ['nominal', 'disturbed']]
        print(f'Start N={n}', flush=True)
        integrate_control_profiles(cfg, profiles, None, out,
            node_acceleration=pulse, node_input_metadata=pulse.metadata(),
            progress=lambda t,end,wall,counts,gap: print(
                f'N={n}: {t/86400:g}/{end/86400:g} d, {wall:.1f} s, transitions={counts}, min gap={gap:.4f} m',
                flush=True))
        dump_normalized_config(cfg, out/'normalized_config.yaml')
        write_json(out/'input_checks.json', dict(passed=True,
            only_changes_from_baseline=['formation.satellite_count','simulation.duration_s','disturbance.enabled'],
            pair_difference='central node acceleration only',
            zero_initial_radius_dispersion=True, zero_B_dispersion=True,
            identical_initial_longitudes=True, unchanged_controller=True,
            input=pulse.metadata(), spec_sha256=hashlib.sha256(SPEC.read_bytes()).hexdigest()))


if __name__ == '__main__': main()
