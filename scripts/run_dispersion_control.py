"""Fixed relay controller versus frozen dispersion profiles; no parameter tuning."""
import argparse, hashlib, json, os, sys
from pathlib import Path
import numpy as np
import yaml
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT/'src'))
from sat_string.config import load_config, config_from_dict, dump_normalized_config
from sat_string.io import write_json
from sat_string.dispersion_control import integrate_control_profiles
from run_nrlmsis_case import make_environment

OUT = Path(os.environ.get('SAT_CONTROL_OUTPUT',ROOT/'outputs/dispersion_control_20261011'))
DOC = Path(os.environ.get('SAT_CONTROL_DOC',ROOT/'docs/dispersion_control_transfer')); DATA = DOC/'data'
SPEC = ROOT/'configs/dispersion_control_transfer.yaml'


def frozen_profiles(spec):
    previous = json.loads((ROOT/spec['screening_profiles']).read_text())['profiles']
    metrics = json.loads((ROOT/'docs/initial_dispersion_screening/data/metrics.json').read_text())['results']
    result = []
    for case, factor, sigma in [('A', 'a', spec['sigma_a_m']), ('B', 'B', spec['sigma_B_fraction'])]:
        candidates = [r for r in metrics if r['factor']==factor and r['nominal_sigma']==sigma
                      and r['kind']=='random' and r['days']==7]
        median = sorted(candidates, key=lambda r:r['max_error_m'])[len(candidates)//2]
        p = next(p for p in previous if p['id']==median['case_id'])
        result.append(dict(p, case=case, control_enabled=True,
                           source_profile=p['id'], selection='median7d', screening_7d_max_m=median['max_error_m']))
    old = json.loads((ROOT/spec['combination_profile']).read_text())['profiles'][0]
    off = np.asarray(old['delta_a_m'])*spec['sigma_a_m']
    eps = np.asarray(old['epsilon'])*(spec['sigma_B_fraction']/.01)
    ab = dict(old, id='AB_a0.1m_B0.3percent_seed0', case='AB', control_enabled=True,
        source_profile=old['id'], selection='scale previously frozen combination, no new draws',
        delta_a_m=off.tolist(), epsilon=eps.tolist(), B_multiplier=(1+eps).tolist(),
        nominal_sigma=spec['sigma_a_m'], realized_sigma_a_m=float(off.std()),
        realized_sigma_B_fraction=float(eps.std()))
    return result+[ab, dict(ab, case='AB_OFF', control_enabled=False)]


def configuration(days):
    spec = yaml.safe_load(SPEC.read_text())
    path=ROOT/spec['baseline_run']/'normalized_config.yaml'
    if not path.exists():
        # A clean Git checkout has the frozen configuration but no huge local runs.
        archive=ROOT/'docs/dispersion_control_transfer/data'
        path=archive/('normalized_config_60d.yaml' if (archive/'normalized_config_60d.yaml').exists() else 'normalized_config_30d.yaml')
    base = load_config(path)
    raw = base.to_dict(); raw['simulation']['duration_s'] = float(days)*86400
    cfg = config_from_dict(raw)
    checked = cfg.to_dict(); checked['simulation']['duration_s'] = base.simulation.duration_s
    if checked != base.to_dict(): raise ValueError('Only duration may change')
    if not (cfg.derived.controller_kp_s2_inv==1/21600**2 and cfg.derived.controller_kd_s_inv==2/21600
            and cfg.controller.dwell_time_s==300 and cfg.controller.slew_time_s==120
            and cfg.controller.control_period_s==10 and cfg.controller.threshold_on_m_s2==1.5e-7
            and cfg.controller.threshold_off_m_s2==1.4e-7 and not cfg.disturbance.enabled):
        raise ValueError('Unexpected controller or atmosphere changes')
    return spec, cfg


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--days', type=int, choices=[30,60], required=True)
    parser.add_argument('--prefix', type=Path)
    args = parser.parse_args(); spec, cfg = configuration(args.days)
    DATA.mkdir(parents=True, exist_ok=True)
    profiles = frozen_profiles(spec)
    checkpoint = None
    if args.prefix:
        if not (args.prefix/'completion.json').exists(): raise ValueError('Prefix incomplete')
        checkpoint = json.loads((args.prefix/'checkpoint.json').read_text())
        provenance = json.loads((args.prefix/'environment.json').read_text())
    if args.days==60 and not args.prefix: raise ValueError('60d requires saved 30d checkpoint and go decision')
    if args.days==60:
        decision = json.loads((DATA/'extension_decision.json').read_text())
        if not decision['extend_to_60d'] or not decision.get('reviewed',False): raise ValueError('No reviewed scientific go decision')
    env, backend, provider = make_environment(ROOT/spec['environment'], ROOT/spec['weather'], spec['cache_seconds'])
    weather_hash = hashlib.sha256((ROOT/spec['weather']).read_bytes()).hexdigest()
    prior = json.loads((ROOT/spec['baseline_run']/'density_environment.json').read_text())
    if prior['environment']!=env or prior['weather_sha256']!=weather_hash or prior['cache_seconds']!=spec['cache_seconds']:
        raise ValueError('Environment differs from frozen active period')
    if checkpoint and (provenance['environment']!=env or provenance['weather_sha256']!=weather_hash
                       or provenance['cache_seconds']!=spec['cache_seconds']):
        raise ValueError('Continuation environment differs')
    if not (DATA/'profiles.json').exists(): write_json(DATA/'profiles.json', dict(profiles=profiles, specification=spec))
    elif json.loads((DATA/'profiles.json').read_text())['profiles']!=profiles: raise ValueError('Frozen profiles changed')
    output = OUT/f'batch_{args.days}d'
    callback = lambda t, end, wall, counts, gap: print(
        f'{t/86400:g}/{end/86400:g}d; {wall:.1f}s; commands A/B/AB/OFF={counts}; min gap={gap:.3f}m', flush=True)
    integrate_control_profiles(cfg, profiles, provider, output, checkpoint=checkpoint, progress=callback)
    dump_normalized_config(cfg, output/'normalized_config.yaml')
    write_json(output/'environment.json', dict(environment=env, weather_sha256=weather_hash,
        cache_seconds=spec['cache_seconds'], pymsis_version='0.12.0', msis_version='2.0',
        density_units='kg/m3', backend_points=backend.points_evaluated,
        maximum_cache_radius_offset_m=provider.maximum_radius_offset_m,
        maximum_cache_phase_offset_rad=provider.maximum_phase_offset_rad,
        prefix=str(args.prefix) if args.prefix else None))
    print('Saved', output, flush=True)


if __name__=='__main__': main()
