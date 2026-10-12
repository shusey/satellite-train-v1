import json
from pathlib import Path
import numpy as np
import pytest
from sat_string.reference_condition import normalize_saved,make_reference,load_frozen
from sat_string.passive_dispersion import initial_profile_state
from sat_string.controller import gap_kinematics

ROOT=Path(__file__).parents[1]


def source():
    return json.loads((ROOT/'docs/dispersion_control_transfer/data/profiles.json').read_text())


def test_saved_shapes_and_exact_sample_std(tmp_path):
    old=source();obj=make_reference(old)
    for name,target,original in [('reference_delta_a_profile',.5,old['profiles'][0]['delta_a_m']),
                                 ('reference_B_profile',.005,old['profiles'][1]['epsilon'])]:
        x=np.asarray(obj[name]);z=np.asarray(original);z-=z.mean()
        assert abs(x.mean())<1e-16
        assert x.std(ddof=1)==pytest.approx(target,rel=1e-15)
        np.testing.assert_allclose(x,target*z/z.std(ddof=1),rtol=1e-15,atol=1e-18)
    path=tmp_path/'frozen.json';path.write_text(json.dumps(obj),encoding='utf-8')
    assert load_frozen(path)==obj


def test_AB_is_same_A_plus_same_B_and_OFF_matches():
    ps=make_reference(source())['profiles']
    assert ps[2]['delta_a_m']==ps[0]['delta_a_m']
    assert ps[2]['epsilon']==ps[1]['epsilon']
    for on,off in zip(ps[:3],ps[3:]):
        assert on['delta_a_m']==off['delta_a_m'] and on['B_multiplier']==off['B_multiplier']
        assert on['control_enabled'] and not off['control_enabled']


def test_changed_frozen_profile_rejected(tmp_path):
    obj=make_reference(source());obj['profiles'][2]['epsilon'][0]+=.001
    p=tmp_path/'bad.json';p.write_text(json.dumps(obj),encoding='utf-8')
    with pytest.raises(ValueError,match='profile mismatch'):load_frozen(p)


@pytest.mark.parametrize('values,target', [([1.]*21,.5),([float('nan')]*21,.5),([0.]*20,.5),([0.]*21,-1)])
def test_invalid_normalization(values,target):
    with pytest.raises(ValueError):normalize_saved(values,target)


def test_reference_fixed_configuration_and_zero_initial_position_error():
    import sys
    import yaml
    sys.path.insert(0,str(ROOT/'scripts'))
    from run_dispersion_control import configuration
    spec=yaml.safe_load((ROOT/'configs/reference_dispersion.yaml').read_text())
    assert (spec['sigma_a_m'],spec['sigma_B_fraction'],spec['std_ddof'])==(.5,.005,1)
    assert spec['cache_seconds']==120
    _,cfg=configuration(30)
    profiles=make_reference(source())['profiles']
    state,_=initial_profile_state(cfg,[p['delta_a_m'] for p in profiles],[p['B_multiplier'] for p in profiles])
    for member in state.reshape(2,6,21).transpose(1,0,2):
        e,_,_=gap_kinematics(*member,reference_radius_m=cfg.derived.reference_radius_m,
            target_spacing_rad=cfg.derived.target_spacing_rad,mu_m3_s2=cfg.constants.mu_m3_s2)
        np.testing.assert_allclose(e,0,rtol=0,atol=1e-8)
    assert cfg.derived.controller_kp_s2_inv==1/21600**2
    assert cfg.derived.controller_kd_s_inv==2/21600
    assert cfg.controller.dwell_time_s==300 and cfg.controller.command_delay_s==0
    assert cfg.simulation.integration_step_s==2 and cfg.simulation.output_step_s==10


def test_reference_RMS_and_thresholds_use_pre_reversal_window():
    from sat_string.config import load_config
    from sat_string.controller import initial_formation
    from sat_string.dispersion_control import derived_histories
    from sat_string.reference_metrics import evaluate_reference
    cfg=load_config(ROOT/'tests/fixtures/regression_24h/config.yaml')
    a,u=initial_formation(21,cfg.derived.reference_radius_m,cfg.derived.target_spacing_rad)
    t=np.arange(0.,31.,10.);aa=np.tile(a,(4,1));uu=np.tile(u,(4,1))
    def history(x):return derived_histories(t,aa.copy(),x,np.full((4,21),.01),
        np.zeros((4,21),dtype='int8'),np.full((4,21),4e-12),cfg)
    on=history(uu.copy());uu[2:,1:]-=100100/cfg.derived.reference_radius_m
    off=history(uu)
    m=evaluate_reference(on,cfg,[],off)
    assert m['matched_pre_allLOW_reversal_improvement']['last_valid_sample_s']==10
    assert m['matched_pre_allLOW_reversal_RMS_improvement']['all_low_m']<1e-8
    assert not m['matched_pre_allLOW_reversal_thresholds'][0]['OFF']['ever_exceeded']
    assert off['gap_error_m'][2,0]<-100000


def test_completed_cycle_rate_counts_completion_time():
    from sat_string.config import load_config
    from sat_string.controller import initial_formation
    from sat_string.dispersion_control import derived_histories
    from sat_string.reference_metrics import evaluate_reference
    cfg=load_config(ROOT/'tests/fixtures/regression_24h/config.yaml')
    a,u=initial_formation(21,cfg.derived.reference_radius_m,cfg.derived.target_spacing_rad)
    t=np.array([0.,10.,20.])
    h=derived_histories(t,np.tile(a,(3,1)),np.tile(u,(3,1)),np.full((3,21),.01),
        np.zeros((3,21),dtype='int8'),np.full((3,21),4e-12),cfg)
    events=[dict(satellite_id=1,command_time=0.,start_time=0.,complete_time=5.,from_mode='L',to_mode='H'),
        dict(satellite_id=1,command_time=10.,start_time=10.,complete_time=20.,from_mode='H',to_mode='L')]
    m=evaluate_reference(h,cfg,events,h)
    assert m['segments'][0]['completed_cycles']==1
    assert m['segments'][0]['cycles_per_day']==pytest.approx(86400/20)
