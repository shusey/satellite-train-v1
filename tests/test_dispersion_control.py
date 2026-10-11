from pathlib import Path
import json
import numpy as np
import pytest
from sat_string.config import load_config, config_from_dict
from sat_string.attitude import AttitudeStateMachine, H, L
from sat_string.dispersion_control import integrate_control_profiles, profile_rhs, derived_histories
from sat_string.controller import gap_kinematics
from sat_string.passive_dispersion import initial_profile_state
from sat_string.simulator import simulate_case, CaseDefinition

BASE = load_config(Path(__file__).parents[1]/'tests/fixtures/regression_24h/config.yaml')


def cfg(duration=1000):
    raw = BASE.to_dict(); raw['simulation']['duration_s']=duration
    raw['atmosphere']['density_model']='nrlmsis2'; raw['disturbance']['enabled']=False
    raw['disturbance']['start_time_s']=0.;raw['disturbance']['duration_s']=1.
    return config_from_dict(raw)


def zero(on=True):
    return dict(delta_a_m=[0.]*21, B_multiplier=[1.]*21, control_enabled=on)


def provider(t, a, u):
    return np.full(a.shape, 4e-12)*(1+.3*np.exp(-(u-.08)**2/.001))


@pytest.mark.parametrize('on', [False, True])
def test_zero_profile_preserves_original_states_areas_modes_events(tmp_path, on):
    c=cfg(); _, events=integrate_control_profiles(c,[zero(on)],provider,tmp_path/'new')
    old=simulate_case(c,CaseDefinition('test',on,False),density_provider=provider)
    for key in ['a_m','lambda_rad','area_m2','mode','density_kg_m3']:
        np.testing.assert_array_equal(np.load(tmp_path/'new'/f'{key}.npy')[:,0],getattr(old,key))
    assert events[0]==[e.to_dict() for e in old.switch_events]


@pytest.mark.parametrize('time', [0.,60.,120.])
def test_B_multiplier_applies_during_LOW_slew_HIGH(time):
    c=cfg(); state,f=initial_profile_state(c,np.zeros(21),np.ones(21))
    act=AttitudeStateMachine(21,.01,.03,120,300)
    act.apply_commands(0,np.full(21,H));act.process_completions(time)
    f=np.linspace(.99,1.01,21)
    base=profile_rhs(time,state,c,provider,np.ones(21),act.areas)
    changed=profile_rhs(time,state,c,provider,f,act.areas)
    np.testing.assert_allclose(changed[0],base[0]*f,rtol=1e-15)
    np.testing.assert_array_equal(changed[1],base[1])


@pytest.mark.parametrize('boundary', [60,130,550])
def test_checkpoint_preserves_active_slew_and_dwell(tmp_path, boundary):
    c=cfg(boundary); p=zero()
    # Trigger immediately with a position-dependent density/velocity initial state.
    p['delta_a_m']=(np.array([-.5,.5]*10+[0.])).tolist()
    integrate_control_profiles(c,[p],provider,tmp_path/'prefix')
    cp=json.loads((tmp_path/'prefix'/'checkpoint.json').read_text())
    integrate_control_profiles(cfg(2000),[p],provider,tmp_path/'suffix',checkpoint=cp)
    _, events=integrate_control_profiles(cfg(2000),[p],provider,tmp_path/'whole')
    for key in ['a_m','lambda_rad','area_m2','mode']:
        np.testing.assert_array_equal(np.load(tmp_path/'suffix'/f'{key}.npy'),np.load(tmp_path/'whole'/f'{key}.npy')[boundary//10:])
    cp2=json.loads((tmp_path/'suffix'/'checkpoint.json').read_text())
    assert cp2['events']==events
    assert cp['events'][0]  # This is an actuator continuation, not an empty all-LOW case.


def test_checkpoint_rejects_changed_profile(tmp_path):
    p=zero();integrate_control_profiles(cfg(100),[p],provider,tmp_path/'first')
    cp=json.loads((tmp_path/'first'/'checkpoint.json').read_text())
    changed=dict(p,control_enabled=False)
    with pytest.raises(ValueError,match='inputs differ'):
        integrate_control_profiles(cfg(200),[changed],provider,tmp_path/'second',checkpoint=cp)


def test_frozen_profile_selection_and_AB_scaling():
    import sys
    sys.path.insert(0,str(Path(__file__).parents[1]/'scripts'))
    from run_dispersion_control import configuration, frozen_profiles
    spec,c=configuration(30); profiles=frozen_profiles(spec)
    assert [p['seed'] for p in profiles]==[4,2,0,0]
    source=json.loads((Path(__file__).parents[1]/spec['combination_profile']).read_text())['profiles'][0]
    np.testing.assert_array_equal(profiles[2]['delta_a_m'],np.asarray(source['delta_a_m'])*.1)
    np.testing.assert_array_equal(profiles[2]['epsilon'],np.asarray(source['epsilon'])*(.003/.01))
    assert c.derived.controller_kp_s2_inv==1/21600**2


def test_batch_controllers_do_not_couple_profiles(tmp_path):
    p=zero();p['delta_a_m']=np.linspace(-.1,.1,21).tolist()
    integrate_control_profiles(cfg(300),[zero(),p],provider,tmp_path/'batch')
    integrate_control_profiles(cfg(300),[p],provider,tmp_path/'single')
    np.testing.assert_array_equal(np.load(tmp_path/'batch'/'a_m.npy')[:,1],np.load(tmp_path/'single'/'a_m.npy')[:,0])
