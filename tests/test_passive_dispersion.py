from pathlib import Path
import numpy as np
import pytest
from sat_string.config import config_from_dict,load_config
from sat_string.controller import gap_kinematics
from sat_string.passive_dispersion import (centered_normal,centered_alternating,
    initial_profile_state,validate_profiles,passive_rhs,theory_initial_gap_rate,integrate_profiles,state_metrics)
from sat_string.integrator import rk4_step

BASE=load_config(Path(__file__).parents[1]/'tests/fixtures/regression_24h/config.yaml')

def configuration(duration=100.):
    raw=BASE.to_dict();raw['atmosphere']['density_model']='nrlmsis2';raw['disturbance']['enabled']=False
    raw['disturbance']['start_time_s']=0.;raw['disturbance']['duration_s']=1.
    raw['simulation']['duration_s']=duration
    return config_from_dict(raw)

def test_centered_normal_reproducible_no_variance_renormalization():
    expected=np.random.default_rng(3).normal(0,1,21);expected-=expected.mean()
    np.testing.assert_array_equal(centered_normal(21,1,3),expected)
    np.testing.assert_allclose(centered_normal(21,.1,3),.1*expected,atol=1e-16)
    assert abs(expected.mean())<1e-15

def test_odd_alternating_preserves_mean_and_adjacent_difference():
    x=centered_alternating(21,1)
    assert abs(x.mean())<1e-15
    np.testing.assert_allclose(np.abs(np.diff(x)),2)

@pytest.mark.parametrize('factors',[[0]*21,[1.01]*21,[float('nan')]*21])
def test_invalid_B_profile_rejected(factors):
    with pytest.raises(ValueError):validate_profiles(np.zeros(21),factors,21)

def test_initial_radius_changes_but_position_does_not():
    cfg=configuration();off=centered_normal(21,1,0)
    state,_=initial_profile_state(cfg,off,np.ones(21))
    e,v,_=gap_kinematics(*state,reference_radius_m=cfg.derived.reference_radius_m,
        target_spacing_rad=cfg.derived.target_spacing_rad,mu_m3_s2=cfg.constants.mu_m3_s2)
    np.testing.assert_allclose(e,0,atol=1e-9)
    np.testing.assert_allclose(v,theory_initial_gap_rate(cfg,off),rtol=2e-6,atol=1e-11)

def test_B_only_constant_density_matches_exact_radius_solution():
    cfg=configuration();epsilon=centered_alternating(21,.01);state,f=initial_profile_state(cfg,np.zeros(21),1+epsilon)
    rho=4e-12;provider=lambda t,a,u:np.full(a.shape,rho)
    y=state.copy()
    for t in range(0,100,2):y=rk4_step(lambda tt,ss:passive_rhs(tt,ss,cfg,provider,f),t,y,2)
    B=cfg.spacecraft.drag_coefficient*cfg.spacecraft.area_low_m2/cfg.spacecraft.mass_kg
    exact=(np.sqrt(state[0])-.5*B*f*np.sqrt(cfg.constants.mu_m3_s2)*rho*100)**2
    np.testing.assert_allclose(y[0],exact,atol=1e-7,rtol=0)
    assert np.diff(y[0])[0]>0  # Higher first-satellite B loses more a.

def test_batch_is_independent_and_zero_profile_matches_original(tmp_path):
    from sat_string.simulator import simulate_case,CaseDefinition
    cfg=configuration();offset=centered_normal(21,1,1)
    profiles=[dict(delta_a_m=[0.]*21,B_multiplier=[1.]*21),dict(delta_a_m=offset.tolist(),B_multiplier=[1.]*21)]
    provider=lambda t,a,u:np.full(a.shape,4e-12)*(1+.1*np.sin(t/100))
    integrate_profiles(cfg,profiles,provider,tmp_path/'batch')
    original=simulate_case(cfg,CaseDefinition('all_low',False,False),density_provider=provider)
    for name in ('a_m','lambda_rad'):
        np.testing.assert_array_equal(np.load(tmp_path/'batch'/f'{name}.npy')[:,0],getattr(original,name))
    integrate_profiles(cfg,[profiles[1]],provider,tmp_path/'single')
    np.testing.assert_array_equal(np.load(tmp_path/'batch'/'a_m.npy')[:,1],np.load(tmp_path/'single'/'a_m.npy')[:,0])

def test_continuation_does_not_repeat_or_change_prefix(tmp_path):
    provider=lambda t,a,u:np.full(a.shape,4e-12)
    p=[dict(delta_a_m=[0.]*21,B_multiplier=[1.]*21)]
    state=integrate_profiles(configuration(50),p,provider,tmp_path/'prefix')
    integrate_profiles(configuration(100),p,provider,tmp_path/'suffix',start_s=50,initial_state=state)
    integrate_profiles(configuration(100),p,provider,tmp_path/'whole')
    np.testing.assert_array_equal(np.load(tmp_path/'suffix'/'a_m.npy'),np.load(tmp_path/'whole'/'a_m.npy')[5:])


def test_saved_state_metrics_rms_crossing_and_gap_distribution():
    cfg=configuration();state,_=initial_profile_state(cfg,np.zeros(21),np.ones(21))
    a=np.tile(state[0],(3,1));u=np.tile(state[1],(3,1));t=np.array([0.,10.,20.])
    u[1,1:]+=110/cfg.derived.reference_radius_m
    u[2,1:]+=600/cfg.derived.reference_radius_m
    result=state_metrics(t,a,u,cfg)
    assert result['max_error_m']==pytest.approx(600,abs=1e-8)
    assert result['rms_error_m']==pytest.approx(np.sqrt((110**2+600**2)/(3*20)),abs=1e-8)
    assert result['gap_rms_m'][0]==pytest.approx(np.sqrt((110**2+600**2)/3),abs=1e-8)
    assert result['first_exceedance_s']=={'100':10.,'500':20.,'1000':None}
    assert result['first_label_order_reversal_s'] is None


def test_label_order_reversal_is_reported_without_reordering():
    cfg=configuration();state,_=initial_profile_state(cfg,np.zeros(21),np.ones(21))
    a=np.tile(state[0],(2,1));u=np.tile(state[1],(2,1))
    u[1,1:]-=100100/cfg.derived.reference_radius_m
    result=state_metrics(np.array([0.,10.]),a,u,cfg)
    assert result['first_label_order_reversal_s']==10
    assert result['min_labelled_spacing_m']==pytest.approx(-100,abs=1e-8)
