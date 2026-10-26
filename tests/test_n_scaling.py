import numpy as np
import pytest
from dataclasses import replace
from sat_string.config import load_config
from sat_string.controller import gap_kinematics
from sat_string.dynamics import mean_motion_rad_s, node_acceleration_radius_rate
from sat_string.dispersion_control import integrate_control_profiles, profile_rhs
from sat_string.n_scaling import CentralNodePulse, scaling_configuration, gap_graph_coordinates, response_summary


@pytest.mark.parametrize('n',[11,21,51,101,201,301])
def test_common_conditions_and_graph(n):
    base = load_config('configs/baseline.yaml')
    cfg = scaling_configuration(base,n,14)
    assert cfg.controller == base.controller
    assert cfg.derived.controller_kp_s2_inv == base.derived.controller_kp_s2_inv
    assert cfg.derived.controller_kd_s_inv == base.derived.controller_kd_s_inv
    assert cfg.derived.target_spacing_rad == base.derived.target_spacing_rad
    pulse=CentralNodePulse(n)
    d=pulse(12600).reshape(2,n)
    assert np.count_nonzero(d[0])==0
    assert np.count_nonzero(d[1])==1
    np.testing.assert_allclose(np.diff(d[1])[n//2-1:n//2+1],[1e-6,-1e-6])
    r,side=gap_graph_coordinates(n)
    np.testing.assert_array_equal(r[:n//2][::-1],np.arange(n//2))
    np.testing.assert_array_equal(r[n//2:],np.arange(n//2))
    assert side[n//2-1]=='left' and side[n//2]=='right'


def test_input_units_and_exact_gap_acceleration_mapping():
    a=np.array([6778130.,6778140.,6778137.]); d=np.array([0.,1e-6,0.])
    mu=3.986004418e14; ar=6778137.
    da=node_acceleration_radius_rate(a,d,reference_radius_m=ar,mu_m3_s2=mu)
    dg_da=-1.5*ar*mean_motion_rad_s(mu,a)/a
    np.testing.assert_allclose(dg_da*da,d,atol=1e-22)
    np.testing.assert_allclose(np.diff(dg_da*da),[1e-6,-1e-6])
    assert da[1]<0
    with pytest.raises(ValueError):
        node_acceleration_radius_rate(a,[1.],reference_radius_m=ar,mu_m3_s2=mu)


def test_pulse_integral_and_boundaries():
    pulse=CentralNodePulse(11)
    assert np.count_nonzero(pulse(10800))==np.count_nonzero(pulse(14400))==0
    t=np.linspace(10800,14400,361)
    integral=np.trapezoid([pulse(x)[16] for x in t],t)
    assert integral==pytest.approx(1e-6*3600/2)
    with pytest.raises(ValueError): CentralNodePulse(10)


def test_zero_hook_preserves_old_integrator(tmp_path):
    base=load_config('configs/baseline.yaml')
    cfg=scaling_configuration(base,11,14)
    cfg=replace(cfg,simulation=replace(cfg.simulation,duration_s=2000.))
    p=[dict(control_enabled=True,delta_a_m=[0.]*11,B_multiplier=[1.]*11)]
    integrate_control_profiles(cfg,p,None,tmp_path/'legacy')
    integrate_control_profiles(cfg,p,None,tmp_path/'zero',node_acceleration=lambda t:np.zeros(11),
                               node_input_metadata=dict(boundaries_s=[]))
    for key in ['a_m','lambda_rad','area_m2','mode','density_kg_m3']:
        np.testing.assert_array_equal(np.load(tmp_path/'legacy'/f'{key}.npy'),np.load(tmp_path/'zero'/f'{key}.npy'))


def test_injected_pair_initial_and_early_rate(tmp_path):
    cfg=scaling_configuration(load_config('configs/baseline.yaml'),11,14)
    cfg=replace(cfg,simulation=replace(cfg.simulation,duration_s=20.))
    p=[dict(control_enabled=True,delta_a_m=[0.]*11,B_multiplier=[1.]*11)]*2
    pulse=CentralNodePulse(11,start_time_s=0,duration_s=20,amplitude_m_s2=1e-6)
    integrate_control_profiles(cfg,p,None,tmp_path/'pair',node_acceleration=pulse,node_input_metadata=pulse.metadata())
    a=np.load(tmp_path/'pair/a_m.npy'); u=np.load(tmp_path/'pair/lambda_rad.npy')
    np.testing.assert_array_equal(a[0,0],a[0,1]); np.testing.assert_array_equal(u[0,0],u[0,1])
    _,v,_=gap_kinematics(a[-1],u[-1],reference_radius_m=cfg.derived.reference_radius_m,
        target_spacing_rad=cfg.derived.target_spacing_rad,mu_m3_s2=cfg.constants.mu_m3_s2)
    np.testing.assert_allclose((v[1]-v[0])[4:6],[1e-5,-1e-5],rtol=1e-4,atol=1e-9)


def test_summary_censoring_terminal_and_signed_peak():
    t=np.arange(6.)
    x=np.array([[0,0],[1,-2],[3,-1],[2,-1],[1,0],[0,0.]])
    m=response_summary(t,x,amplitude=2,terminal_window_s=1,pulse_end_s=1,recovery_band_m=1,hold_s=1)
    assert m['P_max_m']==3 and m['gain_s2']==1.5 and m['max_response_gap_id']==1
    assert m['peak_time_by_gap_s']==[2.,1.]
    assert m['E_ult_m']==1 and m['recovery_candidate_s']==4
    m=response_summary(t,t[:,None],amplitude=1,terminal_window_s=1,pulse_end_s=1)
    assert m['peak_at_observation_end'] and m['recovery_candidate_s'] is None


def test_default_batch_matches_original_hybrid_with_switching(tmp_path):
    from sat_string.simulator import simulate_case, CaseDefinition
    base=load_config('configs/baseline.yaml')
    cfg=replace(base,simulation=replace(base.simulation,duration_s=24000.))
    result=simulate_case(cfg,CaseDefinition('original',True,True))
    p=[dict(control_enabled=True,delta_a_m=[0.]*21,B_multiplier=[1.]*21)]
    _,events=integrate_control_profiles(cfg,p,None,tmp_path/'batch',disturbance_enabled=True)
    assert len(events[0])>0
    assert events[0]==[e.to_dict() for e in result.switch_events]
    for key in ['a_m','lambda_rad','density_kg_m3','area_m2','mode']:
        np.testing.assert_array_equal(np.load(tmp_path/'batch'/f'{key}.npy')[:,0],getattr(result,key))
