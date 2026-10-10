from dataclasses import replace
import numpy as np
import pytest
from sat_string.config import load_config
from sat_string.controller import switching_function
from sat_string.nrlmsis_analysis import pd_components,inactive_prefix_proof,density_statistics,density_to_gap_acceleration,current_radius_arc_error


def test_batched_pd_matches_actual_controller():
    e=np.array([[1.,-2.,3.],[-4.,5.,6.]])
    v=e*.02;p,d=pd_components(e,v,.3,.2)
    for i in range(2):np.testing.assert_allclose(p[i]+d[i],switching_function(e[i],v[i],.3,.2))


def test_inactive_prefix_checks_every_tick_and_rejects_wrong_schedule():
    c=load_config('tests/fixtures/regression_24h/config.yaml')
    t=np.arange(0,1001,10.);e=np.zeros((len(t),20));v=e.copy()
    e[t>=500,-1]=-100
    q=np.asarray([switching_function(x,y,c.derived.controller_kp_s2_inv,c.derived.controller_kd_s_inv) for x,y in zip(e,v)])
    arrays=dict(time_s=t,gap_error_m=e,gap_rate_m_s=v,q_m_s2=q,mode=np.zeros_like(q),area_m2=np.full_like(q,.01))
    proof=inactive_prefix_proof(arrays,c,120)
    assert proof['first_potential_high_time_s']==500
    assert proof['reusable_until_s']==360
    assert proof['all_ticks_checked']
    with pytest.raises(ValueError):inactive_prefix_proof(arrays,replace(c,simulation=replace(c.simulation,output_step_s=20)),120)
    arrays['q_m_s2'][0,0]=1
    with pytest.raises(AssertionError):inactive_prefix_proof(arrays,c,120)


def test_density_statistics_distinguish_common_and_adjacent():
    d=density_statistics(np.array([[1.,2.,3.],[3.,3.,3.]])*1e-12)
    assert d['maximum_fleet_spread_over_mean']==pytest.approx(1.)
    assert d['maximum_adjacent_absolute_difference_kg_m3']==pytest.approx(1e-12)
    assert density_statistics(np.ones((2,3))*1e-12)['rms_differential_density_kg_m3']==0


def test_gap_acceleration_has_correct_sign_and_drag_scaling():
    c=load_config('tests/fixtures/regression_24h/config.yaml');a=np.full(3,c.derived.reference_radius_m)
    rho=np.array([1.,2.,3.])*1e-12;area=np.full(3,.01)
    acc=density_to_gap_acceleration(a,rho,area,c)
    assert np.all(acc>0)
    np.testing.assert_allclose(density_to_gap_acceleration(a,2*rho,area,c),2*acc)
    np.testing.assert_allclose(density_to_gap_acceleration(a,np.ones(3)*1e-12,area,c),0)


def test_common_radius_contraction_is_not_reference_radius_formation_error():
    a0=6778137.;phase=np.array([0.,100000/a0,200000/a0])
    np.testing.assert_allclose(current_radius_arc_error(np.full(3,a0),phase,100000),0,atol=1e-10)
    np.testing.assert_allclose(current_radius_arc_error(np.full(3,.999*a0),phase,100000),-100,atol=1e-9)
