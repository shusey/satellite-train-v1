from pathlib import Path
import runpy
import numpy as np
import pytest
from sat_string.config import load_config,config_from_dict,dump_normalized_config
from sat_string.simulator import simulate_case,STANDARD_CASES
from sat_string.kp_analysis import activation_state,enhance_activations,match_all_low,fixed_state_forecasts,off_before_dwell
from sat_string.long_horizon import decompose_q


def baseline():return load_config('tests/fixtures/regression_24h/config.yaml')


def test_legacy_Kp_omitted_and_derived_unchanged():
    c=baseline();assert c.derived.controller_kp_s2_inv==1/21600.**2
    assert 'kp_s2_inv' not in c.to_dict()['controller']
    assert config_from_dict(c.to_dict())==c


@pytest.mark.parametrize('factor',[.3,1.,3.])
def test_independent_kp_roundtrip_only_requested_changes(tmp_path,factor):
    c=baseline();fn=runpy.run_path(str(Path(__file__).parents[1]/'scripts/run_kp_case.py'))['kp_config']
    new=fn(c,factor,30);expected=c.to_dict()
    expected['controller']['kp_s2_inv']=factor*c.derived.controller_kp_s2_inv
    expected['derived']['controller_kp_s2_inv']=expected['controller']['kp_s2_inv']
    expected['simulation']['duration_s']=30*86400.
    assert new.to_dict()==expected
    assert new.derived.controller_kd_s_inv==c.derived.controller_kd_s_inv
    dump_normalized_config(new,tmp_path/'config.yaml');assert load_config(tmp_path/'config.yaml')==new


@pytest.mark.parametrize('value',[0.,-1.,float('nan'),float('inf'),True])
def test_invalid_kp_rejected(value):
    data=baseline().to_dict();data['controller']['kp_s2_inv']=value
    with pytest.raises(ValueError):config_from_dict(data)


def test_explicit_default_kp_preserves_short_hybrid(short_config):
    data=short_config.to_dict();data['controller']['kp_s2_inv']=short_config.derived.controller_kp_s2_inv
    a=simulate_case(short_config,STANDARD_CASES[0]);b=simulate_case(config_from_dict(data),STANDARD_CASES[0])
    assert [x.to_dict() for x in a.switch_events]==[x.to_dict() for x in b.switch_events]
    for key,value in a.arrays().items():np.testing.assert_array_equal(value,b.arrays()[key])


@pytest.mark.parametrize('satellite,ids',[ (1,[1]),(2,[1,2]),(21,[20]) ])
def test_activation_uses_physical_error_and_incident_gaps(satellite,ids):
    c=baseline();e=np.arange(1,21,dtype=float);v=e*.0001
    r=activation_state(e,v,satellite,c)
    assert r['incident_gap_ids']==ids
    assert r['incident_signed_error_m']==[float(k) for k in ids]
    assert r['incident_max_abs_error_m']==max(ids)
    assert r['max_abs_error_m']==20
    assert r['rms_error_m']==pytest.approx(np.sqrt(np.mean(e*e)))
    assert r['q_P_rate_m_s3']==pytest.approx(c.derived.controller_kp_s2_inv*r['B_rate_m_s'])


@pytest.mark.parametrize('satellite',[0,22])
def test_activation_invalid_satellite(satellite):
    with pytest.raises(ValueError):activation_state(np.zeros(20),np.zeros(20),satellite,baseline())


def test_activation_preserves_signed_state_and_separates_command_from_crossing():
    c=baseline();e=np.array([np.arange(20,dtype=float),-np.arange(20,dtype=float)])
    v=np.ones_like(e)*.0001;t=np.array([10.,20.])
    rows=[dict(satellite_id=2,command_time_s=20.,sampled_crossing_time_s=10.)]
    r=enhance_activations(rows,dict(time_s=t,gap_error_m=e,gap_rate_m_s=v),c)[0]
    assert r['crossing_incident_signed_error_m']==[0.,1.]
    assert r['command_incident_signed_error_m']==[0.,-1.]


def test_same_state_Kp_scales_P_and_slope_not_D_or_physical_error():
    c=baseline();e=np.arange(20,dtype=float);v=np.full(20,.0001)
    p0,d0=decompose_q(e[None,:],v[None,:],c.derived.controller_kp_s2_inv,c.derived.controller_kd_s_inv)
    for factor in [.3,3.]:
        data=c.to_dict();data['controller']['kp_s2_inv']=factor*c.derived.controller_kp_s2_inv
        new=config_from_dict(data);p,d=decompose_q(e[None,:],v[None,:],new.derived.controller_kp_s2_inv,new.derived.controller_kd_s_inv)
        np.testing.assert_allclose(p,p0*factor);np.testing.assert_array_equal(d,d0)
        r=activation_state(e,v,1,new);base=activation_state(e,v,1,c)
        assert r['q_P_rate_m_s3']==pytest.approx(factor*base['q_P_rate_m_s3'])
        assert r['max_abs_error_m']==base['max_abs_error_m']


def test_same_state_Kp_forecast_changes_both_offset_and_slope():
    c=baseline();e=np.zeros(20);v=np.zeros(20);e[0]=1;v[0]=.0001
    rows=fixed_state_forecasts(e,v,0.,np.zeros(21),c)
    assert all(x['satellite_id']==1 for x in rows)
    assert rows[0]['predicted_command_s']>rows[1]['predicted_command_s']>rows[2]['predicted_command_s']
    assert rows[2]['q_P_rate_m_s3']==pytest.approx(10*rows[0]['q_P_rate_m_s3'])
    assert len({x['q_D_m_s2'] for x in rows})==1


def test_all_LOW_match_allows_Kp_but_rejects_Kd():
    c=baseline();data=c.to_dict();data['controller']['kp_s2_inv']=3*c.derived.controller_kp_s2_inv
    assert match_all_low(config_from_dict(data),c)
    data['controller']['kd_s_inv']=3*c.derived.controller_kd_s_inv
    assert not match_all_low(config_from_dict(data),c)


def test_early_off_sample_not_same_as_eligible_return_gate():
    c=baseline();t=np.arange(0,601,10.);q=np.full((len(t),21),2*c.controller.threshold_off_m_s2)
    q[2,0]=0
    event=dict(satellite_id=1,command_time=10.,start_time=10.,complete_time=130.,to_mode='H')
    row=off_before_dwell(dict(time_s=t,q_m_s2=q),c,[event])[0]
    assert row['first_off_after_HIGH_start_s']==20 and row['off_observed_before_gate']
    assert row['earliest_return_gate_s']==430 and row['gate_observed']
    assert row['first_off_to_LOW_command_s'] is None


def test_off_window_right_censored():
    c=baseline();t=np.arange(0,201,10.)
    event=dict(satellite_id=1,command_time=10.,start_time=10.,complete_time=130.,to_mode='H')
    row=off_before_dwell(dict(time_s=t,q_m_s2=np.zeros((len(t),21))),c,[event])[0]
    assert row['off_observed_before_gate'] and not row['gate_observed']


def test_first_off_after_gate_is_recorded_with_actual_return_delay():
    c=baseline();t=np.arange(0,601,10.);q=np.full((len(t),21),2*c.controller.threshold_off_m_s2)
    q[t>=450,0]=0
    events=[dict(satellite_id=1,command_time=10.,start_time=10.,complete_time=130.,to_mode='H'),
        dict(satellite_id=1,command_time=450.,start_time=450.,complete_time=570.,to_mode='L')]
    row=off_before_dwell(dict(time_s=t,q_m_s2=q),c,events)[0]
    assert row['first_off_after_HIGH_start_s']==450
    assert not row['off_observed_before_gate'] and row['gate_observed']
    assert row['first_off_to_LOW_command_s']==0
