from pathlib import Path
import runpy
import numpy as np
import pytest
from sat_string.config import load_config,config_from_dict,dump_normalized_config
from sat_string.simulator import simulate_case,STANDARD_CASES
from sat_string.kd_analysis import velocity_distribution,contribution,match_all_low,threshold_events,return_gates,fixed_state_forecasts,enhance_residual,first_day_metrics
from sat_string.long_horizon import decompose_q


def test_legacy_gain_and_normalization_are_preserved():
    cfg=load_config("tests/fixtures/regression_24h/config.yaml")
    assert cfg.derived.controller_kp_s2_inv==1/21600.**2
    assert cfg.derived.controller_kd_s_inv==2/21600.
    assert "kd_s_inv" not in cfg.to_dict()["controller"]
    assert config_from_dict(cfg.to_dict())==cfg


@pytest.mark.parametrize("factor",[.3,1.,3.])
def test_independent_kd_roundtrip_and_only_requested_changes(tmp_path,factor):
    cfg=load_config("tests/fixtures/regression_24h/config.yaml")
    builder=runpy.run_path(str(Path(__file__).parents[1]/"scripts/run_kd_case.py"))["kd_config"]
    result=builder(cfg,factor,30)
    expected=cfg.to_dict();expected["controller"]["kd_s_inv"]=factor*cfg.derived.controller_kd_s_inv
    expected["derived"]["controller_kd_s_inv"]=expected["controller"]["kd_s_inv"]
    expected["simulation"]["duration_s"]=30*86400.
    assert result.to_dict()==expected
    assert result.derived.controller_kp_s2_inv==cfg.derived.controller_kp_s2_inv
    dump_normalized_config(result,tmp_path/"config.yaml")
    assert load_config(tmp_path/"config.yaml")==result


@pytest.mark.parametrize("value",[0.,-1.,float("nan"),float("inf"),True])
def test_invalid_kd_rejected(value):
    data=load_config("tests/fixtures/regression_24h/config.yaml").to_dict()
    data["controller"]["kd_s_inv"]=value
    with pytest.raises(ValueError):config_from_dict(data)


def test_explicit_default_kd_preserves_short_hybrid_history(short_config):
    data=short_config.to_dict();data["controller"]["kd_s_inv"]=short_config.derived.controller_kd_s_inv
    explicit=config_from_dict(data)
    old=simulate_case(short_config,STANDARD_CASES[0]);new=simulate_case(explicit,STANDARD_CASES[0])
    assert [e.to_dict() for e in old.switch_events]==[e.to_dict() for e in new.switch_events]
    for key,value in old.arrays().items():np.testing.assert_array_equal(value,new.arrays()[key])


def test_velocity_distribution_and_phase_directions():
    e=np.array([-2,-2,2,2,0,1,1.]);v=np.array([-1,1,-1,1,2,0,1.])
    m=velocity_distribution(e,v)
    assert m["V_res_all_m_s"]==2
    assert m["V_res_rms_all_m_s"]==pytest.approx(np.sqrt(9/7))
    assert m["V_res_median_abs_all_m_s"]==1
    assert m["worsening_gap_fraction"]==pytest.approx(3/7)
    assert m["recovering_gap_fraction"]==pytest.approx(2/7)
    assert m["neutral_gap_fraction"]==pytest.approx(2/7)
    assert m["error_squared_mean_rate_m2_s"]==pytest.approx(2/7)
    assert m["worsening_weight_fraction"]==pytest.approx(5/9)
    assert sum(m["quadrant_counts"].values())==7


def test_zero_residual_and_invalid_vectors():
    m=velocity_distribution(np.zeros(3),np.zeros(3))
    assert m["neutral_gap_fraction"]==1 and m["worsening_weight_fraction"] is None
    for e,v in [([],[]),([1],[1,2]),([1],[np.nan])]:
        with pytest.raises(ValueError):velocity_distribution(e,v)


def test_endpoint_distribution_uses_referenced_signed_state_row():
    cfg=load_config("tests/fixtures/regression_24h/config.yaml")
    states=dict(gap_error_m=np.array([[0.,0.],[2.,-3.]]),gap_rate_m_s=np.array([[0.,0.],[-.001,-.002]]))
    report=dict(summary={},local_cycles=[dict(end_state_row=1,Q_D_res_all_m_s2=.003*cfg.derived.controller_kd_s_inv)],group_cycles=[],quiet_cycles=[])
    row=enhance_residual(report,states,cfg)["local_cycles"][0]
    assert row["V_res_all_m_s"]==.002
    assert row["V_res_rms_all_m_s"]==pytest.approx(np.sqrt(2.5e-6))
    assert row["worsening_gap_fraction"]==.5
    assert row["Q_D_res_div_Kd_m_s"]==pytest.approx(.003)


def test_first_day_window_excludes_late_errors_and_commands():
    cfg=load_config("tests/fixtures/regression_24h/config.yaml")
    arrays=dict(time_s=np.array([0.,86400.,86410.]),gap_error_m=np.array([[0.,0.],[3.,4.],[100.,100.]]))
    events=[dict(satellite_id=1,command_time=86410.,start_time=86410.,complete_time=86530.,from_mode="L",to_mode="H")]
    row=first_day_metrics(arrays,cfg,events)
    assert row["max_error_m"]==4 and row["rms_error_m"]==2.5
    assert row["transition_starts"]==0


@pytest.mark.parametrize("p,d,kind",[(3.,1.,"P_dominant"),(1.,3.,"D_dominant"),(1.,1.,"mixed"),(0.,0.,"zero")])
def test_contribution_categories(p,d,kind):
    assert contribution(p,d)["contribution_class"]==kind


def test_contribution_preserves_cancellation_signs():
    m=contribution(-1.,3.)
    assert m["D_magnitude_fraction"]==.75 and m["opposed_signs"]
    assert m["cancellation_ratio"]==.5


def test_same_physical_velocity_cancels_kd_in_P_growth():
    e=np.array([[1.,-2.,3.]]);v=np.array([[.001,-.002,.003]])
    kp=1/21600.**2;kd=2/21600.
    _,base=decompose_q(e,v,kp,kd)
    for factor in [.3,3.]:
        _,d=decompose_q(e,v,kp,kd*factor)
        np.testing.assert_allclose(d,base*factor,rtol=1e-15)
        np.testing.assert_allclose(kp/(kd*factor)*d,kp/kd*base,rtol=1e-15)


def test_all_LOW_match_allows_only_kd():
    cfg=load_config("tests/fixtures/regression_24h/config.yaml")
    data=cfg.to_dict();data["controller"]["kd_s_inv"]=cfg.derived.controller_kd_s_inv*3
    altered=config_from_dict(data)
    assert match_all_low(altered,cfg)
    data["controller"]["dwell_time_s"]=150
    assert not match_all_low(config_from_dict(data),cfg)


def test_crossing_and_command_separated_by_LOW_gate(short_config):
    cfg=short_config
    t=np.arange(0,71,10.)
    e=np.zeros((len(t),4));v=np.zeros_like(e)
    v[:,0]=np.array([0,0,2,2,1,2,2,2.])*cfg.controller.threshold_on_m_s2/cfg.derived.controller_kd_s_inv
    p,d=decompose_q(e,v,cfg.derived.controller_kp_s2_inv,cfg.derived.controller_kd_s_inv)
    a=dict(time_s=t,gap_error_m=e,gap_rate_m_s=v,q_m_s2=p+d)
    events=[dict(satellite_id=1,command_time=60.,start_time=60.,complete_time=70.,to_mode="H")]
    row=threshold_events(a,cfg,events)[0]
    assert row["sampled_crossing_time_s"]==20 and row["command_minus_crossing_s"]==40
    assert row["crossing_contribution_class"]=="D_dominant"


def test_frozen_baseline_return_gates_are_dwell_limited():
    from sat_string.io import load_case_arrays,load_switch_events
    root=Path("tests/fixtures/regression_24h")
    cfg=load_config(root/"config.yaml");a=load_case_arrays(root/"controlled_disturbed")
    # Fixture omits q; reconstruct the saved-geometry signal.
    p,d=decompose_q(a["gap_error_m"],a["gap_rate_m_s"],cfg.derived.controller_kp_s2_inv,cfg.derived.controller_kd_s_inv)
    a["q_m_s2"]=p+d
    rows=return_gates(a,cfg,load_switch_events(root/"controlled_disturbed"))
    assert len(rows)==20 and all(x["observed_hold_s"]==300 for x in rows)
    assert not any(x["off_not_met_at_gate"] for x in rows)


def test_unfinished_HIGH_return_is_not_lost_from_gate_statistics():
    cfg=load_config("tests/fixtures/regression_24h/config.yaml")
    t=np.arange(0,601,10.)
    q=np.zeros((len(t),21));q[:,0]=2*cfg.controller.threshold_off_m_s2
    events=[dict(satellite_id=1,command_time=10.,start_time=10.,complete_time=130.,to_mode="H")]
    row=return_gates(dict(time_s=t,q_m_s2=q),cfg,events)[0]
    assert row["earliest_return_gate_s"]==430 and row["off_not_met_at_gate"]
    assert row["return_right_censored"] and row["low_command_s"] is None
    assert not row["completed_by_tmax"]


def test_unobserved_return_gate_is_censored_rather_than_assumed_blocked():
    cfg=load_config("tests/fixtures/regression_24h/config.yaml")
    t=np.arange(0,201,10.)
    events=[dict(satellite_id=1,command_time=10.,start_time=10.,complete_time=130.,to_mode="H")]
    row=return_gates(dict(time_s=t,q_m_s2=np.zeros((len(t),21))),cfg,events)[0]
    assert not row["gate_observed"] and row["off_not_met_at_gate"] is None


def test_fixed_state_forecast_slope_is_invariant_but_margin_changes():
    cfg=load_config("tests/fixtures/regression_24h/config.yaml")
    e=np.zeros(20);v=np.zeros(20);v[0]=.0001
    predictions=fixed_state_forecasts(e,v,0.,np.zeros(21),cfg)
    assert len(predictions)==3
    assert all(p["satellite_id"]==1 for p in predictions)
    np.testing.assert_allclose([p["q_P_rate_m_s3"] for p in predictions],cfg.derived.controller_kp_s2_inv*.0001)
    assert predictions[0]["predicted_command_s"]>predictions[1]["predicted_command_s"]>predictions[2]["predicted_command_s"]


@pytest.mark.parametrize("actual",[40.,50.,None])
def test_indexed_event_windows_match_legacy_scan_at_boundaries(actual):
    from sat_string.residual import _compare,predict_crossing
    cfg=load_config("tests/fixtures/regression_24h/config.yaml")
    t=np.arange(0,61,10.);v=np.full((len(t),20),.0001);e=v*t[:,None]
    p,d=decompose_q(e,v,cfg.derived.controller_kp_s2_inv,cfg.derived.controller_kd_s_inv)
    arrays=dict(time_s=t,gap_error_m=e,gap_rate_m_s=v,q_m_s2=p+d)
    events=[dict(command_time=0.,complete_time=20.),dict(command_time=20.,complete_time=40.),dict(command_time=40.,complete_time=60.)]
    summary=dict(state_sample_time_s=30.,source_sample_index=3,all_low_at_sample=False)
    state=dict(q_D_m_s2=d[3]);kp,kd=cfg.derived.controller_kp_s2_inv,cfg.derived.controller_kd_s_inv
    pred=predict_crossing(p[3,0],d[3,0],30.,30.,kp,kd,cfg.controller.threshold_on_m_s2,10.)
    legacy=_compare(arrays,events,summary,state,0,pred,actual,cfg,.1)
    index=(np.sort([x["command_time"] for x in events]),np.sort([x["complete_time"] for x in events]))
    indexed=_compare(arrays,events,summary,state,0,pred,actual,cfg,.1,index)
    assert indexed==legacy
