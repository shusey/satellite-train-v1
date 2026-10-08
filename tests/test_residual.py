from dataclasses import replace
from pathlib import Path
import json
import runpy
import numpy as np
import pytest
import yaml
from sat_string.config import load_config
from sat_string.controller import switching_function
from sat_string.io import load_case_arrays,load_switch_events
from sat_string.long_horizon import decompose_q
from sat_string.residual import end_state,predict_crossing,cycle_pairs,residual_analysis,quiet_packets
from sat_string.simulator import simulate_case,STANDARD_CASES


def synthetic(config):
    t=np.arange(0,61,10,dtype=float);n=config.formation.satellite_count
    v=np.tile(np.array([.001,-.002,.003,.004]),(len(t),1))
    e=t[:,None]*v
    q=np.array([switching_function(x,y,config.derived.controller_kp_s2_inv,
        config.derived.controller_kd_s_inv) for x,y in zip(e,v)])
    return dict(time_s=t,gap_error_m=e,gap_rate_m_s=v,q_m_s2=q,mode=np.zeros((len(t),n),dtype=np.int8))


def test_residual_velocity_and_D_scopes(short_config):
    a=synthetic(short_config)
    row,state=end_state(a,30,[2,4],[2],short_config)
    assert row['V_res_all_m_s']==.004
    assert row['V_res_seeds_m_s']==.002
    assert row['V_res_participants_m_s']==.004
    kd=short_config.derived.controller_kd_s_inv
    assert row['Q_D_res_all_m_s2']==pytest.approx(.005*kd)
    assert row['Q_D_res_seeds_m_s2']==pytest.approx(.003*kd)
    assert state['q_D_m_s2'][1]<0
    np.testing.assert_array_equal(state['gap_rate_m_s'],a['gap_rate_m_s'][3])


def test_first_post_completion_sample_is_explicit(short_config):
    row,_=end_state(synthetic(short_config),12,[1],[1],short_config)
    assert row['state_sample_time_s']==20 and row['sample_offset_s']==8


def test_positive_D_crossing_and_control_grid():
    p=predict_crossing(1e-7,1e-8,100,400,2e-9,1e-4,1.5e-7,10)
    assert p['predicted_crossing_s']==pytest.approx(200100)
    assert p['predicted_command_s']==200100
    assert p['q_P_rate_m_s3']==pytest.approx(2e-13)


@pytest.mark.parametrize('D',[0.,-1e-8])
def test_nonpositive_D_has_no_upward_forecast(D):
    assert predict_crossing(1e-7,D,0,0,2e-9,1e-4,1.5e-7,10)['predicted_command_s'] is None


def test_LOW_dwell_gate_and_downward_missed_crossing():
    rising=predict_crossing(1.49e-7,2e-9,100,407,2e-9,1e-4,1.5e-7,10)
    assert rising['predicted_crossing_s']==100 and rising['predicted_command_s']==410
    falling=predict_crossing(2e-7,-1e-8,0,100,1,1,1.5e-7,10)
    assert falling['status']=='falls_below_on_before_allowed' and falling['predicted_command_s'] is None


def test_local_cycle_pairing_preserves_unfinished_cycle():
    events=[dict(satellite_id=1,command_time=t,start_time=t,complete_time=t+10,
                 from_mode='L' if mode=='H' else 'H',to_mode=mode) for t,mode in [(0,'H'),(20,'L'),(40,'H')]]
    pairs=cycle_pairs(events,50)
    assert len(pairs)==2 and pairs[0]['return_event']['complete_time']==30
    assert pairs[1]['return_event'] is None and pairs[1]['local_cycle_number']==2


def test_frozen_24h_endpoint_extraction_and_censoring():
    root=Path(__file__).parent/'fixtures/regression_24h'
    cfg=load_config(root/'config.yaml');a=load_case_arrays(root/'controlled_disturbed')
    p,d=decompose_q(a['gap_error_m'],a['gap_rate_m_s'],cfg.derived.controller_kp_s2_inv,cfg.derived.controller_kd_s_inv)
    a['q_m_s2']=p+d
    report,states=residual_analysis(a,cfg,load_switch_events(root/'controlled_disturbed'),gap_s=80953.6935537842)
    assert report['summary']['local_completed_cycles']==20
    assert report['summary']['dwell_limited_returns']==20
    assert len(report['group_cycles'])==1 and report['group_cycles'][0]['all_low_at_sample']
    group=report['group_cycles'][0]
    assert group['completion_time_s']==62000 and group['right_censored']
    assert group['prediction']['predicted_command_s']>86400
    assert len(states['gap_error_m'])==20+1+len(report['quiet_cycles'])
    assert all(row['all_low_at_sample'] for row in report['quiet_cycles'])
    assert all(row['sample_offset_s']==0 for row in report['local_cycles'])
    assert any('intervening_commands' in row['invalid_reasons'] for row in report['local_cycles'])
    json.dumps(report,allow_nan=False)


@pytest.mark.parametrize('dwell',[150,300,600])
def test_dwell_configuration_readback_changes_only_requested_fields(tmp_path,dwell):
    cfg=load_config(Path(__file__).parents[1]/'tests/fixtures/regression_24h/config.yaml')
    builder=runpy.run_path(str(Path(__file__).parents[1]/'scripts/run_dwell_case.py'))['dwell_config']
    result=builder(cfg,dwell,30)
    path=tmp_path/'config.yaml';path.write_text(yaml.safe_dump(result.to_dict()),encoding='utf-8')
    assert load_config(path)==result
    expected=cfg.to_dict();expected['controller']['dwell_time_s']=dwell;expected['simulation']['duration_s']=30*86400
    assert result.to_dict()==expected
    assert result.derived.controller_kp_s2_inv==cfg.derived.controller_kp_s2_inv
    assert result.derived.controller_kd_s_inv==cfg.derived.controller_kd_s_inv


def test_pre_switch_LOW_prefix_is_dwell_independent(short_config):
    original=replace(short_config,simulation=replace(short_config.simulation,duration_s=14.))
    baseline=simulate_case(original,STANDARD_CASES[0]);assert not baseline.switch_events
    for dwell in [150.,600.]:
        cfg=replace(original,controller=replace(original.controller,dwell_time_s=dwell))
        alternative=simulate_case(cfg,STANDARD_CASES[0])
        assert not alternative.switch_events
        for key,value in baseline.arrays().items():np.testing.assert_array_equal(value,alternative.arrays()[key])


def test_quiet_packets_preserve_short_gaps_and_merge_overlaps():
    pairs=[]
    for sat,start,stop in [(1,10,100),(2,50,150),(1,160,200),(2,200,250)]:
        pairs.append(dict(satellite_id=sat,rise=dict(command_time=start),return_event=dict(complete_time=stop)))
    packets=quiet_packets(pairs,300)
    assert len(packets)==2
    assert packets[0]['end_time_s']==150 and packets[0]['transition_starts']==4
    assert packets[1]['start_time_s']==160 and packets[1]['end_time_s']==250
    unfinished=pairs+[dict(satellite_id=1,rise=dict(command_time=280),return_event=None)]
    assert not quiet_packets(unfinished,300)[-1]['complete']


def test_all_LOW_reference_is_dwell_independent(short_config):
    baseline=simulate_case(short_config,STANDARD_CASES[2])
    for dwell in [150.,600.]:
        cfg=replace(short_config,controller=replace(short_config.controller,dwell_time_s=dwell))
        result=simulate_case(cfg,STANDARD_CASES[2])
        assert not result.switch_events
        for key,value in baseline.arrays().items():np.testing.assert_array_equal(value,result.arrays()[key])
