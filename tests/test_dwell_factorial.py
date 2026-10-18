from pathlib import Path
import numpy as np
import pytest
from sat_string.config import load_config, config_from_dict
from sat_string.attitude import AttitudeStateMachine, H, L
from sat_string.continuation import restore_attitude
from sat_string.dwell_factorial import factorial_config, contrasts


def base():
    return load_config(Path(__file__).parent/'fixtures/regression_24h/config.yaml')


@pytest.mark.parametrize('high,low,restart', [(300,300,840),(150,300,690),(300,150,690),(150,150,540)])
def test_distinct_arrival_gates_and_replay(high, low, restart):
    cfg = factorial_config(base(), high, low)
    m = AttitudeStateMachine(1,.01,.03,120,300,
        dwell_high_time_s=high, dwell_low_time_s=low)
    m.apply_commands(0,[H]); m.process_completions(120)
    assert not m.command_allowed(120+high-1)[0]
    assert m.command_allowed(120+high)[0]
    restored = restore_attitude(cfg, m.events, 120)
    assert restored.next_allowed_time[0] == 120+high
    m.apply_commands(120+high,[L]); m.process_completions(240+high)
    assert not m.command_allowed(restart-1)[0]
    assert m.command_allowed(restart)[0]
    restored = restore_attitude(cfg, m.events, 240+high)
    assert restored.next_allowed_time[0] == restart
    m.apply_commands(restart,[H])
    assert m.events[-1].start_time == restart


def test_legacy_serialization_and_single_override_fallback():
    cfg=base(); raw=cfg.to_dict()
    assert 'dwell_high_time_s' not in raw['controller']
    assert config_from_dict(raw) == cfg
    raw['controller']['dwell_high_time_s']=150
    new=config_from_dict(raw)
    assert new.controller.effective_dwell_high_s==150
    assert new.controller.effective_dwell_low_s==cfg.controller.dwell_time_s


@pytest.mark.parametrize('field', ['dwell_high_time_s','dwell_low_time_s'])
@pytest.mark.parametrize('value', [-1,0,float('nan'),float('inf')])
def test_invalid_override_rejected(field,value):
    raw=base().to_dict();raw['controller'][field]=value
    with pytest.raises(ValueError):config_from_dict(raw)


@pytest.mark.parametrize('dwell', [150,300])
def test_explicit_diagonal_matches_legacy_event_state_and_area(dwell):
    old=AttitudeStateMachine(1,.01,.03,120,dwell)
    new=AttitudeStateMachine(1,.01,.03,120,300,dwell_high_time_s=dwell,dwell_low_time_s=dwell)
    for t in range(0,3000,10):
        for m in [old,new]:
            m.process_completions(t)
            if m.command_allowed(t)[0]: m.apply_commands(t,[H if m.mode[0]==L else L])
        np.testing.assert_array_equal(old.mode,new.mode)
        np.testing.assert_array_equal(old.areas(t),new.areas(t))
        np.testing.assert_array_equal(old.next_allowed_time,new.next_allowed_time)
        assert old.events==new.events


def test_descriptive_contrast_definitions():
    assert contrasts([10,7,12,9])==dict(HIGH=-3,LOW=2,interaction=0)
    assert contrasts([10,10,10,4])==dict(HIGH=-3,LOW=-3,interaction=-6)


@pytest.mark.parametrize('high,low', [(150,300),(300,150)])
def test_residual_analysis_uses_high_return_gate_and_low_forecast_gate(high,low):
    from sat_string.residual import residual_analysis
    from sat_string.long_horizon import decompose_q
    cfg=factorial_config(base(),high,low)
    t=np.arange(0.,2001.,10.)
    e=np.zeros((len(t),20));v=e.copy();v[:,0]=.002
    p,d=decompose_q(e,v,cfg.derived.controller_kp_s2_inv,cfg.derived.controller_kd_s_inv)
    mode=np.zeros((len(t),21),dtype=np.int8)
    mode[t<120,0]=1;mode[(t>=120)&(t<120+high),0]=2
    mode[(t>=120+high)&(t<240+high),0]=3
    events=[dict(satellite_id=1,command_time=0.,start_time=0.,complete_time=120.,from_mode='L',to_mode='H'),
            dict(satellite_id=1,command_time=120.+high,start_time=120.+high,complete_time=240.+high,from_mode='H',to_mode='L')]
    arrays=dict(time_s=t,gap_error_m=e,gap_rate_m_s=v,q_m_s2=p+d,mode=mode)
    report,_=residual_analysis(arrays,cfg,events,groups=[])
    cycle=report['local_cycles'][0]
    assert cycle['earliest_low_command_s']==120+high
    assert cycle['high_equivalent_time_s']==120+high
    assert cycle['return_driver']=='dwell_limited'
    assert cycle['prediction']['predicted_command_s']==240+high+low
