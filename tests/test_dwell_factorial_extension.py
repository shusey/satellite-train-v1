from copy import deepcopy
from pathlib import Path
import numpy as np
import pytest
from sat_string.config import load_config
from sat_string.dwell_factorial import factorial_config,extension_config,period_cycle_diagnostics


def fixture(high=150,low=300):
    base=load_config(Path(__file__).parent/'fixtures/regression_24h/config.yaml')
    cfg=factorial_config(base,high,low)
    profile=dict(delta_a_m=[0.]*21,B_multiplier=[1.]*21)
    state=[np.full(21,cfg.derived.reference_radius_m).tolist(),(np.arange(21)*cfg.derived.target_spacing_rad).tolist()]
    return dict(configuration=cfg.to_dict(),profiles=[profile],disturbance_enabled=False,time_s=30*86400.,state=state),profile


@pytest.mark.parametrize('high,low',[(150,300),(300,150)])
def test_extension_changes_only_duration_and_preserves_source(high,low):
    cp,profile=fixture(high,low);before=deepcopy(cp)
    result=extension_config(cp,profile)
    expected=deepcopy(cp['configuration']);expected['simulation']['duration_s']=60*86400.
    assert result.to_dict()==expected and cp==before


@pytest.mark.parametrize('fault',['time','profile','disturbance','nonfinite','reversal','radius','shape','diagonal'])
def test_extension_rejects_invalid_endpoint(fault):
    cp,p=fixture()
    if fault=='time':cp['time_s']=7*86400.
    elif fault=='profile':cp['profiles']=[]
    elif fault=='disturbance':cp['disturbance_enabled']=True
    elif fault=='nonfinite':cp['state'][0][0]=float('nan')
    elif fault=='reversal':cp['state'][1][1]=cp['state'][1][0]-1
    elif fault=='radius':cp['state'][0][0]=0.
    elif fault=='shape':cp['state']=[[1.],[2.]]
    else:cp['configuration']['controller']['dwell_low_time_s']=150.
    with pytest.raises(ValueError):extension_config(cp,p)


def test_late_cycle_window_preserves_completion_boundary_and_censoring():
    rows=[dict(satellite_id=1,cycle=k,end_s=t,next_HIGH_s=t+10 if t<60 else None,
        next_wait_s=10 if t<60 else None,pulse_g_impulse_m_s=k/1000,V_incident_after_m_s=k/2000,
        next_qD_neighbor=2e-7,next_qD_natural=1e-10) for k,t in enumerate([29,30,60,61],1)]
    edges=[dict(satellite_id=1,cycle=k,velocity_after_m_s=k/2000,worsening_after=True) for k in range(1,5)]
    result=period_cycle_diagnostics(rows,edges,30,60)
    assert result['completed_cycles']==2
    assert result['re_trigger']['right_censored']==1
    assert result['re_trigger']['observed_wait_s']['n']==1
    assert result['residual_abs_incident_gap_velocity_m_s']['median']==pytest.approx(.00125)
    assert result['neighbor_exceeds_natural_fraction']==1


def test_empty_late_window_does_not_invent_residual_or_retrigger():
    result=period_cycle_diagnostics([],[],30,60)
    assert result['completed_cycles']==0
    assert result['residual_abs_incident_gap_velocity_m_s']['median'] is None
    assert result['neighbor_exceeds_natural_fraction'] is None
