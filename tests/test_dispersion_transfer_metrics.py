import numpy as np
import pytest
from sat_string.config import load_config
from sat_string.controller import initial_formation
from sat_string.dispersion_control import derived_histories
from sat_string.dispersion_transfer_metrics import evaluate_transfer, quiet_candidate
from sat_string.evaluation import EvaluationSettings

C=load_config('tests/fixtures/regression_24h/config.yaml')


def data():
    a,u=initial_formation(21,C.derived.reference_radius_m,C.derived.target_spacing_rad)
    t=np.arange(0.,21.,10.);aa=np.tile(a,(3,1));uu=np.tile(u,(3,1))
    return derived_histories(t,aa,uu,np.full((3,21),.01),np.zeros((3,21),dtype='int8'),np.full((3,21),4e-12),C)


def test_uncontrolled_zero_cost_and_signed_gap():
    a=data();r=evaluate_transfer(a,C,[],a)
    assert r['transitions']==r['completed_cycles']==0
    assert r['extra_mean_altitude_loss_m']==r['high_equivalent_total_s']==0
    assert r['min_signed_reference_gap_m']==pytest.approx(100000)
    assert r['first_label_order_reversal_s'] is None


def test_extra_altitude_uses_same_profile_OFF_not_initial_altitude():
    off=data();on=data();on['a_m'][-1]-=2
    r=evaluate_transfer(on,C,[],off)
    assert r['extra_mean_altitude_loss_m']==2
    assert r['extra_loss_per_satellite_m']==[2]*21


def test_mismatched_initial_state_rejected():
    off=data();on=data();on['a_m'][0,0]+=.1
    with pytest.raises(ValueError,match='initial a'):evaluate_transfer(on,C,[],off)


def test_reversal_preserves_labels_and_flags_invalid_precision():
    a=data();a['lambda_rad'][-1,1:]-=100100/C.derived.reference_radius_m
    a=derived_histories(*[a[k] for k in ('time_s','a_m','lambda_rad','area_m2','mode','density_kg_m3')],C)
    r=evaluate_transfer(a,C,[],a)
    assert r['first_label_order_reversal_s']==20
    assert r['allLOW_order_reversal_s']==20
    assert r['min_signed_reference_gap_m']==pytest.approx(-100,abs=1e-8)
    assert r['min_physical_pair_chord_m']>=0
    assert not r['full_horizon_precision_reduction_valid']
    assert r['matched_pre_allLOW_reversal_improvement']['last_valid_sample_s']==10


def test_accepted_error_requires_ALL_LOW_and_full_verification_window():
    t=np.arange(0.,51.,10.);e=np.zeros((6,20));v=e.copy();mode=np.zeros((6,21),dtype='int8')
    mode[2,0]=2
    r=quiet_candidate(t,e,v,mode,EvaluationSettings(50,.001,20))
    assert r['candidate_start_s']==30
    assert r['verification_horizon_s']==20
    assert quiet_candidate(t,e,v,mode,EvaluationSettings(50,.001,30))['candidate_start_s'] is None


def test_exceedance_RMS_and_percent_keep_physical_units():
    off=data();on=data();off['lambda_rad'][1:,1:]+=600/C.derived.reference_radius_m
    off=derived_histories(*[off[k] for k in ('time_s','a_m','lambda_rad','area_m2','mode','density_kg_m3')],C)
    r=evaluate_transfer(off,C,[],off)
    assert r['max_error_m']==pytest.approx(600,abs=1e-8)
    assert r['max_error_percent_of_100km']==pytest.approx(.6,abs=1e-9)
    assert r['rms_error_m']==pytest.approx(np.sqrt(2*600**2/60),abs=1e-8)
    assert r['exceedances'][0]['exceedance_fraction']==.5
