from copy import deepcopy
import numpy as np
import pytest

from sat_string.config import config_from_dict
from sat_string.disturbance import earth_track_increment, fwhm_to_sigma
from sat_string.simulator import simulate_case, CaseDefinition
from sat_string.relative_motion import simulate_all_low_exposure, evaluate_motion


def test_center_independent_of_train_and_no_time_window():
    kwargs=dict(reference_radius_m=10,center_at_epoch_m=20,propagation_speed_m_s=3,
                sigma_m=4,amplitude_fraction=.49)
    for t in (-100,0,10000):
        longitude=(20+3*t)/10
        assert earth_track_increment([longitude],t,**kwargs)[0]==pytest.approx(.49)
        assert earth_track_increment([longitude+.4],t,**kwargs)[0]==pytest.approx(.49*np.exp(-.5))
        assert earth_track_increment([longitude],t,enabled=False,**kwargs)[0]==0


def test_fwhm_is_half_peak():
    sigma=fwhm_to_sigma(404000)
    val=earth_track_increment([202000],0,reference_radius_m=1,center_at_epoch_m=0,
        propagation_speed_m_s=0,sigma_m=sigma,amplitude_fraction=.49)
    assert val[0]==pytest.approx(.245)


@pytest.mark.parametrize("new_model",[False,True])
def test_passive_quadratures_preserve_physical_trajectory(short_config,new_model):
    d=deepcopy(short_config.to_dict())
    if new_model: d["disturbance"].update(type="earth_track_gaussian",center_mode="manual",
                                        center_s_m=300000,propagation_speed_m_s=730)
    config=config_from_dict(d)
    standard=simulate_case(config,CaseDefinition("all_low",False,True))
    logged,exposure=simulate_all_low_exposure(config)
    for key in standard.arrays():
        np.testing.assert_array_equal(standard.arrays()[key],logged.arrays()[key])
    assert np.all(np.diff(exposure["J_rho"],axis=0)>=0)
    assert np.all(exposure["J_drag"][-1]<0)
    metrics=evaluate_motion(logged.arrays(),exposure,config)
    assert metrics["exposure"]["total_exposure_radius_identity_max_error_m"]<1e-7


def test_manual_epoch_required_and_duration_not_gate(short_config):
    d=short_config.to_dict()
    d["disturbance"]["type"]="earth_track_gaussian"
    with pytest.raises(ValueError,match="manual"): config_from_dict(d)
    d["disturbance"]["center_mode"]="manual"
    d["disturbance"]["duration_s"]=1e9
    assert config_from_dict(d).disturbance.duration_s==1e9


def test_frozen_gaussian_passages_have_equal_integrals():
    # Analytic full Gaussian area; no simulation assumptions about decay.
    t=np.arange(-1000,1000,.25)
    offset=np.array([-100000,0,100000])
    speed=7000
    rho=earth_track_increment(speed*t[:,None]+offset,t[:,None],
        reference_radius_m=1,center_at_epoch_m=0,propagation_speed_m_s=0,
        sigma_m=171000,amplitude_fraction=.49)
    integrals=np.trapezoid(rho,t,axis=0)
    np.testing.assert_allclose(integrals,.49*np.sqrt(2*np.pi)*171000/speed,rtol=1e-12)


def test_comoving_not_mislabeled_as_post_passage(short_config):
    d=short_config.to_dict()
    d["disturbance"].update(type="earth_track_gaussian",center_mode="manual",center_s_m=0,
        propagation_speed_m_s=short_config.derived.reference_radius_m*short_config.derived.reference_mean_motion_rad_s)
    cfg=config_from_dict(d)
    result,exposure=simulate_all_low_exposure(cfg)
    m=evaluate_motion(result.arrays(),exposure,cfg)
    assert m["motion"]["fwhm_encounter_time_s"] is None
    assert m["post_passage"]=={}
    assert m["frozen_rate_forecast_30d"]["max_abs_error_m"] is None


def test_disabled_exposure_has_no_nan_ratio(short_config):
    d=short_config.to_dict()
    d["disturbance"].update(type="earth_track_gaussian",center_mode="manual",
                            center_s_m=0,enabled=False)
    cfg=config_from_dict(d)
    result,exposure=simulate_all_low_exposure(cfg)
    metrics=evaluate_motion(result.arrays(),exposure,cfg)
    assert metrics["exposure"]["J_delta_spread_kg_s_m3"]==0
    assert metrics["exposure"]["relative_J_delta_spread"] is None
