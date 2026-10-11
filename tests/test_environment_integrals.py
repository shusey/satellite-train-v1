"""Synthetic analytic and missing-support tests; no expensive physical sweep."""
import importlib.util
from pathlib import Path
import numpy as np
import pytest
from sat_string.environment_integrals import all_low_integrals,gap_rate_from_exposure,supported_integral


def test_constant_density_orbit_identity():
    mu,cd,area,mass,a0=3.986004418e14,2.2,.01,4.,6778137.
    rho=np.array([3e-12,4e-12,5e-12]);t=86400.
    a=(np.sqrt(a0)-cd*area*np.sqrt(mu)/(2*mass)*rho*t)**2
    jr,jd=all_low_integrals(a,a0,mu,cd,area,mass)
    # Input radius rounding at6.8Mm limits reconstruction to about1e-17 SI.
    np.testing.assert_allclose(jr,rho*t,rtol=2e-12,atol=1e-17)
    assert np.all(jd<0)
    np.testing.assert_allclose(jd,np.sqrt(mu)*(1/np.sqrt(a0)-1/np.sqrt(a)),rtol=5e-12,atol=1e-11)
    predicted=gap_rate_from_exposure(jr,a0,a0,mu,cd,area,mass)
    np.testing.assert_allclose(predicted,a0*np.diff(np.sqrt(mu/a**3)),rtol=5e-9,atol=1e-12)
    assert np.all(predicted>0)  # Higher drag lowers a and increases mean motion.


def test_equal_common_density_no_gap_rate():
    np.testing.assert_array_equal(gap_rate_from_exposure([1e-6]*3,6778137,6778137,3.986e14,2.2,.01,4),[0,0])


def test_initial_integrals_are_zero():
    jr,jd=all_low_integrals(np.full((2,3),6778137.),6778137.,3.986e14,2.2,.01,4)
    assert not jr.any() and not jd.any()


def test_missing_support_does_not_bridge_or_claim_full_duration():
    r=supported_integral([0,10,20,30,100,110],[1,1,999,1,1,1],[True,True,False,True,True,True],15,
                         uncertainty=[.1,.1,999,.1,.1,.1])
    assert r['integral_observed_support']==20
    assert r['coverage_s']==20 and r['max_contiguous_coverage_s']==10
    assert len(r['segments'])==2
    assert r['fully_correlated_error_envelope']==2


def test_signed_integral_cancellation_and_absolute_area():
    r=supported_integral([0,1,2,3],[1,1,-1,-1],[True]*4,1)
    assert r['integral_observed_support']==0 and r['absolute_interval_area']==2


@pytest.mark.parametrize('times',[[0,0],[1,0]])
def test_nonincreasing_times_rejected(times):
    with pytest.raises(ValueError):supported_integral(times,[1,1],[True,True],20)


def test_invalid_physics_rejected():
    with pytest.raises(ValueError):all_low_integrals([0],1,1,1,1,1)


@pytest.mark.parametrize('kp,days,dst,expected',[(4.3,0,-27,'quiet'),(8.3,7,-128,'active'),(9,10,-406,'strong'),(5,0,-40,'unclassified'),(8,1,-150,'unclassified'),(8,2,-250,'strong')])
def test_weather_only_classification(kp,days,dst,expected):
    spec=importlib.util.spec_from_file_location('prepare_environment_comparison',Path(__file__).parents[1]/'scripts/prepare_environment_comparison.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    assert module.weather_class(kp,days,dst)==expected
