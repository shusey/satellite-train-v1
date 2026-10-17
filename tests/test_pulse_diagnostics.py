import numpy as np
import pytest
from sat_string.pulse_diagnostics import (acceleration_terms, cumulative_trapezoid,
                                         local_control_split, node_map, reversal_flags)


def test_chain_rule_and_three_way_split():
    a=np.array([6.7e6,6.8e6]);rho=np.array([3e-12,5e-12]);ar=np.array([.01,.03]);eps=np.array([.005,-.005])
    mu=3.986e14;r=6.778e6;cd=2.2;m=4.
    terms=acceleration_terms(a,rho,ar,eps,r_ref=r,mu=mu,cd=cd,mass=m,area_low=.01)
    adot=-rho*cd*ar/m*(1+eps)*np.sqrt(mu*a)
    expected=-1.5*r*np.sqrt(mu)*a**(-2.5)*adot
    np.testing.assert_allclose(sum(terms.values()),expected,rtol=1e-14)
    assert terms['control'][0]==0 and terms['mismatch'][1]<0


def test_equal_adjacent_pulses_cancel():
    own,other=local_control_split(np.array([.002,.002,.002]),1)
    np.testing.assert_allclose(own+other,0,atol=1e-20)
    np.testing.assert_allclose(own,[.002,-.002])
    assert node_map(own)[1]==pytest.approx(-.004)


def test_boundary_satellite_and_incidence():
    own,other=local_control_split(np.array([.002,0.,0.]),0)
    np.testing.assert_allclose(own,[-.002,0]);assert not other.any()
    np.testing.assert_allclose(node_map(own),[-.002,.002,0])


def test_integral_affine_nonuniform_samples():
    t=np.array([0.,2.,5.,9.]);y=np.column_stack((2*t+3,-t))
    np.testing.assert_allclose(cumulative_trapezoid(y,t),np.column_stack((t*t+3*t,-t*t/2)))
    with pytest.raises(ValueError):cumulative_trapezoid(y,[0,2,2,9])


def test_reversal_distinguished_from_growth_and_neighbor_cancellation():
    flags=reversal_flags(np.array([1.,1.,0.]),np.array([-.5,2.,1.]),np.array([-2.,-2.,-2.]))
    assert flags['actual_sign_reversal'].tolist()==[True,False,False]
    assert flags['actual_abs_growth'].tolist()==[False,True,True]
    assert flags['own_alone_crosses_zero'].tolist()==[True,True,False]
