import importlib.util
from pathlib import Path
import numpy as np
import pytest
from sat_string.observed_density import (gps_2003_to_utc,normalize_to_400,
    separation_differences,fit_gaussian_event,gaussian,geocentric_radius,read_igg)


def test_gps_utc_and_guard():
    assert gps_2003_to_utc(['2003-11-20T00:00:13'])[0]==np.datetime64('2003-11-20')
    with pytest.raises(ValueError):gps_2003_to_utc(['2010-01-01'])


def test_normalization_preserves_relative_anomaly():
    delta,rho=normalize_to_400([2.,6.],[1.,4.],[5.,8.])
    np.testing.assert_allclose(delta,[1.,.5]);np.testing.assert_allclose(rho,[10.,12.])


@pytest.mark.parametrize('bad',[0.,-1.,np.nan,np.inf])
def test_normalization_rejects_invalid(bad):
    with pytest.raises(ValueError):normalize_to_400(bad,1.,1.)


def test_d100_interpolation_and_no_flag_gap_bridging():
    s=np.arange(10)*75.;d=s*.001;t=np.datetime64('2003-10-25T01')+np.arange(10)*np.timedelta64(10,'s')
    valid=np.ones(10,bool);valid[4]=False
    diff,lag=separation_differences(s,d,t,valid,100)
    assert diff[0]==pytest.approx(.1);assert lag[0]==pytest.approx(100/7.5)
    assert np.isnan(diff[2]) and np.isnan(diff[4])


def test_driver_jump_not_classified_as_spatial_gradient():
    s=np.arange(6)*75.;t=np.datetime64('2003-10-25T02:59:30')+np.arange(6)*np.timedelta64(10,'s')
    d=np.array([0,0,0,10,10,10.]);diff,_=separation_differences(s,d,t,np.ones(6,bool),100)
    assert np.nanmax(np.abs(diff))==0
    assert np.isnan(diff[1]) and np.isnan(diff[2])


def test_gaussian_recovery_with_sloping_regional_background():
    pytest.importorskip('scipy')
    s=np.arange(-1500,1501,75.);delta=gaussian(s,.7,60,230,.2,.00005)
    fit=fit_gaussian_event(s,delta,0)
    assert fit['accepted'];assert fit['amplitude_delta']==pytest.approx(.7,rel=1e-4)
    assert fit['fwhm_km']==pytest.approx(2*np.sqrt(2*np.log(2))*230,rel=1e-4)
    assert fit['amplitude_fraction_local_bg']==pytest.approx(.7/1.2,rel=1e-4)


def test_linear_trend_not_accepted_as_a_pulse():
    pytest.importorskip('scipy')
    s=np.arange(-1500,1501,75.)
    assert not fit_gaussian_event(s,.4+s*.0001,0)['accepted']


def test_unresolved_gaussian_not_accepted():
    pytest.importorskip('scipy')
    s=np.arange(-1500,1501,75.)
    assert not fit_gaussian_event(s,gaussian(s,.5,0,42.5,0,0),0)['accepted']


def test_geocentric_radius_poles_and_equator():
    radius=geocentric_radius([0,90],[400000,400000])
    assert radius[0]==pytest.approx(6778137)
    assert radius[1]==pytest.approx(6378137*(1-1/298.257222101)+400000)


def test_real_geographic_msis_is_2_0_repeatable_and_si(monkeypatch):
    msis=pytest.importorskip('pymsis')
    from sat_string.observed_density import evaluate_nrlmsis_geographic
    class FrozenWeather:
        def at(self,t):return np.tile([150,150,4,4,4,4,4,4,4],(len(t),1))
    original=msis.calculate;calls=[]
    def recording(*a,**kw):
        calls.append(kw);return original(*a,**kw)
    monkeypatch.setattr(msis,'calculate',recording)
    dates=np.array(['2003-11-20T00:00:00','2003-11-20T00:00:00'],dtype='datetime64[s]')
    a=evaluate_nrlmsis_geographic(dates,[10,10],[50,50],[400,400],FrozenWeather())
    b=evaluate_nrlmsis_geographic(dates,[10,10],[50,50],[400,400],FrozenWeather())
    np.testing.assert_array_equal(a,b);assert a[0]==a[1]
    assert np.all((a>1e-14)&(a<1e-10)) # total mass density in kg/m3
    assert all(c['version']=='2.0' and c['geomagnetic_activity']==-1 for c in calls)


def test_actual_netcdf_flags_and_gps_units(tmp_path):
    nc=pytest.importorskip('netCDF4');path=tmp_path/'tiny.nc'
    with nc.Dataset(path,'w') as ds:
        ds.createDimension('time',20)
        def var(name,values,unit):
            v=ds.createVariable(name,'f8',('time',));v[:]=values;v.units=unit;return v
        v=var('time',np.arange(20)/6,'minutes since 2003-10-25 00:00:00');v.long_name='GPS time'
        for name,value,unit in [('lat',0,'degrees_north'),('lon',10,'degrees_east'),('alt',400000,'m'),
            ('density',4e-12,'kg/m^3'),('arglat',np.arange(20)*.6,'degree'),('lst',12,'h'),('shadow',0,'1')]:var(name,value,unit)
        gap=ds.createGroup('gap');g=gap.createVariable('orbit','i4',('time',));g[:]=0;g[5]=1
        flags=ds.createGroup('flag');v=flags.createVariable('path','i4',('time',));v[:]=1
        v=flags.createVariable('thrust','i4',('time',));v[:]=0;v[8]=1
    data,qc=read_igg(path,'2003-10-24T23:59:00','2003-10-25T01:00:00')
    assert data['utc'][0]==np.datetime64('2003-10-24T23:59:47')
    assert qc['accepted']==18;assert not data['valid'][5] and not data['valid'][8]
    assert data['valid'][0] # orbit direction 1 is not a rejection flag


def test_grace_interpolation_no_rejected_bin_or_gap():
    path=Path(__file__).resolve().parents[1]/'scripts/compare_grace_density.py'
    import sys
    sys.path.insert(0,str(path.parent))
    spec=importlib.util.spec_from_file_location('grace_comparison',path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    t=np.array([3600.,3650.,3700.,4000.]);valid=np.array([True,True,False,True])
    result=module.interpolate_good(t,np.arange(4.),valid,np.array([3625.,3675.,3800.,4100.]))
    assert result[0]==.5;assert np.all(np.isnan(result[1:]))
