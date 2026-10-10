from dataclasses import replace
from pathlib import Path
import numpy as np
import pytest
from sat_string.config import load_config,config_from_dict
from sat_string.dynamics import density_at_state,orbital_rates
from sat_string.nrlmsis import OrbitGeometry,SpaceWeather,NRLMSIS2,TrajectoryDensityCache,WGS84_A,WGS84_F,ERA_RATE

WEATHER=Path(__file__).parent/'fixtures/nrlmsis2/drivers.csv'


def geometry(inc=51.6):return OrbitGeometry('2024-01-01T00:00:00',inc,0.,0.)


def backend():
    pytest.importorskip('pymsis')
    return NRLMSIS2(geometry(),SpaceWeather(WEATHER,geometry().epoch_utc))


def test_equatorial_units_and_rotation():
    g=geometry(0);lon,lat,alt=g.coordinates([0,10],np.full(2,WGS84_A+400000),np.zeros(2))
    np.testing.assert_allclose(alt,400,atol=1e-9)
    np.testing.assert_allclose(lat,0,atol=1e-12)
    np.testing.assert_allclose(np.deg2rad(lon[1]-lon[0]),-10*ERA_RATE,atol=1e-12)


def test_geodetic_inverse_matches_inclined_cartesian_position():
    g=geometry();phase=np.linspace(-3,3,100)
    a=np.full(100,WGS84_A+400000);lon,lat,alt=g.coordinates(0,a,phase)
    lat=np.deg2rad(lat);alt=alt*1000;e2=WGS84_F*(2-WGS84_F)
    n=WGS84_A/np.sqrt(1-e2*np.sin(lat)**2)
    p=(n+alt)*np.cos(lat);z=(n*(1-e2)+alt)*np.sin(lat)
    np.testing.assert_allclose(z,a*np.sin(np.deg2rad(51.6))*np.sin(phase),atol=1e-7)
    np.testing.assert_allclose(p,np.sqrt(a*a-z*z),atol=1e-7)
    assert np.ptp(alt)>10000  # geodetic altitude is not a - spherical Earth radius


def test_coordinate_wrap_is_physically_continuous():
    g=geometry();phase=np.linspace(-4*np.pi,4*np.pi,2001)
    lon,lat,alt=g.coordinates(0,np.full_like(phase,WGS84_A+400000),phase)
    assert np.max(np.abs(np.diff(np.unwrap(np.deg2rad(lon)))))<.04
    assert np.max(np.abs(np.diff(lat)))<1
    assert np.all(np.isfinite(alt))


def test_pole_altitude_and_invalid_geometry():
    _,lat,alt=geometry(90).coordinates(0,WGS84_A+400000,np.pi/2)
    assert lat==pytest.approx(90)
    assert alt==pytest.approx((WGS84_A+400000-WGS84_A*(1-WGS84_F))/1000)
    with pytest.raises(ValueError):geometry(181)


def test_weather_is_three_hour_piecewise_constant():
    w=SpaceWeather(WEATHER,geometry().epoch_utc)
    np.testing.assert_array_equal(w.at([0]),w.at([10799]))
    np.testing.assert_array_equal(w.at([10800]),w.values[1:2])
    with pytest.raises(ValueError):w.at([-1])
    with pytest.raises(ValueError):w.at([w.seconds[-1]+1])


def test_legacy_normalization_does_not_add_density_model():
    c=load_config('tests/fixtures/regression_24h/config.yaml')
    assert 'density_model' not in c.to_dict()['atmosphere']
    assert config_from_dict(c.to_dict()).to_dict()==c.to_dict()
    c=replace(c,atmosphere=replace(c.atmosphere,density_model='nrlmsis2'))
    with pytest.raises(ValueError):density_at_state(0,[6778137],[0],c,disturbance_enabled=False)
    with pytest.raises(ValueError):density_at_state(0,[6778137],[0],c,disturbance_enabled=True,density_provider=lambda *x:[1e-12])


def test_actual_nrlmsis2_is_reproducible_and_physical():
    b=backend();a=np.full(5,6778137.);u=np.linspace(-1,1,5)
    x=b(3600,a,u);y=b(3600,a,u)
    np.testing.assert_array_equal(x,y)
    assert np.all((x>1e-14)&(x<1e-9))
    assert np.ptp(x)>0


def test_explicit_version_and_mass_units(monkeypatch):
    b=backend();import pymsis
    seen={}
    def fake(dates,lons,lats,alts,f107s,f107as,aps,**kw):
        seen.update(kw);seen['alt_km']=alts;seen['ap']=aps
        result=np.zeros((len(dates),11));result[:,0]=1e-12;return result
    monkeypatch.setattr(pymsis,'calculate',fake)
    rho=b(0,np.array([6778137.]),np.array([0.]))
    assert seen['version']=='2.0' and seen['geomagnetic_activity']==-1
    assert seen['alt_km'][0]==pytest.approx(400)
    assert seen['ap'].shape==(1,7) and rho[0]==1e-12
    da,_=orbital_rates([6778137],rho,[.01],mu_m3_s2=398600441800000,drag_coefficient=2.2,mass_kg=4)
    assert da[0]==pytest.approx(-1e-12*2.2*.01/4*np.sqrt(398600441800000.*6778137.))


@pytest.mark.parametrize('step',[10,30,60,120])
def test_cached_density_matches_direct_along_actual_decay(step):
    b=backend();cache=TrajectoryDensityCache(b,step,398600441800000)
    a0=np.full(3,6778137.);u0=np.array([-.15,0,.15]);n=np.sqrt(398600441800000/a0**3)
    for time in range(0,241,2):
        a=a0-.001*time;u=u0+n*time
        np.testing.assert_allclose(cache(time,a,u),b(time,a,u),rtol=.001,atol=0)


def test_cache_midnight_and_weather_jump_do_not_smear_drivers():
    b=backend()
    for time in [10799.,10800.,86399.,86400.]:
        cache=TrajectoryDensityCache(b,120,398600441800000)
        a=np.full(3,6778137.);u=np.array([-.15,0,.15])
        np.testing.assert_allclose(cache(time,a,u),b(time,a,u),rtol=1e-6,atol=0)


def test_optional_density_hook_preserves_orbital_core_with_identical_force():
    from sat_string.simulator import simulate_case,CaseDefinition
    c=load_config('tests/fixtures/regression_24h/config.yaml')
    c=replace(c,simulation=replace(c.simulation,duration_s=120))
    external=replace(c,atmosphere=replace(c.atmosphere,density_model='nrlmsis2'))
    case=CaseDefinition('same_force',True,False)
    old=simulate_case(c,case)
    new=simulate_case(external,case,density_provider=lambda t,a,u:density_at_state(t,a,u,c,disturbance_enabled=False))
    for key,value in old.arrays().items():np.testing.assert_array_equal(value,new.arrays()[key])
    assert old.switch_events==new.switch_events
