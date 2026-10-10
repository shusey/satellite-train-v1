"""NRLMSIS 2.0 only: geographic geometry, explicit historical indices, local cache.

No gravity/control/drag changes. Input a is geocentric radius, lambda is orbital
phase (NOT geographic longitude). WGS84 supplies geodetic altitude/latitude.
"""
from dataclasses import dataclass
from pathlib import Path
import csv
import os
import numpy as np

WGS84_A = 6378137.0
WGS84_F = 1 / 298.257223563
ERA_RATE = 2*np.pi*1.00273781191135448/86400


@dataclass(frozen=True)
class OrbitGeometry:
    epoch_utc: str
    inclination_deg: float
    raan_deg: float
    initial_phase_deg: float

    def __post_init__(self):
        if not 0 <= self.inclination_deg <= 180 or not np.all(np.isfinite([self.inclination_deg,self.raan_deg,self.initial_phase_deg])):
            raise ValueError('Invalid orbit geometry')
        if np.isnat(np.datetime64(self.epoch_utc)):
            raise ValueError('Invalid epoch')

    @property
    def era0(self):
        days=(np.datetime64(self.epoch_utc,'us')-np.datetime64('2000-01-01T12:00:00','us'))/np.timedelta64(1,'D')
        return 2*np.pi*((.7790572732640+1.00273781191135448*float(days))%1)

    def coordinates(self,time_s,radius_m,phase_rad):
        t,a,phase=np.broadcast_arrays(np.asarray(time_s,float),np.asarray(radius_m,float),np.asarray(phase_rad,float))
        if not np.all(np.isfinite(t+a+phase)) or np.any(a <= WGS84_A):
            raise ValueError('Finite above-Earth position required')
        u=phase+np.deg2rad(self.initial_phase_deg);inc=np.deg2rad(self.inclination_deg)
        p=a*np.sqrt(np.cos(u)**2+(np.cos(inc)*np.sin(u))**2)
        z=a*np.sin(inc)*np.sin(u)
        lon=np.arctan2(np.cos(inc)*np.sin(u),np.cos(u))+np.deg2rad(self.raan_deg)-self.era0-ERA_RATE*t
        lon=(np.rad2deg(lon)+180)%360-180
        e2=WGS84_F*(2-WGS84_F)
        lat=np.arctan2(z,p*(1-e2))
        for _ in range(7):
            n=WGS84_A/np.sqrt(1-e2*np.sin(lat)**2)
            lat=np.arctan2(z+e2*n*np.sin(lat),p)
        n=WGS84_A/np.sqrt(1-e2*np.sin(lat)**2)
        # Stable also at the poles; project radial vector onto ellipsoid normal.
        alt=p*np.cos(lat)+z*np.sin(lat)-n*(1-e2*np.sin(lat)**2)
        return lon,np.rad2deg(lat),alt/1000


class SpaceWeather:
    """Frozen CSV of EXACT seven-component MSIS drivers every three hours."""
    def __init__(self,path,epoch_utc):
        with Path(path).open(encoding='utf-8-sig') as f:rows=list(csv.DictReader(f))
        dates=np.asarray([r['utc'] for r in rows],dtype='datetime64[s]')
        self.seconds=np.asarray((dates-np.datetime64(epoch_utc,'s'))/np.timedelta64(1,'s'),float)
        self.values=np.asarray([[float(r[k]) for k in ['f107_previous_day','f107_center81']+[f'ap{i}' for i in range(7)]] for r in rows])
        if len(rows)<2 or np.any(np.diff(self.seconds)!=10800) or not np.all(np.isfinite(self.values)) or np.any(self.values<0):
            raise ValueError('Incomplete or nonfinite frozen space weather')
        if np.any(self.values[:,:2]<=0):raise ValueError('Positive solar flux required')

    def at(self,time_s):
        t=np.atleast_1d(time_s).astype(float)
        if np.any(t<self.seconds[0]) or np.any(t>self.seconds[-1]):raise ValueError('Weather outside frozen window')
        return self.values[np.searchsorted(self.seconds,t,side='right')-1]


class NRLMSIS2:
    def __init__(self,geometry,weather):
        import pymsis
        from pymsis import msis
        if pymsis.__version__!='0.12.0':raise RuntimeError('This adapter is tested with pymsis 0.12.0')
        # The Windows Fortran wrapper truncates long parameter directories. Use
        # an existing short relative path; never modify package/model files.
        path=Path(pymsis.__file__).resolve().parent
        relative=os.path.relpath(path,Path.cwd()).replace('\\','/')+'/'
        if len(relative.encode('ascii'))>=100 or not (path/'msis2.0.parm').is_file():
            raise RuntimeError('NRLMSIS parameter path needs a shorter working directory')
        msis._MSIS_PARAMETER_PATH=relative
        self.geometry,self.weather=geometry,weather
        self.points_evaluated=0

    def evaluate(self,time_s,radius_m,phase_rad,*,drivers=None):
        import pymsis
        t,a,u=np.broadcast_arrays(np.asarray(time_s,float),np.asarray(radius_m,float),np.asarray(phase_rad,float))
        lon,lat,alt=self.geometry.coordinates(t,a,u)
        dates=np.datetime64(self.geometry.epoch_utc,'us')+np.rint(t.ravel()*1e6).astype('timedelta64[us]')
        w=self.weather.at(t.ravel()) if drivers is None else np.broadcast_to(drivers,t.shape+(9,)).reshape(-1,9)
        rho=pymsis.calculate(dates,lon.ravel(),lat.ravel(),alt.ravel(),w[:,0],w[:,1],w[:,2:],
            version='2.0',geomagnetic_activity=-1)[:,pymsis.Variable.MASS_DENSITY].astype(float).reshape(t.shape)
        self.points_evaluated+=t.size
        if np.any(rho<=0) or not np.all(np.isfinite(rho)):raise FloatingPointError('Invalid kg/m3 density')
        return rho

    __call__=evaluate


class TrajectoryDensityCache:
    """Cubic log-density along short ballistic forecast + radial correction.

    Reanchor to each satellite's actual state at every time block. A circular
    constant-a forecast is used ONLY to tabulate the density, never to propagate
    physical state. Actual a supplies first-order log-density radial correction.
    Four equally spaced nodes determine cubic time interpolation. No cache
    block crosses three-hour driver jumps; both sides use their own drivers.
    """
    def __init__(self,backend,interval_s,mu_m3_s2):
        if interval_s not in (10,30,60,120):raise ValueError('Use a validated cache candidate')
        self.backend,self.interval,self.mu=backend,float(interval_s),mu_m3_s2
        self.end=-np.inf;self.start=-np.inf;self.blocks=0
        self.maximum_radius_offset_m=0.;self.maximum_phase_offset_rad=0.

    def __call__(self,time_s,radius_m,phase_rad):
        a,u=np.asarray(radius_m,float),np.asarray(phase_rad,float)
        if time_s < self.start-1e-9:raise ValueError('Cache requires monotone RK stages')
        if time_s>=self.end-1e-9:
            self.start=float(time_s);self.end=(np.floor(time_s/self.interval)+1)*self.interval
            self.end=min(self.end,(np.floor(time_s/10800)+1)*10800)
            self.a0,self.u0=a.copy(),u.copy();self.n=np.sqrt(self.mu/a**3)
            length=self.end-self.start
            offsets=np.rint(np.linspace(0,length,4))
            # pymsis 0.12 passes integer UTC seconds and integer day-of-year.
            # Keep the last pre-midnight node on that UTC day, then reanchor
            # on the next day rather than smearing its DOY step backwards.
            if self.end % 86400 == 0:
                offsets[-1]=length-1
            if len(np.unique(offsets))<4:
                self.end=-np.inf
                return self.backend.evaluate(time_s,a,u)
            ts=self.start+offsets[:,None]+np.zeros_like(a)[None,:]
            aa=np.broadcast_to(a,ts.shape).copy();uu=u+self.n*offsets[:,None]
            drivers=self.backend.weather.at([self.start])[0]
            # One batched Fortran call: four time nodes, two radial probes.
            rt=np.concatenate([ts.ravel(),np.full(a.size*2,self.start)])
            ra=np.concatenate([aa.ravel(),a-100,a+100])
            ru=np.concatenate([uu.ravel(),u,u])
            log=np.log(self.backend.evaluate(rt,ra,ru,drivers=drivers))
            self.coeff=np.linalg.solve(np.vander(offsets/length,4,increasing=True),log[:4*a.size].reshape(4,a.size))
            self.radial_gradient=(log[5*a.size:]-log[4*a.size:5*a.size])/200
            self.blocks+=1
        dt=time_s-self.start;z=dt/(self.end-self.start)
        da=a-self.a0;du=u-(self.u0+self.n*dt)
        self.maximum_radius_offset_m=max(self.maximum_radius_offset_m,float(np.max(np.abs(da))))
        self.maximum_phase_offset_rad=max(self.maximum_phase_offset_rad,float(np.max(np.abs(du))))
        if np.max(np.abs(da))>100 or np.max(np.abs(du))>1e-6:
            raise FloatingPointError('Cache trajectory moved outside validated spatial neighborhood')
        return np.exp(self.coeff[0]+z*(self.coeff[1]+z*(self.coeff[2]+z*self.coeff[3]))+self.radial_gradient*da)
