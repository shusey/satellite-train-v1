"""Observation postprocessing only; does not change the simulation dynamics.

D_L is a single-track, short-time-lag estimate, not a simultaneous spatial
measurement. All interpolation is limited to contiguous quality-approved data.
"""
import numpy as np


def gps_2003_to_utc(gps_dates):
    dates=np.asarray(gps_dates,dtype='datetime64[us]')
    if np.any(dates<np.datetime64('1999-01-01')) or np.any(dates>=np.datetime64('2006-01-01')):
        raise ValueError('This explicit historical conversion supports 1999-2005 only')
    return dates-np.timedelta64(13,'s')


def read_igg(path,start_utc,end_utc,eclipse_guard_s=30):
    from netCDF4 import Dataset
    with Dataset(path) as ds:
        tv=ds['time'];unit=tv.units
        if not unit.startswith('minutes since ') or tv.long_name!='GPS time':raise ValueError('Unexpected time specification')
        origin=unit.split('since ')[1].split()[0]
        parts=origin.split('-');epoch=np.datetime64(f'{int(parts[0]):04d}-{int(parts[1]):02d}-{int(parts[2]):02d}','us')
        gps=epoch+np.rint(np.asarray(tv[:])*60e6).astype('timedelta64[us]')
        utc=gps_2003_to_utc(gps)
        idx=np.flatnonzero((utc>=np.datetime64(start_utc))&(utc<np.datetime64(end_utc)))
        if ds['density'].units!='kg/m^3' or ds['alt'].units!='m':raise ValueError('Unexpected density/altitude unit')
        data={key:np.asarray(np.ma.filled(ds[key][idx],np.nan),float) for key in ['lat','lon','alt','density','arglat','lst','shadow']}
        valid=np.all(np.column_stack([np.isfinite(v) for v in data.values()]),axis=1)&(data['density']>0)
        qc={}
        for name,v in ds.groups['gap'].variables.items():
            bad=np.ma.filled(v[idx],-1)!=0;qc['gap_'+name]=int(np.sum(bad));valid &= ~bad
        # flag/path is an orbit direction, not a quality flag; never reject it.
        bad=np.ma.filled(ds.groups['flag']['thrust'][idx],-1)!=0
        qc['thrust']=int(np.sum(bad));valid &= ~bad
        shadow=data['shadow'];edge=(shadow>0)&(shadow<1)
        edge[1:] |= np.abs(np.diff(shadow))>0.01
        n=int(np.ceil(eclipse_guard_s/10))
        guarded=np.convolve(edge.astype(int),np.ones(2*n+1,int),mode='same')>0
        qc['eclipse_guard']=int(np.sum(guarded));valid &= ~guarded
        radius=geocentric_radius(data['lat'],data['alt'])
        phase=np.unwrap(np.deg2rad(data['arglat']))
        distance=np.r_[0,np.cumsum(np.diff(phase)*(radius[1:]+radius[:-1])/2/1000)]
        if np.any(np.diff(distance)<=0):raise ValueError('Nonmonotone orbital track')
        data.update(utc=utc[idx],s_km=distance,valid=valid)
        qc.update(rows=len(idx),accepted=int(np.sum(valid)),cadence_s=float(np.median(np.diff(utc[idx])/np.timedelta64(1,'s'))),
            altitude_min_km=float(np.min(data['alt'])/1000),altitude_max_km=float(np.max(data['alt'])/1000))
        return data,qc


def geocentric_radius(latitude_deg,altitude_m):
    # GRS80 ellipsoid recorded in the actual TND-IGG crs metadata.
    lat=np.deg2rad(latitude_deg);a=6378137.;f=1/298.257222101;e2=f*(2-f)
    n=a/np.sqrt(1-e2*np.sin(lat)**2)
    return np.hypot((n+altitude_m)*np.cos(lat),(n*(1-e2)+altitude_m)*np.sin(lat))


def contiguous_segments(valid,dates,max_step_s=11,split_weather=True):
    valid=np.asarray(valid,bool);dates=np.asarray(dates,dtype='datetime64[us]')
    seconds=(dates-np.datetime64('1970-01-01'))/np.timedelta64(1,'s')
    cuts=(np.diff(seconds)>max_step_s)|(np.diff(seconds)<=0)
    if split_weather:cuts |= np.floor(seconds[1:]/10800)!=np.floor(seconds[:-1]/10800)
    indices=np.flatnonzero(valid)
    if not len(indices):return []
    borders=(np.diff(indices)!=1)|cuts[indices[:-1]]
    return [x for x in np.split(indices,np.flatnonzero(borders)+1) if len(x)>1]


def separation_differences(s_km,delta,dates,valid,separation_km):
    """Return signed delta(s+L,t+lag)-delta(s,t); no gap/driver jump crossing."""
    s=np.asarray(s_km,float);d=np.asarray(delta,float)
    difference=np.full(len(s),np.nan);lag=np.full(len(s),np.nan)
    sec=(np.asarray(dates,dtype='datetime64[us]')-np.datetime64('1970-01-01'))/np.timedelta64(1,'s')
    for ix in contiguous_segments(valid,dates):
        target=s[ix]+separation_km;inside=target<=s[ix[-1]];use=ix[inside]
        difference[use]=np.interp(target[inside],s[ix],d[ix])-d[use]
        lag[use]=np.interp(target[inside],s[ix],sec[ix])-sec[use]
    return difference,lag


def gaussian(s,amplitude,center,sigma,offset,slope):
    return offset+slope*(s-center)+amplitude*np.exp(-0.5*((s-center)/sigma)**2)


def fit_gaussian_event(s_km,delta,center_km,sign=1,*,minimum_amplitude=.1,
                       minimum_improvement=.25,minimum_r_squared=.8,minimum_fwhm_samples=4):
    from scipy.optimize import least_squares
    s=np.asarray(s_km,float)-center_km;y=np.asarray(delta,float)
    if len(s)<12 or not np.all(np.isfinite(y)):raise ValueError('Insufficient complete fit data')
    offset=float(np.median(y));amp=sign*max(.1,float(np.max(sign*(y-offset))))
    bounds=([0 if sign>0 else -20,-300,30,-20,-.02], [20 if sign>0 else 0,300,1200,20,.02])
    fit=least_squares(lambda p:gaussian(s,*p)-y,[amp,0,300,offset,0],bounds=bounds,loss='soft_l1',f_scale=.03)
    a,c,sigma,b,m=fit.x;pred=gaussian(s,*fit.x);ss=float(np.sum((y-np.mean(y))**2))
    r2=1-float(np.sum((pred-y)**2))/ss if ss>0 else None
    linear=np.polyval(np.polyfit(s,y,1),s)
    linear_ss=float(np.sum((linear-y)**2))
    pulse_improvement=1-float(np.sum((pred-y)**2))/linear_ss if linear_ss>0 else 0.
    fwhm=float(2*np.sqrt(2*np.log(2))*sigma)
    step=float(np.median(np.diff(s)))
    return dict(amplitude_delta=float(a),center_km=float(center_km+c),sigma_km=float(sigma),fwhm_km=fwhm,
        offset_delta=float(b),slope_delta_per_km=float(m),r_squared=r2,rmse_delta=float(np.sqrt(np.mean((pred-y)**2))),
        improvement_vs_linear=pulse_improvement,
        fwhm_samples=fwhm/step,amplitude_fraction_local_bg=float(a/(1+b)) if b>-1 else None,
        accepted=bool(fit.success and abs(a)>=minimum_amplitude and pulse_improvement>=minimum_improvement and r2 is not None and r2>=minimum_r_squared and fwhm/step>=minimum_fwhm_samples and sigma<1199))


def normalize_to_400(rho_obs,rho_bg_obs,rho_bg_400):
    obs,bg,bg400=np.broadcast_arrays(rho_obs,rho_bg_obs,rho_bg_400)
    if not np.all(np.isfinite(obs+bg+bg400)) or np.any(obs<=0) or np.any(bg<=0) or np.any(bg400<=0):raise ValueError('Finite positive densities required')
    delta=obs/bg-1
    return delta,bg400*(1+delta)


def evaluate_nrlmsis_geographic(dates,lon,lat,alt_km,weather):
    """Explicit 2.0, storm Ap history, observed-date/location/height; kg/m3."""
    import pymsis
    from .nrlmsis import NRLMSIS2,OrbitGeometry
    # Reuse pinned-version/coefficient-path checks; no orbit geometry is used below.
    NRLMSIS2(OrbitGeometry('2003-10-25T00:00:00',0,0,0),weather)
    dates=np.asarray(dates,dtype='datetime64[us]')
    seconds=(dates-np.datetime64('2003-10-25T00:00:00'))/np.timedelta64(1,'s')
    drivers=weather.at(seconds)
    out=pymsis.calculate(dates,np.asarray(lon),np.asarray(lat),np.asarray(alt_km),drivers[:,0],drivers[:,1],drivers[:,2:],
        version='2.0',geomagnetic_activity=-1)
    rho=out[:,pymsis.Variable.MASS_DENSITY].astype(float)
    if not np.all(np.isfinite(rho)) or np.any(rho<=0):raise FloatingPointError('Invalid NRLMSIS density')
    return rho
