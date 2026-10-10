"""Selected CHAMP/GRACE events: NRLMSIS2 residuals, scales, local fits.

No fleet simulation or physical/control changes. Raw and full processed data
stay local; compact statistics and representative segments can be published.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
import yaml
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from sat_string.nrlmsis import SpaceWeather
from sat_string.observed_density import (read_igg,normalize_to_400,evaluate_nrlmsis_geographic,
    contiguous_segments,separation_differences,fit_gaussian_event)


def write_csv(path,rows):
    if path.exists():raise FileExistsError(path)
    if not rows:return
    with path.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def freeze_weather(raw,out):
    with raw.open(encoding='utf-8-sig') as f:
        reader=csv.DictReader(f);columns=reader.fieldnames;rows={r['DATE']:r for r in reader}
    epoch=np.datetime64('2003-10-25');end=np.datetime64('2003-11-23')
    support=[]
    for day in np.arange(epoch-np.timedelta64(40,'D'),end+np.timedelta64(41,'D')):
        r=rows[str(day)]
        if r['F10.7_DATA_TYPE']!='OBS':raise ValueError('Unobserved flux: no filling')
        support.append(r)
    write_csv(out/'observed_support.csv',support)
    def ap_at(t):
        day=t.astype('datetime64[D]');slot=int((t-day)/np.timedelta64(3,'h'))+1
        return float(rows[str(day)][f'AP{slot}'])
    drivers=[]
    for t in np.arange(epoch.astype('datetime64[h]'),end.astype('datetime64[h]')+np.timedelta64(1,'h'),np.timedelta64(3,'h')):
        day=t.astype('datetime64[D]');r=rows[str(day)]
        history=[float(r['AP_AVG'])]+[ap_at(t-np.timedelta64(h,'h')) for h in [0,3,6,9]]
        history += [float(np.mean([ap_at(t-np.timedelta64(h,'h')) for h in range(12,34,3)])),float(np.mean([ap_at(t-np.timedelta64(h,'h')) for h in range(36,58,3)]))]
        drivers.append(dict(utc=str(t.astype('datetime64[s]')),f107_previous_day=float(rows[str(day-np.timedelta64(1,'D'))]['F10.7_OBS']),
            f107_center81=float(r['F10.7_OBS_CENTER81']),**{f'ap{i}':v for i,v in enumerate(history)}))
    write_csv(out/'drivers.csv',drivers)
    import pymsis.utils as utils
    utils._F107_AP_PATH=out/'observed_support.csv';utils._DATA=None
    # November 4 radio-burst flux (560.9 sfu) is outside every selected window.
    # pymsis' automatic reader substitutes its 81-day mean for flux >400.
    # Our explicit OBS inputs remain unchanged; verify only actually used dates.
    selected=[r for r in drivers if r['utc'][:10] in
        {str(d) for first,last in [('2003-10-25','2003-11-01'),('2003-11-20','2003-11-23')]
         for d in np.arange(np.datetime64(first),np.datetime64(last))}]
    a,b,c=utils.get_f107_ap(np.array([r['utc'] for r in selected],dtype='datetime64[s]'))
    ours=np.array([[r['f107_previous_day'],r['f107_center81']]+[r[f'ap{i}'] for i in range(7)] for r in selected])
    if not np.array_equal(ours,np.column_stack([a,b,c])):raise ValueError('Driver history differs from pymsis')
    return SpaceWeather(out/'drivers.csv','2003-10-25T00:00:00')


def summary(x):
    x=np.asarray(x);x=x[np.isfinite(x)]
    if not len(x):return dict(count=0,median=None,p90=None,p95=None,p99=None,max=None)
    return dict(count=len(x),median=float(np.median(x)),p90=float(np.percentile(x,90)),p95=float(np.percentile(x,95)),p99=float(np.percentile(x,99)),max=float(np.max(x)))


def extract_events(data,window,mission,settings):
    from scipy.ndimage import median_filter
    from scipy.signal import find_peaks
    s=data['s_km'];delta=data['delta'];events=[]
    for ix in contiguous_segments(data['valid'],data['utc']):
        if s[ix[-1]]-s[ix[0]]<3000:continue
        ds=float(np.median(np.diff(s[ix])))
        width=int(round(settings['background_window_km']/ds));width+=1-width%2
        local=delta[ix]-median_filter(delta[ix],size=width,mode='nearest')
        for sign in [1,-1]:
            peaks,_=find_peaks(sign*local,prominence=settings['minimum_prominence_delta'],distance=max(1,int(settings['minimum_peak_distance_km']/ds)))
            for peak in peaks:
                at=ix[peak];center=s[at];half=settings['fit_window_halfwidth_km']
                if center-s[ix[0]]<half or s[ix[-1]]-center<half:continue
                fitix=ix[(s[ix]>=center-half)&(s[ix]<=center+half)]
                fit=fit_gaussian_event(s[fitix],delta[fitix],center,sign,
                    minimum_amplitude=settings['minimum_fitted_amplitude_delta'],
                    minimum_improvement=settings['minimum_improvement_vs_linear'],
                    minimum_r_squared=settings['acceptable_r_squared'],
                    minimum_fwhm_samples=settings['minimum_resolved_fwhm_samples'])
                # A high R2 alone cannot guarantee a single isolated pulse.
                neighbors,_=find_peaks(sign*local,prominence=max(.1,.25*abs(fit['amplitude_delta'])))
                close=neighbors[np.abs(s[ix[neighbors]]-fit['center_km'])<fit['fwhm_km']]
                multiplet=len(close)>1
                if multiplet:fit['accepted']=False
                fit['shape']='multiple_peak' if multiplet else ('gaussian_candidate' if fit['accepted'] else 'poor_or_unresolved_fit')
                nearest=int(np.argmin(np.abs(s-fit['center_km'])))
                sensitivity=[]
                for wkm in settings['sensitivity_windows_km']:
                    n=int(round(wkm/ds));n+=1-n%2
                    baseline=median_filter(delta[ix],size=n,mode='nearest')
                    sensitivity.append(float(delta[at]-baseline[peak]))
                pair=data['d100'][fitix];pair=pair[np.isfinite(pair)]
                event=dict(window=window,mission=mission,event_number=len(events)+1,utc=str(data['utc'][nearest].astype('datetime64[s]')),
                    sign=sign,peak_delta_raw=float(delta[at]),local_peak_delta=float(local[peak]),
                    local_peak_sensitivity_min=float(min(sensitivity)),local_peak_sensitivity_max=float(max(sensitivity)),
                    background_density_kg_m3=float(data['rho_bg'][nearest]),rho400_kg_m3=float(data['rho400'][nearest]),
                    altitude_km=float(data['alt'][nearest]/1000),latitude_deg=float(data['lat'][nearest]),longitude_deg=float(data['lon'][nearest]),
                    local_solar_time_h=float(data['lst'][nearest]),d100_local_max=float(np.max(np.abs(pair))) if len(pair) else None,
                    apparent_passage_time_s=float((data['utc'][fitix[-1]]-data['utc'][fitix[0]])/np.timedelta64(1,'s')*fit['fwhm_km']/(s[fitix[-1]]-s[fitix[0]])),
                    structure_lifetime_s=None,**fit)
                events.append(event)
    return events


def omni_subset(raw,windows):
    data=np.loadtxt(raw);rows=[]
    for r in data:
        date=np.datetime64(f'{int(r[0])}-01-01T00')+np.timedelta64(int(r[1])-1,'D')+np.timedelta64(int(r[2]),'h')
        name=next((w['name'] for w in windows if np.datetime64(w['start_utc'])<=date<np.datetime64(w['end_utc'])),None)
        if name is None:continue
        rows.append(dict(window=name,utc=str(date.astype('datetime64[s]')),kp=float(r[38]/10) if r[38]!=99 else None,
            dst_nT=float(r[40]) if r[40]!=99999 else None,ae_nT=float(r[41]) if r[41]!=9999 else None,
            ap_3h=float(r[49]) if r[49]!=999 else None,f107_adjusted_sfu=float(r[50]) if r[50]!=999.9 else None))
    return rows


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',type=Path,default=Path('configs/observed_density_design.yaml'))
    p.add_argument('--raw-dir',type=Path,required=True);p.add_argument('--weather-raw',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    if args.output.exists():raise FileExistsError(args.output)
    args.output.mkdir(parents=True);out=args.output
    config=yaml.safe_load(args.config.read_text(encoding='utf-8'))
    weather=freeze_weather(args.weather_raw,out)
    statistics=[];scales=[];events=[];qc=[]
    for window in config['windows']:
        month=window['start_utc'][:7]
        for mission in ['champ','grace']:
            path=args.raw_dir/f'{mission}_TND-IGG_RL01_{month}.nc'
            data,quality=read_igg(path,window['start_utc'],window['end_utc'])
            valid=data['valid'];dates=data['utc'][valid];lat=data['lat'][valid];lon=data['lon'][valid]
            bg=evaluate_nrlmsis_geographic(dates,lon,lat,data['alt'][valid]/1000,weather)
            bg400=evaluate_nrlmsis_geographic(dates,lon,lat,np.full(len(dates),400.),weather)
            delta,rho400=normalize_to_400(data['density'][valid],bg,bg400)
            for key,values in [('rho_bg',bg),('rho_bg400',bg400),('rho400',rho400),('delta',delta)]:
                data[key]=np.full(len(valid),np.nan);data[key][valid]=values
            for scale in config['scales_km']:
                diff,lag=separation_differences(data['s_km'],data['delta'],data['utc'],valid,scale)
                if scale==100:data['d100']=diff;data['lag100_s']=lag
                scales.append(dict(window=window['name'],mission=mission,separation_km=scale,median_lag_s=summary(lag)['median'],
                    **summary(np.abs(diff)),interpretation='single-track lagged estimate; resolution-limited'))
            statistics.append(dict(window=window['name'],mission=mission,delta_min=float(np.min(delta)),delta_max=float(np.max(delta)),
                delta_median=float(np.median(delta)),rho_obs_min=float(np.min(data['density'][valid])),rho_obs_max=float(np.max(data['density'][valid])),
                rho400_min=float(np.min(rho400)),rho400_max=float(np.max(rho400)),d100_max=summary(np.abs(data['d100']))['max']))
            qc.append(dict(window=window['name'],mission=mission,**quality))
            events.extend(extract_events(data,window['name'],mission,config['local_event_method']))
            np.savez_compressed(out/f"{mission}_{window['name']}.npz",**data)
            print(window['name'],mission,quality['accepted'],'accepted; max delta',np.max(delta),flush=True)
    write_csv(out/'window_statistics.csv',statistics);write_csv(out/'scale_statistics.csv',scales)
    write_csv(out/'local_events.csv',events);write_csv(out/'quality_control.csv',qc)
    weatherrows=omni_subset(args.raw_dir/'omni2_2003.dat',config['windows'])
    write_csv(out/'space_weather_hourly.csv',weatherrows)
    for event in events:
        t=np.datetime64(event['utc'],'h');r=next(r for r in weatherrows if np.datetime64(r['utc'],'h')==t)
        event.update(kp=r['kp'],dst_nT=r['dst_nT'],ap_3h=r['ap_3h'],ae_nT=r['ae_nT'])
    write_csv(out/'local_events_with_weather.csv',events)
    import pymsis,scipy,netCDF4
    provenance=dict(config=config,pymsis=pymsis.__version__,scipy=scipy.__version__,netCDF4=netCDF4.__version__,
        background_version='2.0',mass_density_unit='kg/m3',gps_minus_utc_s=13,normalization='obs/bg_obs * bg_400',
        source_files={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in args.raw_dir.glob('*.nc')},
        weather_raw_sha256=hashlib.sha256(args.weather_raw.read_bytes()).hexdigest(),no_missing_value_filling=True,
        cadence_is_not_effective_resolution=True,NRLMSIS_driver_boundary_pairs_excluded=True)
    (out/'provenance.json').write_text(json.dumps(provenance,indent=2),encoding='utf-8')
    print('Events',len(events),'accepted',sum(e['accepted'] for e in events),flush=True)


if __name__=='__main__':main()
