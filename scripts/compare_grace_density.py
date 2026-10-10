"""Sutton v2.3 A/B, coarse three-degree bins: NOT 100-km resolved data."""
import argparse
import io
import json
import tarfile
from pathlib import Path
import sys
import numpy as np
import yaml
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from sat_string.observed_density import gps_2003_to_utc,evaluate_nrlmsis_geographic
from sat_string.nrlmsis import SpaceWeather
from analyze_observed_density import write_csv,summary


def read_selected(path,start,end):
    rows=[]
    with tarfile.open(path) as archive:
        for member in archive:
            if not member.isfile():continue
            day=int(member.name.split('_')[-1].split('.')[0])
            date=np.datetime64('2003-01-01')+np.timedelta64(day-1,'D')
            if not np.datetime64(start)-np.timedelta64(1,'D')<=date<np.datetime64(end)+np.timedelta64(1,'D'):continue
            rows.append(np.loadtxt(io.BytesIO(archive.extractfile(member).read()),skiprows=2))
    a=np.concatenate(rows)
    gps=np.datetime64('2003-01-01','us')+(np.rint(((a[:,1]-1)*86400+a[:,2])*1e6)).astype('timedelta64[us]')
    dates=gps_2003_to_utc(gps)
    use=(dates>=np.datetime64(start))&(dates<np.datetime64(end))
    a=a[use];dates=dates[use]
    # Column 18 differs in release PDF vs actual header: exclude all nonzero.
    good=np.all(np.isfinite(a),axis=1)&(a[:,11]>0)&(a[:,15]>0)&(a[:,16]>0)&(a[:,17]==0)
    return a,dates,good


def positions(a):
    lat=np.deg2rad(a[:,4]);lon=np.deg2rad(a[:,5]);h=a[:,6];f=1/298.257223563;e2=f*(2-f)
    n=6378.137/np.sqrt(1-e2*np.sin(lat)**2)
    return np.column_stack([(n+h)*np.cos(lat)*np.cos(lon),(n+h)*np.cos(lat)*np.sin(lon),(n*(1-e2)+h)*np.sin(lat)])


def interpolate_good(seconds,values,valid,target):
    """No rejected bin, >100-s gap, extrapolation or 3-hour driver jump."""
    j=np.searchsorted(seconds,target,side='right')-1;k=j+1
    inside=(j>=0)&(k<len(seconds));j=np.clip(j,0,len(seconds)-1);k=np.clip(k,0,len(seconds)-1)
    keep=inside&valid[j]&valid[k]&(seconds[k]-seconds[j]<=100)&(np.floor(seconds[j]/10800)==np.floor(seconds[k]/10800))
    f=np.divide(target-seconds[j],seconds[k]-seconds[j],out=np.zeros(len(target)),where=seconds[k]!=seconds[j])
    if values.ndim>1:f=f[:,None]
    result=values[j]+f*(values[k]-values[j]);result[~keep]=np.nan
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--raw-dir',type=Path,required=True)
    p.add_argument('--drivers',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--config',type=Path,default=Path('configs/observed_density_design.yaml'));args=p.parse_args()
    if args.output.exists():raise FileExistsError(args.output)
    args.output.mkdir(parents=True)
    weather=SpaceWeather(args.drivers,'2003-10-25T00:00:00');results=[];samples=[]
    windows=yaml.safe_load(args.config.read_text(encoding='utf-8'))['windows']
    for window in windows:
        series=[]
        for mission in ['A','B']:
            a,dates,good=read_selected(args.raw_dir/f'Density_grace{mission}_3deg_2003_v2.3.ascii.tar.gz',window['start_utc'][:10],window['end_utc'][:10])
            bg=np.full(len(a),np.nan);bg[good]=evaluate_nrlmsis_geographic(dates[good],a[good,5],a[good,4],a[good,6],weather)
            delta=a[:,11]/bg-1;error=a[:,15]/bg
            seconds=(dates-np.datetime64('1970-01-01'))/np.timedelta64(1,'s')
            series.append((a,dates,good,bg,delta,error,seconds,positions(a)))
        a,dates,good,bg,delta,error,t,pos=series[0];b,bd,bgood,bbg,bdelta,berror,bt,bpos=series[1]
        bsim=interpolate_good(bt,bdelta,bgood,t);psim=interpolate_good(bt,bpos,bgood,t)
        separation=np.linalg.norm(pos-psim,axis=1)
        # Project along local B trajectory to estimate nearest geographic passage.
        velocity=interpolate_good(bt,np.gradient(bpos,bt,axis=0),bgood,t)
        lag=np.sum((pos-psim)*velocity,axis=1)/np.sum(velocity**2,axis=1)
        bmatch=interpolate_good(bt,bdelta,bgood,t+lag)
        pmatch=interpolate_good(bt,bpos,bgood,t+lag)
        mismatch=np.linalg.norm(pos-pmatch,axis=1)
        valid=good&np.isfinite(bsim)&np.isfinite(bmatch)&(np.abs(lag)<60)&(mismatch<50)
        # Summaries are complementary coarse-bin checks, not pure space/time signals.
        simultaneous=delta-bsim;matched=delta-bmatch
        results.append(dict(window=window['name'],a_rows=len(a),b_rows=len(b),pairs=int(np.sum(valid)),
            median_separation_km=summary(separation[valid])['median'],median_time_shift_s=summary(lag[valid])['median'],
            median_matched_position_error_km=summary(mismatch[valid])['median'],
            simultaneous_median_abs=summary(np.abs(simultaneous[valid]))['median'],simultaneous_p95_abs=summary(np.abs(simultaneous[valid]))['p95'],
            simultaneous_max_abs=summary(np.abs(simultaneous[valid]))['max'],
            matched_median_abs=summary(np.abs(matched[valid]))['median'],matched_p95_abs=summary(np.abs(matched[valid]))['p95'],
            median_relative_density_uncertainty_a=summary(a[good,15]/a[good,11])['median'],
            median_relative_density_uncertainty_b=summary(b[bgood,15]/b[bgood,11])['median']))
        # Full selected pairs remain local NPZ; compact CSV contains every tenth pair.
        np.savez_compressed(args.output/f"pairs_{window['name']}.npz",utc=dates,valid=valid,delta_a=delta,delta_b_same_time=bsim,
            delta_b_matched=bmatch,separation_km=separation,lag_s=lag,position_mismatch_km=mismatch)
        for i in np.flatnonzero(valid)[::10]:
            samples.append(dict(window=window['name'],utc=str(dates[i]),delta_a=float(delta[i]),delta_b_same_time=float(bsim[i]),
                delta_b_matched=float(bmatch[i]),separation_km=float(separation[i]),lag_s=float(lag[i]),position_mismatch_km=float(mismatch[i])))
        print(results[-1],flush=True)
    write_csv(args.output/'grace_ab_statistics.csv',results);write_csv(args.output/'grace_ab_sample.csv',samples)
    (args.output/'method.json').write_text(json.dumps(dict(product='Sutton v2.3 3-degree geographic bins',background='NRLMSIS 2.0',
        gps_minus_utc_s=13,max_interpolation_gap_s=100,max_matching_error_km=50,
        same_location_method='local tangent projection; approximate geographic match, no complete separation of space/time',
        resolution_warning='300-350 km averaging; not direct evidence for 100-km structures'),indent=2),encoding='utf-8')


if __name__=='__main__':main()
