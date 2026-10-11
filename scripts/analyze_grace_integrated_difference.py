"""Coarse GRACE A/B integrated density on quality-controlled observed support.

No 21-satellite field reconstruction. Cached raw data and weather are reused;
three one-day windows are selected by weather before density is evaluated.
This method is conditional on model-based common-altitude normalization.
"""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from sat_string.environment_integrals import supported_integral
from sat_string.nrlmsis import SpaceWeather
from sat_string.observed_density import evaluate_nrlmsis_geographic
from sat_string.io import prepare_output_directory,write_json
from compare_grace_density import read_selected,positions,interpolate_good
from analyze_observed_density import write_csv

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/real_environment_20261011/grace_integrals'
DOC=ROOT/'docs/real_environment_applicability/data'


def select_windows(path):
    """Weather-only first eligible24h within the already frozen driver support.

    No density/result read. These are method examples, not climatology or
    independent storm samples; the disturbed day is storm recovery.
    """
    data=np.loadtxt(path)
    times=np.array([np.datetime64(f'{int(r[0])}-01-01')+np.timedelta64(int(r[1])-1,'D')+np.timedelta64(int(r[2]),'h') for r in data])
    candidates=[]
    for day in np.arange(np.datetime64('2003-10-25'),np.datetime64('2003-11-23')):
        use=(times>=day)&(times<day+np.timedelta64(1,'D'));a=data[use]
        if len(a)!=24 or np.any(a[:,40]==99999) or np.any(a[:,38]==99):continue
        kp=float(a[:,38].max()/10);dst=float(a[:,40].min())
        group='severe' if dst<=-250 else ('disturbed' if kp>=5 and dst<=-50 else ('quiet' if kp<5 and dst>-50 else 'unclassified'))
        candidates.append(dict(name=group,start_utc=str(day)+'T00:00:00',end_utc=str(day+np.timedelta64(1,'D'))+'T00:00:00',
             max_kp=kp,min_dst_nT=dst,mean_f107_adjusted=float(a[:,50].mean())))
    windows=[next(w for w in candidates if w['name']==g) for g in ('quiet','disturbed','severe')]
    write_json(DOC/'grace_window_selection.json',dict(selected=windows,candidates=candidates,
      rule='First eligible24h in frozen2003-10-25..11-23 support: quiet Kp<5/Dst>-50; disturbed Kp>=5/-250<Dst<=-50; severe Dst<=-250; no density read',
      qualifications='Geomagnetic quiet only, not low solar activity. Disturbed sample is recovery from Halloween storm. No population quantiles.',
      source='https://spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni/omni2_2003.dat'))
    return windows


def main():
    prepare_output_directory(OUT);DOC.mkdir(parents=True,exist_ok=True)
    raw=ROOT/'outputs/observed_density_20261010/raw'
    drivers=ROOT/'outputs/observed_density_20261010/analysis_v2/drivers.csv'
    weather=SpaceWeather(drivers,'2003-10-25T00:00:00')
    windows=select_windows(raw/'omni2_2003.dat')
    summaries=[];segments=[]
    for window in windows:
        series=[]
        for mission in ('A','B'):
            a,dates,good=read_selected(raw/f'Density_grace{mission}_3deg_2003_v2.3.ascii.tar.gz',window['start_utc'][:10],window['end_utc'][:10])
            background=np.full(len(a),np.nan);bg400=background.copy()
            background[good]=evaluate_nrlmsis_geographic(dates[good],a[good,5],a[good,4],a[good,6],weather)
            bg400[good]=evaluate_nrlmsis_geographic(dates[good],a[good,5],a[good,4],np.full(good.sum(),400.),weather)
            factor=bg400/background
            series.append((a,dates,good,bg400,a[:,11]*factor,a[:,15]*factor,
                           np.asarray((dates-np.datetime64('1970-01-01'))/np.timedelta64(1,'s'),float),positions(a)))
        a,dates,good,bga,ra,ea,t,pos=series[0]
        b,bd,bgood,bgb,rb,eb,bt,bpos=series[1]
        rb_sim=interpolate_good(bt,rb,bgood,t);eb_sim=interpolate_good(bt,eb,bgood,t)
        bgb_sim=interpolate_good(bt,bgb,bgood,t);psim=interpolate_good(bt,bpos,bgood,t)
        velocity=interpolate_good(bt,np.gradient(bpos,bt,axis=0),bgood,t)
        lag=np.sum((pos-psim)*velocity,axis=1)/np.sum(velocity**2,axis=1)
        rb_match=interpolate_good(bt,rb,bgood,t+lag);eb_match=interpolate_good(bt,eb,bgood,t+lag)
        pmatch=interpolate_good(bt,bpos,bgood,t+lag);mismatch=np.linalg.norm(pos-pmatch,axis=1)
        valid=good&np.isfinite(rb_sim)&np.isfinite(rb_match)&(np.abs(lag)<60)&(mismatch<50)
        # Conservative removal at driver boundaries avoids integrating a jump
        # introduced by the model-based height normalization across missing bins.
        valid[1:] &= np.floor(t[1:]/10800)==np.floor(t[:-1]/10800)
        observed=ra-rb_sim;background=bga-bgb_sim;matched=ra-rb_match
        common_integral=supported_integral(t,.5*(ra+rb_sim),valid,100)['integral_observed_support']
        result=dict(window=window['name'],rows=len(t),paired_rows=int(valid.sum()),
          median_separation_km=float(np.nanmedian(np.linalg.norm(pos[valid]-psim[valid],axis=1))),
          median_time_shift_s=float(np.median(lag[valid])),median_altitude_a_km=float(np.median(a[good,6])),
          median_altitude_b_km=float(np.median(b[bgood,6])),
          simultaneous_abs_difference_p95_kg_m3=float(np.percentile(np.abs(observed[valid]),95)),
          matched_abs_difference_p95_kg_m3=float(np.percentile(np.abs(matched[valid]),95)),
          modeled_abs_difference_p95_kg_m3=float(np.percentile(np.abs(background[valid]),95)))
        for label,values,error in [('simultaneous_obs',observed,ea+eb_sim),('matched_obs',matched,ea+eb_match),
                                  ('simultaneous_model',background,None),('obs_minus_model',observed-background,ea+eb_sim)]:
            info=supported_integral(t,values,valid,100,error)
            for n,segment in enumerate(info.pop('segments'),1):
                segments.append(dict(window=window['name'],comparison=label,segment=n,**segment))
            result[label]=info
        result['constant_relative_calibration_offset_equivalent']=abs(result['simultaneous_obs']['integral_observed_support'])/common_integral
        result['calibration_offset_definition']='Magnitude of constant A/B relative density offset that would produce this net integral on the same observed support; no bias removed'
        summaries.append(result)
        np.savez_compressed(OUT/f'{window["name"]}.npz',utc=dates,time_s=t,valid=valid,
                            observed_difference_400=observed,modeled_difference_400=background,matched_difference_400=matched,
                            difference_error_envelope_400=ea+eb_sim,separation_km=np.linalg.norm(pos-psim,axis=1),lag_s=lag)
        print(window['name'],result['simultaneous_obs'],flush=True)
    write_csv(DOC/'grace_segments.csv',segments)
    write_json(DOC/'grace_integrated_summary.json',dict(results=summaries,
      source='https://zenodo.org/records/13287614',product='Sutton v2.3, three-degree averages; NOT100km spatial resolution',
      sign='rho400_A-rho400_B; simultaneous or local tangent time-match',
      normalization='rho400=rho_obs*NRLMSIS2(400)/NRLMSIS2(observation height); preserve model-relative discrepancy; assumption altitude-independent discrepancy',
      support='only consecutive valid paired samples within100s; no rejected-bin/gap bridging; no full24h or60day force reconstruction',
      uncertainty='integrated sum of supplied individual density errors, fully correlated envelope; not confidence interval; missing-data/NRLMSIS errors not included',
      event_selection='First eligible weather-only quiet/disturbed/severe24h in cached driver support; quiet only geomagnetic; disturbed is storm recovery; no inferred climatology',
      raw_sha256={p.name:hashlib.file_digest(p.open('rb'),'sha256').hexdigest() for p in raw.glob('Density_grace*_3deg_2003_v2.3.ascii.tar.gz')},
      weather_sha256=hashlib.sha256(drivers.read_bytes()).hexdigest(),new_fleet_simulation=False))


if __name__=='__main__':main()
