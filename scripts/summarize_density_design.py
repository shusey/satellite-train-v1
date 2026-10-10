"""Reuse density arrays; select track-based candidates and four Japanese figures."""
import argparse
import csv
import json
from pathlib import Path
import sys
import numpy as np
import yaml
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from sat_string.observed_density import gaussian,separation_differences
from analyze_observed_density import extract_events,write_csv,summary


def gaussian_d100(amplitude,sigma_km):
    s=np.linspace(-6*sigma_km,6*sigma_km,12001)
    g=lambda z:amplitude*np.exp(-.5*(z/sigma_km)**2)
    return float(np.max(np.abs(g(s+100)-g(s))))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--analysis',type=Path,required=True)
    p.add_argument('--grace-ab',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--config',type=Path,default=Path('configs/observed_density_design.yaml'));args=p.parse_args()
    if args.output.exists():raise FileExistsError(args.output)
    args.output.mkdir(parents=True);out=args.output
    config=yaml.safe_load(args.config.read_text(encoding='utf-8'));events=[];tracks={};physical=[]
    for window in config['windows']:
        name=window['name']
        for mission in ['champ','grace']:
            d=dict(np.load(args.analysis/f'{mission}_{name}.npz'));tracks[(mission,name)]=d
            # Refit only: no new NRLMSIS calls and no duplicate density calculation.
            ev=extract_events(d,name,mission,config['local_event_method'])
            for e in ev:
                e['fit_d100_fraction_local_bg']=gaussian_d100(e['amplitude_fraction_local_bg'],e['sigma_km']) if e['amplitude_fraction_local_bg'] is not None else None
            events.extend(ev)
            diff,_=separation_differences(d['s_km'],d['rho400'],d['utc'],d['valid'],100)
            relative=2*np.abs(diff)/(2*d['rho400']+diff)
            physical.append(dict(window=name,mission=mission,**summary(relative),
                definition='abs(rho400(s+100)-rho400(s))/mean(two rho400); lagged not simultaneous'))
    write_csv(out/'local_events_refined.csv',events);write_csv(out/'rho400_pair_relative_difference.csv',physical)
    candidates=[e for e in events if e['mission']=='champ' and e['accepted'] and e['sign']==1]
    context=[e for e in candidates if e['window']=='prestorm_october']
    nominal=min(context,key=lambda e:abs(e['fit_d100_fraction_local_bg']-np.median([x['fit_d100_fraction_local_bg'] for x in context])))
    # Select a relatively narrow, well fitted Halloween pulse without manual timing.
    storm=[e for e in candidates if e['window']=='halloween' and e['r_squared']>=.9 and e['fwhm_km']<=600]
    disturbed=max(storm,key=lambda e:e['fit_d100_fraction_local_bg'])
    severe=max(candidates,key=lambda e:e['fit_d100_fraction_local_bg'])
    chosen=[];segments=[]
    for case,e in zip(['context_representative','strong_observed','largest_fitted_D100'],[nominal,disturbed,severe]):
        d=tracks[('champ',e['window'])];i=np.argmin(abs(d['s_km']-e['center_km']))
        regional_background=float(d['rho_bg400'][i]*(1+e['offset_delta']))
        chosen.append(dict(case=case,**e,regional_background_400_kg_m3=regional_background,
            fitted_d100_density_kg_m3=regional_background*e['fit_d100_fraction_local_bg'],
            lifetime_status='not_identified',duration_s=None,propagation_speed_status='requires frame/projection/lifetime validation',
            runnable_as_observation_validated_disturbance=False))
        # Include only contiguous, QC-approved data used by the ±1500 km fit.
        from sat_string.observed_density import contiguous_segments
        ix=next(ix for ix in contiguous_segments(d['valid'],d['utc']) if i in ix)
        use=ix[np.abs(d['s_km'][ix]-e['center_km'])<=1500]
        for j in use:
            segments.append(dict(case=case,utc=str(d['utc'][j]),distance_from_fit_center_km=float(d['s_km'][j]-e['center_km']),
                delta=float(d['delta'][j]),rho_obs_kg_m3=float(d['density'][j]),rho_background_kg_m3=float(d['rho_bg'][j]),
                rho400_kg_m3=float(d['rho400'][j]),altitude_km=float(d['alt'][j]/1000),latitude_deg=float(d['lat'][j]),longitude_deg=float(d['lon'][j])))
    (out/'design_candidates.json').write_text(json.dumps(chosen,ensure_ascii=False,indent=2),encoding='utf-8')
    write_csv(out/'representative_segments.csv',segments)
    weather=list(csv.DictReader((args.analysis/'space_weather_hourly.csv').open()))
    for e in events:
        w=next(w for w in weather if w['utc'][:13]==e['utc'][:13] and w['window']==e['window'])
        e.update(kp=float(w['kp']),ap_3h=float(w['ap_3h']),dst_nT=float(w['dst_nT']),ae_nT=float(w['ae_nT']))
    write_csv(out/'local_events_with_weather.csv',events)
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'Meiryo','axes.unicode_minus':False,'font.size':10})
    names={'prestorm_october':'10月25–27日（事前比較）','halloween':'10月29–31日','november':'11月20–22日'}
    colors=['#3675a9','#d08b14','#ae3450']
    fig,ax=plt.subplots(3,1,figsize=(10,8),layout='constrained')
    for row,(name,title) in enumerate(names.items()):
        d=tracks[('champ',name)];dt=(d['utc']-d['utc'][0])/np.timedelta64(1,'h');y=d['delta'].copy();y[~d['valid']]=np.nan
        ax[row].plot(dt,y*100,lw=.6,label='CHAMP：観測／NRLMSIS 2.0 − 1')
        ax[row].set(title=f'2003年 {title}',ylabel='背景モデル比の偏差 [%]',xlabel='区間開始からの時間 [h]');ax[row].grid(alpha=.3)
    ax[0].legend(fontsize=9);fig.savefig(out/'01_density_residual.png',dpi=160);plt.close(fig)
    rows=list(csv.DictReader((args.analysis/'scale_statistics.csv').open()))
    fig,ax=plt.subplots(figsize=(8,5),layout='constrained')
    for name,c in zip(names,colors):
        r=[r for r in rows if r['mission']=='champ' and r['window']==name]
        x=[float(r['separation_km']) for r in r];ax.plot(x,[100*float(r['p95']) for r in r],'o-',color=c,label=names[name])
    ax.axvspan(0,150,color='gray',alpha=.12,label='少数サンプルでの推定（解像度に注意）')
    ax.set(xlabel='沿軌道の距離差 [km]',ylabel='偏差の差の95百分位 [ポイント]',title='CHAMP単一軌道：短い時間差を含む差分（同時観測ではない）',xlim=(0,1050))
    ax.legend(fontsize=9);ax.grid(alpha=.3);fig.savefig(out/'02_scale_differences.png',dpi=160);plt.close(fig)
    fig,ax=plt.subplots(3,1,figsize=(9,9),layout='constrained')
    labels=['事前比較の代表候補','強い観測候補','100 km差が大きいGaussian候補']
    for a,e,label in zip(ax,chosen,labels):
        segment=[r for r in segments if r['case']==e['case']];s=np.array([r['distance_from_fit_center_km'] for r in segment]);y=np.array([r['delta'] for r in segment])
        a.plot(s,100*y,'o-',ms=3,lw=.6,label='観測偏差');a.plot(s,100*gaussian(s,e['amplitude_delta'],0,e['sigma_km'],e['offset_delta'],e['slope_delta_per_km']),label='Gaussian＋地域的背景')
        a.plot(s,100*(e['offset_delta']+s*e['slope_delta_per_km']),'--',label='地域的背景')
        a.set(title=f"{label}：{e['utc']} UTC／幅 {e['fwhm_km']:.0f} km／決定係数 {e['r_squared']:.3f}",xlabel='ピーク中心からの沿軌道距離 [km]',ylabel='背景モデル比の偏差 [%]');a.grid(alpha=.3)
    ax[0].legend(fontsize=9);fig.savefig(out/'03_gaussian_candidates.png',dpi=160);plt.close(fig)
    stats=list(csv.DictReader((args.grace_ab/'grace_ab_statistics.csv').open()))
    fig,ax=plt.subplots(figsize=(9,5),layout='constrained');x=np.arange(3)
    ax.bar(x-.17,[100*float(r['simultaneous_p95_abs']) for r in stats],.34,label='同時刻A/B（空間差を含む）')
    ax.bar(x+.17,[100*float(r['matched_p95_abs']) for r in stats],.34,label='約28–30秒ずらした近傍位置比較')
    ax.set(xticks=x,xticklabels=[names[r['window']] for r in stats],ylabel='偏差の差の95百分位 [ポイント]',title='GRACE A/B：約216–226 km間隔、緯度3度の平均データ')
    ax.legend(fontsize=9);ax.grid(axis='y',alpha=.3);fig.savefig(out/'04_grace_comparison.png',dpi=160);plt.close(fig)
    diagnostics=dict(total_fit_candidates=len(events),accepted=sum(e['accepted'] for e in events),
        champ_accepted=sum(e['accepted'] and e['mission']=='champ' for e in events),
        champ_positive_candidates=len(candidates),candidate_selection='context median fitted D100; documented narrow Halloween pulse; maximum fitted D100 across selected windows',
        artificial_gaussian_d100=gaussian_d100(.5,100/(2*np.sqrt(2*np.log(2)))),
        cached_input_directory=str(args.analysis),satellite_simulations_executed=0)
    (out/'assessment.json').write_text(json.dumps(diagnostics,indent=2),encoding='utf-8')
    print(json.dumps(diagnostics,indent=2));print([(e['case'],e['utc'],e['amplitude_fraction_local_bg'],e['fwhm_km'],e['fit_d100_fraction_local_bg']) for e in chosen])


if __name__=='__main__':main()
