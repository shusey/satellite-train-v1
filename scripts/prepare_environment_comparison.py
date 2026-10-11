"""Preregister weather-only selection, then freeze three environments.

No formation output is read here. Original cached CelesTrak snapshot is reused.
Operational classes are study definitions, not official Dst storm categories.
"""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'outputs/real_environment_20261011'


def weather_class(max_kp, storm_days, minimum_dst):
    if minimum_dst<=-250: return 'strong'
    if storm_days>=2 and minimum_dst<=-100: return 'active'
    if max_kp<5 and minimum_dst>-50: return 'quiet'
    return 'unclassified'


def freeze(rows, epoch, folder):
    folder.mkdir(parents=True, exist_ok=False)
    end = epoch + np.timedelta64(60, 'D')
    support = np.arange(epoch-np.timedelta64(40,'D'), end+np.timedelta64(41,'D'))
    selected = [rows[str(d)] for d in support]
    for r in selected:
        if r['F10.7_DATA_TYPE'] != 'OBS' or not 0 < float(r['F10.7_OBS']):
            raise ValueError('Missing/interpolated/predicted solar flux: reject period')
        if not all(np.isfinite(float(r[k])) and float(r[k]) >= 0 for k in
                   ['AP_AVG','F10.7_OBS_CENTER81']+[f'AP{i}' for i in range(1,9)]):
            raise ValueError('Invalid Ap/centered flux: reject period')
    def ap(t):
        day=t.astype('datetime64[D]')
        slot=int((t-day)/np.timedelta64(3,'h'))+1
        return float(rows[str(day)][f'AP{slot}'])
    drivers=[]
    for t in np.arange(epoch.astype('datetime64[h]'),end.astype('datetime64[h]')+np.timedelta64(1,'h'),np.timedelta64(3,'h')):
        day=t.astype('datetime64[D]'); r=rows[str(day)]
        history=[float(r['AP_AVG'])]+[ap(t-np.timedelta64(h,'h')) for h in (0,3,6,9)]
        history += [np.mean([ap(t-np.timedelta64(h,'h')) for h in range(12,34,3)]),
                    np.mean([ap(t-np.timedelta64(h,'h')) for h in range(36,58,3)])]
        drivers.append(dict(utc=str(t.astype('datetime64[s]')),f107_previous_day=float(rows[str(day-np.timedelta64(1,'D'))]['F10.7_OBS']),
                            f107_center81=float(r['F10.7_OBS_CENTER81']),**{f'ap{i}':float(v) for i,v in enumerate(history)}))
    for name,data in [('observed_support.csv',selected),('drivers.csv',drivers)]:
        with (folder/name).open('w',encoding='utf-8',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(data[0])); w.writeheader();w.writerows(data)
    import pymsis.utils as utils
    utils._F107_AP_PATH=folder/'observed_support.csv';utils._DATA=None
    a,b,c=utils.get_f107_ap(np.array([r['utc'] for r in drivers],dtype='datetime64[s]'))
    ours=np.array([[r['f107_previous_day'],r['f107_center81']]+[r[f'ap{i}'] for i in range(7)] for r in drivers])
    if not np.array_equal(ours,np.column_stack([a,b,c])):
        raise ValueError('Weather history differs from pinned wrapper')


def main():
    if (OUT/'protocol.json').exists():raise FileExistsError('Frozen protocol exists: do not overwrite/reselect')
    raw=ROOT/'outputs/nrlmsis2_20261010/SW-All.csv'
    with raw.open(encoding='utf-8-sig') as f: rows={r['DATE']:r for r in csv.DictReader(f)}
    omni=OUT/'raw/omni2_2024.dat';data=np.loadtxt(omni)
    times=np.array([np.datetime64(f'{int(r[0])}-01-01')+np.timedelta64(int(r[1])-1,'D')+np.timedelta64(int(r[2]),'h') for r in data])
    candidates=[]
    for month in ('01','03','05','07','09','11'):
        epoch=np.datetime64('2024-'+month+'-01');end=epoch+np.timedelta64(60,'D')
        days=[rows[str(d)] for d in np.arange(epoch,end)]
        max_kp=max(float(r[f'KP{i}'])/10 for r in days for i in range(1,9))
        storm_days=sum(max(float(r[f'KP{i}'])/10 for i in range(1,9))>=5 for r in days)
        m=(times>=epoch)&(times<end); dst=data[m,40]
        if len(dst)!=1440 or np.any(dst==99999): raise ValueError('Missing hourly Dst')
        minimum=float(dst.min())
        group=weather_class(max_kp,storm_days,minimum)
        support=[rows[str(d)] for d in np.arange(epoch-np.timedelta64(40,'D'),end+np.timedelta64(41,'D'))]
        # The support snapshot extends beyond actual MSIS calls. July30 flux400.7
        # is outside the May run; do not reject May for a never-used daily input.
        eligible=all(r['F10.7_DATA_TYPE']=='OBS' and 0<float(r['F10.7_OBS']) for r in support)
        eligible &= all(0<float(rows[str(d-np.timedelta64(1,'D'))]['F10.7_OBS'])<=400
                        and 0<float(rows[str(d)]['F10.7_OBS_CENTER81'])<=400 for d in np.arange(epoch,end+np.timedelta64(1,'D')))
        candidates.append(dict(start=str(epoch),end_exclusive=str(end),group=group,eligible=eligible,max_kp=max_kp,
             min_dst_nT=minimum,days_kp_ge5=storm_days,max_daily_ap=max(float(r['AP_AVG']) for r in days),
             mean_f107_observed=float(np.mean([float(r['F10.7_OBS']) for r in days]))))
    chosen={g:next(r for r in candidates if r['group']==g and r['eligible']) for g in ('quiet','active','strong')}
    protocol=dict(selection_rule='2024 calendar-anchored starts Jan/Mar/May/Jul/Sep/Nov 1, exactly60d; earliest eligible in each weather class; no formation results read',
      class_definitions=dict(quiet='max Kp<5 AND min Dst>-50 nT',active='at least2 Kp>=5 days AND -250<min Dst<=-100 nT',strong='min Dst<=-250 nT; overrides other classes'),
      retreat_criteria=dict(A='km scale: return to control research',B='several100m to1km: requirements/control boundary',C='100 toseveral100m: applicability/environment/lifetime before optimization',D='strong/severe nearbelow100m: withdraw blanket control necessity claim'),
      phase1_decision='No new A(t) model without measured envelope change over4-5min; wave period/phase movement alone insufficient',
      raw_weather_sha256=hashlib.sha256(raw.read_bytes()).hexdigest(),omni_sha256=hashlib.sha256(omni.read_bytes()).hexdigest(),
      sources=['https://celestrak.org/SpaceData/SW-All.csv','https://spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni/omni2_2024.dat'],
      candidates=candidates,selected=chosen,missing_policy='reject; no imputation',fixed_before_new_formation_results=True,
      confounders='Different solar activity, date/season and Earth rotation phase; not a pure geomagnetic-only causal experiment')
    (OUT/'protocol.json').write_text(json.dumps(protocol,indent=2),encoding='utf-8')
    original=yaml.safe_load((ROOT/'configs/nrlmsis2_environment.yaml').read_text(encoding='utf-8'))
    configs=ROOT/'configs/real_environment_applicability';configs.mkdir(exist_ok=False)
    for group,choice in chosen.items():
        freeze(rows,np.datetime64(choice['start']),OUT/'weather'/group)
        spec=dict(original,epoch_utc=choice['start']+'T00:00:00',selection_rule=protocol['selection_rule'],epoch_reason=f'First eligible {group} weather class; preregistered protocol.json')
        (configs/(group+'.yaml')).write_text(yaml.safe_dump(spec,sort_keys=False,allow_unicode=True),encoding='utf-8')
    # The quiet driver series must exactly reproduce the saved January case.
    previous=ROOT/'outputs/nrlmsis2_20261010/weather_observed/drivers.csv'
    assert (OUT/'weather/quiet/drivers.csv').read_bytes()==previous.read_bytes()
    print(json.dumps(chosen,indent=2))


if __name__=='__main__': main()
