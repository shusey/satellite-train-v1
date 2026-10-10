#!/usr/bin/env python3
"""Freeze observed historical CelesTrak drivers; reject imputation/prediction."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from sat_string.io import prepare_output_directory,write_json


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--raw',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();out=prepare_output_directory(args.output)
    with args.raw.open(encoding='utf-8-sig') as f:
        reader=csv.DictReader(f);columns=reader.fieldnames;rows={r['DATE']:r for r in reader}
    epoch=np.datetime64('2024-01-01');end=epoch+np.timedelta64(60,'D')
    support=np.arange(epoch-np.timedelta64(40,'D'),end+np.timedelta64(41,'D'))
    subset=[]
    for date in support:
        row=rows[str(date)]
        if row['F10.7_DATA_TYPE']!='OBS':raise ValueError('Non-observed F10.7; do not silently select another period')
        flux=float(row['F10.7_OBS'])
        if not 0<flux<=400:raise ValueError('Flux outside pymsis valid range; no replacement allowed')
        for name in [f'AP{i}' for i in range(1,9)]+['AP_AVG','F10.7_OBS_CENTER81']:
            if not row[name] or not np.isfinite(float(row[name])) or float(row[name])<0:raise ValueError('Missing/invalid observed driver')
        subset.append(row)
    with (out/'observed_support.csv').open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=columns);w.writeheader();w.writerows(subset)
    def ap_at(time):
        day=str(time.astype('datetime64[D]'));slot=int((time-time.astype('datetime64[D]'))/np.timedelta64(3,'h'))+1
        return float(rows[day][f'AP{slot}'])
    drivers=[]
    for time in np.arange(epoch.astype('datetime64[h]'),end.astype('datetime64[h]')+np.timedelta64(1,'h'),np.timedelta64(3,'h')):
        day=time.astype('datetime64[D]');row=rows[str(day)]
        history=[float(row['AP_AVG'])]+[ap_at(time-np.timedelta64(h,'h')) for h in [0,3,6,9]]
        history += [float(np.mean([ap_at(time-np.timedelta64(h,'h')) for h in range(12,34,3)])),
                    float(np.mean([ap_at(time-np.timedelta64(h,'h')) for h in range(36,58,3)]))]
        drivers.append(dict(utc=str(time.astype('datetime64[s]')),
            f107_previous_day=float(rows[str(day-np.timedelta64(1,'D'))]['F10.7_OBS']),
            f107_center81=float(row['F10.7_OBS_CENTER81']),**{f'ap{i}':v for i,v in enumerate(history)}))
    with (out/'drivers.csv').open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(drivers[0]));w.writeheader();w.writerows(drivers)
    # Independently compare to the pinned wrapper's history implementation.
    import pymsis.utils as utils
    utils._F107_AP_PATH=out/'observed_support.csv';utils._DATA=None
    a,b,c=utils.get_f107_ap(np.asarray([r['utc'] for r in drivers],dtype='datetime64[s]'))
    ours=np.asarray([[r['f107_previous_day'],r['f107_center81']]+[r[f'ap{i}'] for i in range(7)] for r in drivers])
    if not np.array_equal(ours,np.column_stack([a,b,c])):raise ValueError('Driver construction differs from pymsis')
    write_json(out/'provenance.json',dict(source_url='https://celestrak.org/SpaceData/SW-All.csv',
        source_snapshot_sha256=hashlib.sha256(args.raw.read_bytes()).hexdigest(),
        support_start=str(support[0]),support_end=str(support[-1]),epoch_utc=str(epoch),end_utc=str(end),
        rows=len(drivers),missing_value_policy='reject; no local filling',F107_qualifiers=['OBS'],
        pymsis_history_exact_match=True,selection_changed_after_results=False,
        source_description='CelesTrak consolidated GFZ geomagnetic / Canadian observed solar flux; see documentation',
        geomagnetic_activity=-1,inputs='previous-day observed F107; current centered81; seven Ap entries'))
    print('Frozen historical drivers:',len(drivers),'epoch',epoch,'end',end,flush=True)


if __name__=='__main__':main()
