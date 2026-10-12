"""Saved OFF trajectories only: compare dispersion increments above zero baseline."""
import json
import numpy as np
import yaml
from run_reference_dispersion import ROOT,SPEC,DATA
from sat_string.io import write_json


def main():
    spec=yaml.safe_load(SPEC.read_text());out=ROOT/spec['output']
    profile=json.loads((DATA/'profiles.json').read_text())
    with np.load(ROOT/spec['baseline_run']/'results.npz') as file:
        t0=file['time_s'];e0=file['gap_error_m']
    results=[]
    for c,index,oldid,scale in [('A',3,profile['source_A_id'],.5/profile['source_A_sample_std']),
                                ('B',4,profile['source_B_id'],.005/profile['source_B_sample_std'])]:
        # Early response: both old and new OFF retain their original ordering.
        path=ROOT/'outputs/initial_dispersion_20261011/initial_7d'
        meta=json.loads((path/'metadata.json').read_text());j=next(i for i,p in enumerate(meta['profiles']) if p['id']==oldid)
        r0=meta['configuration']['derived']['reference_radius_m']
        target=meta['configuration']['derived']['target_spacing_rad']
        t=np.load(path/'time_s.npy');u=np.load(path/'lambda_rad.npy',mmap_mode='r')[:,j]
        old=r0*(np.diff(u,axis=1)-target)
        rt=np.load(out/'batch_30d/time_s.npy');ru=np.load(out/'batch_30d/lambda_rad.npy',mmap_mode='r')
        ix=np.searchsorted(rt,t);zero_ix=np.searchsorted(t0,t)
        if not np.array_equal(rt[ix],t) or not np.array_equal(t0[zero_ix],t):raise ValueError('Time grids differ')
        ref=r0*(np.diff(ru[ix,index],axis=1)-target)
        zero=e0[zero_ix]
        for days in (1,7):
            mask=t<=days*86400;effect=ref[mask]-zero[mask];prediction=scale*(old[mask]-zero[mask])
            difference=effect-prediction
            results.append(dict(case=c,days=days,actual_saved_profile_scale=scale,
                old_max_error_m=float(np.abs(old[mask]).max()),reference_max_error_m=float(np.abs(ref[mask]).max()),
                max_incremental_scaling_difference_m=float(np.abs(difference).max()),
                RMS_incremental_scaling_difference_m=float(np.sqrt(np.mean(difference**2))),
                RMS_relative_incremental_scaling_difference=float(np.sqrt(np.mean(difference**2)/np.mean(effect**2)))))
    write_json(DATA/'early_OFF_proportionality.json',dict(results=results,
        definition='Compare e_ref-e_zero with profile_scale*(e_old-e_zero); no new physical run',
        interpretation='early fixed-label windows before ordering reversal; not a long-term linearity guarantee'))
    print(results)


if __name__=='__main__':main()
