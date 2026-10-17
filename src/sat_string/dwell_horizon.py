"""Join saved continuation parts without duplicating checkpoint observations."""
import numpy as np


def joined_saved_parts(parts,member,keys):
    """Exact state/actuator boundary; density permits existing cache reanchoring."""
    paths=list(parts);times=[np.load(p/'time_s.npy',mmap_mode='r') for p in paths]
    for previous,current in zip(times,times[1:]):
        if previous[-1]!=current[0]:raise ValueError('Noncontiguous continuation times')
    t=np.concatenate([x if i==0 else x[1:] for i,x in enumerate(times)])
    if np.any(np.diff(t)<=0):raise ValueError('Duplicate or reversed observations')
    arrays={};checks=[]
    for key in keys:
        pieces=[np.load(p/f'{key}.npy',mmap_mode='r')[:,member] for p in paths]
        for i,(previous,current) in enumerate(zip(pieces,pieces[1:])):
            left=np.asarray(previous[-1]);right=np.asarray(current[0])
            density=key=='density_kg_m3'
            same=np.allclose(left,right,rtol=1e-8,atol=0) if density else np.array_equal(left,right)
            if not same:raise ValueError('Checkpoint state changed: '+key)
            checks.append(dict(key=key,time_s=float(times[i][-1]),
                max_abs_difference=float(np.max(abs(left.astype(float)-right.astype(float)))),
                max_relative_difference=float(np.max(abs(left-right)/abs(left))) if density else None,
                tolerance='density rtol=1e-8' if density else 'exact'))
        arrays[key]=np.concatenate([x if i==0 else x[1:] for i,x in enumerate(pieces)])
    return t,arrays,checks
