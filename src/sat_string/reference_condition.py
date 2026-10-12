"""Frozen sensitivity profiles, not estimates of a spacecraft population."""
import json
from pathlib import Path
import numpy as np

STATEMENT = ('0.5 m・0.5%は実機母集団から推定された統計的nominal値ではなく、'
             '文献で確認された現実的な桁を踏まえた工学的Reference sensitivity conditionである')


def normalize_saved(values, target, ddof=1):
    """Preserve the saved spatial pattern; center and set sample std (ddof=1)."""
    z = np.asarray(values, dtype=float)
    if z.shape != (21,) or not np.all(np.isfinite(z)) or target < 0:
        raise ValueError('Expected 21 finite values and nonnegative target')
    z = z-z.mean()
    if z.std(ddof=ddof) == 0:
        if target: raise ValueError('Cannot normalize a constant profile')
        return z
    return z*(target/z.std(ddof=ddof))


def make_reference(source):
    old = source['profiles']
    pa = next(p for p in old if p['case']=='A')
    pb = next(p for p in old if p['case']=='B')
    if pa['id']!='a_0.1_random_4' or pb['id']!='B_0.003_random_2' or pb['B_seed']!=10002:
        raise ValueError('Unexpected saved source profile')
    da = normalize_saved(pa['delta_a_m'], .5)
    eps = normalize_saved(pb['epsilon'], .005)
    zero = np.zeros(21)
    profiles = []
    for case, a, b in [('A',da,zero),('B',zero,eps),('AB',da,eps)]:
        profiles.append(dict(id=case+'_ref', case=case, control_enabled=True,
            delta_a_m=a.tolist(), epsilon=b.tolist(), B_multiplier=(1+b).tolist(),
            realized_sigma_a_m=float(a.std(ddof=1)),
            realized_sigma_B_fraction=float(b.std(ddof=1)), ddof=1))
    profiles += [dict(p, case=p['case']+'_OFF', control_enabled=False) for p in profiles[:3]]
    return dict(statement=STATEMENT, std_definition='sample standard deviation; ddof=1',
        reference_delta_a_profile=da.tolist(), reference_B_profile=eps.tolist(),
        source_A_id=pa['id'], source_B_id=pb['id'], A_rng_seed=4, B_rng_seed=10002,
        source_A_sample_std=float(np.std(pa['delta_a_m'],ddof=1)),
        source_B_sample_std=float(np.std(pb['epsilon'],ddof=1)), profiles=profiles)


def load_frozen(path):
    obj = json.loads(Path(path).read_text(encoding='utf-8'))
    ps = obj['profiles']
    if len(ps)!=6: raise ValueError('Expected three ON and three OFF profiles')
    da=np.asarray(obj['reference_delta_a_profile']); eps=np.asarray(obj['reference_B_profile'])
    for x,target in [(da,.5),(eps,.005)]:
        if not np.isclose(x.mean(),0,atol=1e-14) or not np.isclose(x.std(ddof=1),target,rtol=1e-14,atol=0):
            raise ValueError('Frozen mean/std changed')
    for i,p in enumerate(ps):
        a = da if i%3 in (0,2) else np.zeros(21)
        b = eps if i%3 in (1,2) else np.zeros(21)
        if p['delta_a_m']!=a.tolist() or p['epsilon']!=b.tolist() or p['B_multiplier']!=(1+b).tolist():
            raise ValueError('A/B/AB profile mismatch')
        if p['control_enabled'] != (i<3): raise ValueError('Control flag mismatch')
    return obj
