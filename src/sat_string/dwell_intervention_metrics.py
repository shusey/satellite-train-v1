"""Small offline statistics for a finite-horizon dwell intervention."""
import numpy as np


def distribution(values):
    x=np.asarray(list(values),dtype=float)
    x=x[np.isfinite(x)]
    if not len(x):return dict(n=0,min=None,p25=None,median=None,p75=None,p90=None,max=None)
    qs=np.quantile(x,[0,.25,.5,.75,.9,1])
    return dict(n=int(len(x)),**dict(zip(['min','p25','median','p75','p90','max'],map(float,qs))))


def followup_statistics(rows, end_s):
    """Do not count an unobserved re-trigger as infinite wait or as a failure.

    A threshold outcome is known if a re-trigger is observed within that
    threshold or if follow-up reaches the threshold. Also return the simpler
    observed-waits-only proportion, labelled as such.
    """
    waits=[r['next_wait_s'] for r in rows if r['next_wait_s'] is not None]
    output=dict(observed_wait_s=distribution(waits),completed_cycles=len(rows),
                right_censored=sum(r['next_wait_s'] is None for r in rows))
    for threshold in (1800.,3600.):
        yes=sum(r['next_wait_s'] is not None and r['next_wait_s']<=threshold for r in rows)
        known=sum((r['next_wait_s'] is not None and r['next_wait_s']<=threshold)
                  or end_s-r['end_s']>=threshold for r in rows)
        output[f'within_{int(threshold)}s']=dict(events=yes,known_outcomes=known,
            fraction=yes/known if known else None,unresolved=len(rows)-known,
            observed_waits_only_fraction=yes/len(waits) if waits else None)
    return output
