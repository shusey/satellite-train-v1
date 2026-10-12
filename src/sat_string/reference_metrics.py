"""Reference-specific additions to existing transfer metrics; no force changes."""
import numpy as np
from .dispersion_transfer_metrics import evaluate_transfer
from .necessity import improvement, exceedance_metrics
from .residual import cycle_pairs


def evaluate_reference(arrays, config, events, low_arrays):
    result=evaluate_transfer(arrays,config,events,low_arrays)
    t=arrays['time_s']; e=arrays['gap_error_m']; le=low_arrays['gap_error_m']
    bad=np.flatnonzero(np.any(low_arrays['gap_distance_m']<=0,axis=1))
    cut=int(bad[0]) if len(bad) else len(t)
    if cut>=2:
        result['matched_pre_allLOW_reversal_RMS_improvement']=improvement(
            float(np.sqrt(np.mean(e[:cut]**2))),float(np.sqrt(np.mean(le[:cut]**2))))
        result['matched_pre_allLOW_reversal_thresholds']=[dict(threshold_m=x,
            ON=exceedance_metrics(t[:cut],np.abs(e[:cut]).max(axis=1),x),
            OFF=exceedance_metrics(t[:cut],np.abs(le[:cut]).max(axis=1),x)) for x in (100,500,1000)]
    result['evaluation_thresholds_note']='100/500/1000 m are strict/reference/loose evaluation thresholds, not mission requirements'
    result['first_reversal_resolution_s']=float(np.max(np.diff(t)))
    pairs=cycle_pairs(events,float(t[-1]))
    completions=[p['return_event']['complete_time'] for p in pairs
                 if p['return_event'] and p['return_event']['complete_time']<=t[-1]]
    for segment in result['segments']:
        lo=segment['start_day']*86400;hi=segment['end_day']*86400
        last=hi==t[-1]
        count=sum(lo<=x and (x<=hi if last else x<hi) for x in completions)
        segment['completed_cycles']=count
        segment['cycles_per_day']=count/((hi-lo)/86400)
        segment['high_equivalent_satellite_hours_per_day']=segment['high_equivalent_satellite_s_per_day']/3600
    result['cycle_interval_counting']='Completed L-H-L counted at LOW completion; half-open intervals except final endpoint'
    return result
