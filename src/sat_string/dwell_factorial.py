"""Fixed two-level dwell intervention and descriptive deterministic contrasts."""
from .config import config_from_dict
import numpy as np
from .dwell_intervention_metrics import distribution,followup_statistics

CASE_ORDER = ('H300-L300', 'H150-L300', 'H300-L150', 'H150-L150')


def factorial_config(base, high, low, *, days=30):
    if high not in (150, 300) or low not in (150, 300):
        raise ValueError('Only the predeclared 150/300 second levels are allowed')
    raw = base.to_dict()
    raw['controller']['dwell_high_time_s'] = float(high)
    raw['controller']['dwell_low_time_s'] = float(low)
    raw['simulation']['duration_s'] = days * 86400.
    result = config_from_dict(raw)
    before, after = base.to_dict(), result.to_dict()
    for item in (before, after):
        item['simulation'].pop('duration_s')
        for key in ('dwell_high_time_s', 'dwell_low_time_s'):
            item['controller'].pop(key, None)
    if before != after:
        raise ValueError('Unexpected non-dwell intervention')
    return result


def contrasts(values):
    """Input is in CASE_ORDER; no statistical uncertainty is implied."""
    a, b, c, d = values
    return dict(HIGH=(b+d-a-c)/2, LOW=(c+d-a-b)/2, interaction=d-b-c+a)


def extension_config(checkpoint, profile):
    """Accept a valid fixed-profile 30-day endpoint; change only duration to 60d."""
    source=config_from_dict(checkpoint['configuration'])
    if checkpoint['time_s'] != 30*86400 or source.simulation.duration_s != 30*86400:
        raise ValueError('Expected the completed 30-day checkpoint')
    if checkpoint['profiles'] != [profile] or checkpoint['disturbance_enabled']:
        raise ValueError('Profile or disturbance changed')
    high,low=source.controller.effective_dwell_high_s,source.controller.effective_dwell_low_s
    if (high,low) not in [(150,300),(300,150)]:
        raise ValueError('Extend only the two cross conditions')
    state=np.asarray(checkpoint['state'],dtype=float)
    if state.shape!=(2,source.formation.satellite_count) or not np.isfinite(state).all():
        raise ValueError('Invalid endpoint state')
    if np.any(state[0]<=source.constants.earth_radius_m) or np.any(np.diff(state[1])<=0):
        raise ValueError('Cannot extend a reversed or invalid formation')
    raw=source.to_dict();raw['simulation']['duration_s']=60*86400.
    return config_from_dict(raw)


def period_cycle_diagnostics(cycles,edges,start_s,end_s):
    """Summarize cycles assigned by LOW completion in an inclusive window."""
    selected=[r for r in cycles if start_s<=r['end_s']<=end_s]
    ids={(r['satellite_id'],r['cycle']) for r in selected}
    incident=[r for r in edges if (r['satellite_id'],r['cycle']) in ids]
    followed=[r for r in selected if r['next_HIGH_s'] is not None]
    qualified=[r for r in followed if abs(r['next_qD_neighbor'])>abs(r['next_qD_natural'])]
    return dict(start_s=start_s,end_s=end_s,completed_cycles=len(selected),
        pulse_velocity_contribution_m_s=distribution(r['pulse_g_impulse_m_s'] for r in selected),
        residual_abs_incident_gap_velocity_m_s=distribution(abs(r['velocity_after_m_s']) for r in incident),
        residual_max_incident_velocity_per_cycle_m_s=distribution(r['V_incident_after_m_s'] for r in selected),
        residual_worsening_fraction=float(np.mean([r['worsening_after'] for r in incident])) if incident else None,
        re_trigger=followup_statistics(selected,end_s),
        neighbor_exceeds_natural_fraction=len(qualified)/len(followed) if followed else None,
        neighbor_qualifying_intervals_per_day=len(qualified)/((end_s-start_s)/86400),
        neighbor_abs_D_change_median_m_s2=distribution(abs(r['next_qD_neighbor']) for r in followed)['median'])
