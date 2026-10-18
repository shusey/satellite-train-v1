"""Fixed two-level dwell intervention and descriptive deterministic contrasts."""
from .config import config_from_dict

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
