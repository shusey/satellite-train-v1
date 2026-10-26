"""Finite open-path N scaling under one explicit generalized node pulse.

No geographic/density pulse, random initial offsets or B mismatch are used.
The original density law, actuator, controller and RK4 remain in use.
"""
from dataclasses import dataclass, asdict
import numpy as np
from .config import config_from_dict


@dataclass(frozen=True)
class CentralNodePulse:
    count: int
    amplitude_m_s2: float = 1e-6
    start_time_s: float = 10800.
    duration_s: float = 3600.

    def __post_init__(self):
        if self.count < 3 or self.count % 2 != 1:
            raise ValueError('Odd N>=3 required')
        if not np.all(np.isfinite([self.amplitude_m_s2, self.start_time_s, self.duration_s])):
            raise ValueError('Finite pulse parameters required')
        if self.duration_s <= 0 or self.start_time_s < 0 or self.amplitude_m_s2 == 0:
            raise ValueError('Nonzero amplitude, positive duration, nonnegative start required')

    def __call__(self, t):
        """Flat nominal/disturbed batch: only member 1 has central node input."""
        values = np.zeros(2 * self.count)
        if self.start_time_s < t < self.start_time_s + self.duration_s:
            phase = 2*np.pi*(t-self.start_time_s)/self.duration_s
            values[self.count + self.count//2] = self.amplitude_m_s2 * .5*(1-np.cos(phase))
        return values

    def metadata(self):
        c = self.count//2
        return dict(**asdict(self), units='m/s^2', waveform='raised_cosine',
            definition='d_j = d/dt[a_ref*n(a_j)]_external; w_i=d_(i+1)-d_i',
            members=['nominal', 'disturbed'], central_satellite_id=c+1,
            positive_edge_id=c, negative_edge_id=c+1,
            boundaries_s=[self.start_time_s, self.start_time_s+self.duration_s])


def scaling_configuration(base, count, days):
    """Only N, duration and disabling the legacy density pulse differ from base."""
    raw = base.to_dict()
    raw['formation']['satellite_count'] = count
    raw['simulation']['duration_s'] = days*86400.
    raw['disturbance']['enabled'] = False
    return config_from_dict(raw)


def gap_graph_coordinates(count):
    """Incident edges are r=0 on either side of source node, IDs are 1-based."""
    c = count//2
    edges = np.arange(count-1)
    left = edges < c
    distance = np.where(left, c-1-edges, edges-c)
    return distance, np.where(left, 'left', 'right')


def response_summary(time_s, response, *, amplitude, terminal_window_s,
                     pulse_end_s, recovery_band_m=5., hold_s=86400.):
    """Full saved 10-s response, not downsampled plot data. Final peak is censored."""
    t = np.asarray(time_s); x = np.asarray(response)
    if x.ndim != 2 or len(t) != len(x) or len(t)<2 or np.any(np.diff(t)<=0):
        raise ValueError('Ordered time and matching time-by-gap response required')
    if amplitude == 0 or not np.all(np.isfinite(x)):
        raise ValueError('Finite response and nonzero amplitude required')
    absolute = np.abs(x)
    peak_indices = np.argmax(absolute, axis=0)
    peaks = absolute[peak_indices, np.arange(x.shape[1])]
    envelope = np.max(absolute, axis=1)
    maxgap = int(np.argmax(peaks))
    terminal = t >= t[-1]-terminal_window_s
    after = t >= pulse_end_s
    bad = np.flatnonzero(after & (envelope > recovery_band_m))
    start = int(np.searchsorted(t, pulse_end_s)) if len(bad)==0 else int(bad[-1]+1)
    recovery = float(t[start]) if start<len(t) and t[-1]-t[start]>=hold_s else None
    return dict(P_max_m=float(peaks[maxgap]), gain_s2=float(peaks[maxgap]/abs(amplitude)),
        max_response_gap_id=maxgap+1, peak_time_s=float(t[peak_indices[maxgap]]),
        peak_at_observation_end=bool(peak_indices[maxgap]==len(t)-1),
        E_ult_m=float(np.max(envelope[terminal])),
        final_max_response_m=float(envelope[-1]),
        rms_response_m=float(np.sqrt(np.mean(x*x))),
        recovery_band_m=recovery_band_m, recovery_hold_s=hold_s,
        recovery_candidate_s=recovery,
        verification_horizon_s=None if recovery is None else float(t[-1]-recovery),
        peak_by_gap_m=peaks.tolist(), peak_time_by_gap_s=t[peak_indices].tolist())
