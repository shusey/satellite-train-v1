"""Post-processing and rigorously inactive-prefix reuse for NRLMSIS cases."""
import numpy as np


def pd_components(error, rate, kp, kd):
    """Exact finite-path controller B, batched over saved times."""
    def couple(edge):
        return np.concatenate([edge[...,:1],np.diff(edge,axis=-1),-edge[...,-1:]],axis=-1)
    return kp*couple(error), kd*couple(rate)


def inactive_prefix_proof(arrays, config, cache_seconds):
    """First potential HIGH command on all-LOW trajectory, at every control tick.

    Reuse is valid only for initially LOW machines, no transitions, and exactly
    coincident output/control ticks. Integral 2 s steps already end at these
    ticks, so enabling a controller which never commands cannot alter RK steps.
    End the reusable prefix at an earlier cache anchor to preserve cache phase.
    This proves equality for this numerical input, NOT stability or future
    inactivity. Beyond the returned prefix the real controller must be run.
    """
    t=np.asarray(arrays['time_s']);s=config.simulation;c=config.controller
    if (s.output_step_s!=c.control_period_s or
        not np.allclose(np.diff(t),c.control_period_s,rtol=0,atol=1e-9) or
        not np.isclose(t[0],0,rtol=0,atol=1e-9) or
        c.control_period_s/s.integration_step_s!=round(c.control_period_s/s.integration_step_s)):
        raise ValueError('Inactive reuse requires coincident ticks and integral RK step schedule')
    if np.any(arrays['mode']!=0) or np.any(arrays['area_m2']!=config.spacecraft.area_low_m2):
        raise ValueError('Inactive reuse requires all LOW throughout')
    # Preserve exact operation order in switching_function for independent check.
    edge=(config.derived.controller_kp_s2_inv*arrays['gap_error_m']+
          config.derived.controller_kd_s_inv*arrays['gap_rate_m_s'])
    q=np.concatenate([edge[:,:1],np.diff(edge,axis=1),-edge[:,-1:]],axis=1)
    np.testing.assert_array_equal(q,arrays['q_m_s2'])
    ix=np.flatnonzero(np.any(q>=c.threshold_on_m_s2,axis=1))
    crossing=float(t[ix[0]]) if ix.size else None
    if crossing is None: endpoint=float(t[-1])
    else: endpoint=max(0.,(np.floor(crossing/cache_seconds)-1)*cache_seconds)
    prefix=t<=endpoint
    if np.any(q[prefix]>=c.threshold_on_m_s2):raise ValueError('Prefix contains a possible command')
    return dict(first_potential_high_time_s=crossing,reusable_until_s=endpoint,
        reusable_control_ticks=int(prefix.sum()),maximum_q_in_reusable_prefix_m_s2=float(q[prefix].max()),
        threshold_on_m_s2=c.threshold_on_m_s2,
        minimum_margin_m_s2=float(c.threshold_on_m_s2-q[prefix].max()),
        all_ticks_checked=True, exactly_recomputed_q=True,
        justification='No possible LOW->HIGH command; same LOW areas, RHS, RK boundaries and density cache anchors')


def density_statistics(rho):
    rho=np.asarray(rho,float)
    if rho.ndim!=2 or np.any(rho<=0) or not np.all(np.isfinite(rho)):
        raise ValueError('Positive finite time-by-satellite density required')
    mean=rho.mean(axis=1);spread=np.ptp(rho,axis=1);adjacent=np.diff(rho,axis=1)
    return dict(min_kg_m3=float(rho.min()),max_kg_m3=float(rho.max()),
        time_mean_fleet_mean_kg_m3=float(mean.mean()),
        maximum_fleet_spread_kg_m3=float(spread.max()),
        maximum_fleet_spread_over_mean=float((spread/mean).max()),
        rms_differential_density_kg_m3=float(np.sqrt(np.mean((rho-mean[:,None])**2))),
        maximum_adjacent_absolute_difference_kg_m3=float(np.abs(adjacent).max()),
        rms_adjacent_difference_kg_m3=float(np.sqrt(np.mean(adjacent**2))),
        maximum_adjacent_difference_over_fleet_mean=float((np.abs(adjacent)/mean[:,None]).max()),
        definition='fleet spread over 21 satellites spanning 2000 km; adjacent differences across 100 km gaps')


def density_to_gap_acceleration(radius, rho, area, config):
    """Exact derivative of a0*diff(n(a)) under the existing averaged drag ODE.

    da/dt=-rho*Cd*A/m*sqrt(mu*a); dn/dt=-3*n*da/(2*a).
    This is averaged mean-motion gap acceleration, not Cartesian acceleration.
    """
    c=config.constants;s=config.spacecraft
    return (1.5*config.derived.reference_radius_m*s.drag_coefficient*c.mu_m3_s2/s.mass_kg
            *np.diff(np.asarray(rho)*np.asarray(area)/np.asarray(radius)**2,axis=-1))


def current_radius_arc_error(radius, phase, target_spacing_m):
    """Auxiliary near-circular arc proxy; never changes original e or controller.

    Adjacent mean radius times angular separation, minus a fixed physical target.
    This includes common orbit contraction. Mean orbital states do not support a
    claim of precise instantaneous real-space distance or geodesic separation.
    """
    a,u=np.asarray(radius,float),np.asarray(phase,float)
    if a.shape!=u.shape or a.shape[-1]<2:raise ValueError('Matching satellite states required')
    return .5*(a[...,:-1]+a[...,1:])*np.diff(u,axis=-1)-target_spacing_m
