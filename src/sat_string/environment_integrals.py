"""Passive environmental diagnostics. No controller or dynamics modifications."""
import numpy as np


def all_low_integrals(radius_m, initial_radius_m, mu, cd, area, mass):
    """ODE identities for the existing circular averaged ALL-LOW model.

    sqrt(a)=sqrt(a0)-Cd*A*sqrt(mu)/(2*m)*integral(rho dt).
    J_drag=integral(-rho*Cd*A*mu/(2*m*a) dt)
          =sqrt(mu)*(1/sqrt(a0)-1/sqrt(a)).
    Reconstruct effective integrals consistent with saved physical trajectory;
    not independent observational evidence. Discontinuity/RK errors remain.
    Stable algebra avoids subtracting nearly equal square roots/reciprocals.
    Units: J_rho kg*s/m^3; signed J_drag m/s (negative for orbit decay).
    """
    a=np.asarray(radius_m,float);a0=np.asarray(initial_radius_m,float)
    if np.any(a<=0) or np.any(a0<=0) or not all(np.isfinite(x) and x>0 for x in (mu,cd,area,mass)):
        raise ValueError('Positive radii and physical constants required')
    if not np.all(np.isfinite(a+a0)):raise ValueError('Finite radii required')
    root,root0=np.sqrt(a),np.sqrt(a0)
    root_difference=(a0-a)/(root0+root)
    jr=root_difference/(cd*area*np.sqrt(mu)/(2*mass))
    jd=-np.sqrt(mu)*root_difference/(root0*root)
    return jr,jd


def supported_integral(time_s, values, valid, max_gap_s, uncertainty=None):
    """Trapezoids ONLY between consecutive valid samples within max_gap.

    Do not bridge rejected samples or gaps. Sum is over observed support only,
    not an estimate of the unobserved full interval. Correlated uncertainty
    envelope integrates absolute supplied per-sample error, NOT a confidence
    interval or assumed independent-error reduction.
    absolute_interval_area sums absolute signed trapezoid areas; it is not
    integral(abs(y)) when y crosses zero within a sample interval.
    """
    t=np.asarray(time_s,float);y=np.asarray(values,float);good=np.asarray(valid,bool)
    if (t.ndim!=1 or len(t)<2 or y.shape!=t.shape or good.shape!=t.shape
            or not np.all(np.isfinite(t)) or np.any(np.diff(t)<=0)
            or not np.isfinite(max_gap_s) or max_gap_s<=0):
        raise ValueError('Matching increasing series required')
    dt=np.diff(t);keep=good[:-1]&good[1:]&np.isfinite(y[:-1]+y[1:])&(dt<=max_gap_s)
    pieces=np.where(keep,.5*(y[:-1]+y[1:])*dt,0.)
    edges=np.diff(np.r_[False,keep,False].astype(int))
    starts,stops=np.flatnonzero(edges==1),np.flatnonzero(edges==-1)
    segments=[dict(start_s=float(t[s]),end_s=float(t[e]),duration_s=float(t[e]-t[s]),integral=float(pieces[s:e].sum())) for s,e in zip(starts,stops)]
    envelope=None
    if uncertainty is not None:
        u=np.asarray(uncertainty,float)
        if u.shape!=t.shape or np.any(u[good]<0) or not np.all(np.isfinite(u[good])):raise ValueError('Valid nonnegative uncertainties required')
        envelope=float(np.sum(.5*(u[:-1][keep]+u[1:][keep])*dt[keep]))
    return dict(integral_observed_support=float(pieces.sum()),absolute_interval_area=float(np.abs(pieces).sum()),
      coverage_s=float(dt[keep].sum()),coverage_fraction=float(dt[keep].sum()/(t[-1]-t[0])),
      max_contiguous_coverage_s=max((s['duration_s'] for s in segments),default=0.),
      fully_correlated_error_envelope=envelope,segments=segments)


def gap_rate_from_exposure(j_rho, initial_radius_m, reference_radius_m, mu, cd, area, mass):
    """Map the model's density impulse into mean-motion gap rate, not delta-v."""
    root=np.sqrt(initial_radius_m)-cd*area*np.sqrt(mu)/(2*mass)*np.asarray(j_rho)
    if np.any(root<=0):raise ValueError('Invalid reconstructed orbit')
    return reference_radius_m*np.diff(np.sqrt(mu)/root**3,axis=-1)
