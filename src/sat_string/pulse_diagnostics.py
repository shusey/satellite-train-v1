"""Offline force attribution on an observed trajectory, never a counterfactual.

For g_j=r_ref*n(a_j), v_i=g_(i+1)-g_i and
gdot_j=1.5*r_ref*mu*rho_j*(Cd/m)*A_j*(1+eps_j)/a_j**2.
Splitting this into common LOW, fixed B mismatch, and excess area is exact
on the saved trajectory. Integrals use saved-sample trapezoidal quadrature.
Changing an actual command would also change future states and density;
these additive terms do not simulate that intervention.
"""
import numpy as np


def node_map(edges):
    """Open-path incidence used by q: [e0, e1-e0, ..., -e_last]."""
    x = np.asarray(edges)
    return np.concatenate((x[..., :1], np.diff(x, axis=-1), -x[..., -1:]), axis=-1)


def acceleration_terms(radius, density, area, epsilon, *, r_ref, mu, cd, mass, area_low):
    a, rho, ar, eps = map(np.asarray, (radius, density, area, epsilon))
    if np.any(a <= 0) or np.any(rho < 0) or np.any(1+eps <= 0):
        raise ValueError('Invalid radius, density, or drag multiplier')
    common = 1.5*r_ref*mu*rho*(cd/mass)/a**2
    return dict(low=common*area_low, mismatch=common*area_low*eps,
                control=common*(ar-area_low)*(1+eps))


def cumulative_trapezoid(values, time):
    y, t = np.asarray(values), np.asarray(time)
    if len(y) != len(t) or np.any(np.diff(t) <= 0):
        raise ValueError('Matching values and strictly increasing time required')
    dt = np.diff(t).reshape((-1,)+(1,)*(y.ndim-1))
    return np.concatenate((np.zeros_like(y[:1]), np.cumsum((y[:-1]+y[1:])*dt/2, axis=0)))


def local_control_split(node_impulse, satellite):
    """Own and other satellite contributions to every gap during a local pulse."""
    x = np.asarray(node_impulse)
    own = np.zeros_like(x)
    own[satellite] = x[satellite]
    return np.diff(own), np.diff(x-own)


def reversal_flags(before, after, own_impulse):
    """Kinematic diagnostics; zero velocity is NOT the PD controller's target."""
    b, a, u = map(np.asarray, (before, after, own_impulse))
    return dict(actual_sign_reversal=(b*a < 0), actual_abs_growth=(abs(a) > abs(b)),
                own_opposes_initial=(b*u < 0),
                own_alone_crosses_zero=(b*(b+u) < 0))
