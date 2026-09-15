# Repository understanding (before analysis changes)

Source: default branch main, commit 9449a78 (num1). All application modules,
four scripts, tests, current config and saved run configurations were inspected.
No notebooks or other application configuration files were found.
The tracked .venv is third-party runtime material, not the research model.

## Representative experiment and provenance

Use configs/baseline.yaml without physical modifications. It matches the saved
runs/experiment_13_1.5/normalized_config.yaml. The older runs/standard_suite is
72 h, dt=1 s, thresholds 1e-7/5e-8; it is NOT the current baseline.
Other saved runs include N=51 / 140 h and alternative hysteresis values.
Their NPZ/CSV/PNG artifacts were not committed, so only saved JSON/logs can be
compared until histories are regenerated. Never overwrite those runs.

| Parameter | Implemented baseline | Source |
|---|---:|---|
| N / graph | 21 / finite undirected open path | configs/baseline.yaml |
| target spacing | 100000 m, measured using fixed reference radius | controller.gap_kinematics |
| reference altitude / radius | 400000 / 6778137 m | config.derived |
| mass / CD | 4 kg / 2.2 | spacecraft |
| low / high area | 0.01 / 0.03 m² | spacecraft |
| density at reference altitude | 3.584968e-12 kg/m³ | atmosphere |
| scale height | 58310 m | atmosphere |
| design time | 21600 s | controller |
| Kp / Kd | 2.143347050754458e-9 s⁻² / 9.259259259259259e-5 s⁻¹ | config_from_dict |
| on / off | 1.5e-7 / 1.4e-7 m/s² (both positive) | controller |
| control / command delay | 10 / 0 s | controller |
| slew / dwell after completion | 120 / 300 s | attitude |
| duration / RK4 max step / output | 86400 / 2 / 10 s | simulation |
| recovery tolerances / hold | 50 m / 0.001 m/s / 1800 s | recovery |

Initial a_i=a0, lambda_i=(i-(N+1)/2)d/a0 (1-based); all attitudes LOW.
No pinned leader, no periodic wrapping, no relabeling on order reversal.
End satellites use their sole neighbor. No direct absolute-position feedback.

## Exact implemented equations

Define n(a)=sqrt(mu/a³), a0=R_E+h0, d=target spacing.
State: a_i (near-circular averaged semi-major axis) and unwrapped mean longitude lambda_i.

    da_i/dt = -rho_i CD A_i/m sqrt(mu a_i)
    d lambda_i/dt = n(a_i)
    h_i = a_i - R_E
    e_i = a0 (lambda_(i+1)-lambda_i) - d,  i=1,...,N-1
    de_i/dt = a0 [n(a_(i+1))-n(a_i)]
    p_i = Kp e_i + Kd de_i/dt
    q_1 = p_1
    q_i = p_i - p_(i-1), i=2,...,N-1
    q_N = -p_(N-1)
    Kp = 1/design_time², Kd = 2/design_time

q is a switching signal, NOT an applied acceleration. Every control tick selects
all commands from the same pre-update snapshot:

    LOW and q_i >= h_on -> HIGH
    HIGH and q_i <= h_off -> LOW
    otherwise hold

Commands also require stable mode, no pending command, and expired dwell time.
Modes are 0=LOW, 1=LOW_TO_HIGH, 2=HIGH, 3=HIGH_TO_LOW.
During slew, S(z)=3z²-2z³, z clipped to [0,1], interpolates actual area.
Dwell starts at slew completion, not command issue. With baseline settings,
a LOW->HIGH->LOW episode has at least 120+300+120 s of elevated area.
The stable HIGH portion is at least 300 s; its area-equivalent exposure is at
least 420 s. This minimum pulse is consequential even if q briefly exceeds h_on.

    rho_i = rho0 exp[-(a_i-a0)/H] (1+0.5 G_s(s_i,t) G_t(t))
    s_i = a0 [lambda_i - n(a0)t]
    s_c(t) = 0
    G_s = exp[-(s_i-s_c)²/(2 sigma²)]
    sigma = 100000 / (2 sqrt(2 ln 2)) m
    G_t = [1-cos(2 pi (t-10800)/3600)]/2 within [10800,14400], else 0

The external pulse is spatially Gaussian (noncompact support), so the nearest
neighbors receive direct forcing too. A observed first-switch wave is not by
itself proof of exclusively neighbor-mediated causality. center_of_train sets
s_c=0 in the fixed reference rotating frame; it does not follow the decaying
train's instantaneous center.

## Integrator / existing evaluation

Classical RK4 reevaluates density and interpolated area at all four stages.
Steps stop at control/output/delayed-start/completion events.
At each event: pending starts, completions, simultaneous controller, output.
For this baseline pulse endpoints are aligned to the step grid.
Stored histories: time, a, lambda, h, actual area, mode, density, q,
gap error/rate/distance; switch records include command/start/completion times.

Existing metrics include per-gap/global peak error, minimum separation,
order reversal, peak-error amplification, peak-error speed regressions,
recovery, switching count, sampled HIGH duration, excess area integral,
relative altitude spread and altitude difference against the all-low case.
Existing switching_wave.py gives first LOW->HIGH starts and side regressions
excluding the nearest three satellites. It lacks aggregate cost analysis,
range/censoring interpretation, energy, and cross-case trade-off reporting.
The new analysis will add these without altering simulation dynamics.

## Planned bounded analysis

1. Execute the original three baseline cases before edits to the core.
2. Compare with saved experiment_13_1.5 metrics and event log.
3. Add offline event-exact HIGH duration and smooth-slew equivalent exposure,
   signed changes in a/h/specific and total orbital energy, and cost decomposition.
4. Use seven h_on values already in the existing threshold grid:
   1.0, 1.3, 1.38, 1.5, 1.6, 1.8, 2.0 (times 1e-7).
   Preserve existing width h_on-h_off=1e-8; include the baseline unchanged.
   These are experiment design choices, not new physical assumptions.
5. Verify dt=2 vs dt=1 s for baseline and any selected transition case.
   Compare unchanged-core histories before/after analysis additions.
6. Keep mathematical continuous-input surrogate guarantees separate from
   unproven binary, unilateral, sampled, slew-limited hybrid dynamics.

