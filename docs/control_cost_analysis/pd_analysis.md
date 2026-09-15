# PD control analysis: implementation, surrogate, and limits

## 1. Exact graph correspondence

Let D be the (N-1)-by-N forward incidence matrix, with row i containing
-1 at satellite i and +1 at satellite i+1. Let
x_i=a0[lambda_i-lambda_b(t)]-(i-(N+1)/2)d for an arbitrary common reference
longitude lambda_b(t). Then e=Dx, edot=D xdot and L=D^T D is the open-path
Laplacian with end diagonal entries 1 and interior entries 2.

The implemented controller is exactly

    q = -D^T (Kp e + Kd edot) = -Kp L x - Kd L xdot.

This is unnormalized symmetric nearest-neighbor distributed PD. The end law
is degree one, not a missing-neighbor approximation or a periodic link.
The graph is connected, L 1=0, and sum_i q_i=0. These statements are algebraic
properties of the actual code, tested against switching_function.

Kp acts on relative position errors measured with a fixed reference radius;
Kd acts on the difference of orbital mean motions multiplied by that radius.
It is not Cartesian radial velocity, nor a finite difference of saved position.

The structure corresponds to second-order consensus / bidirectional platoon
feedback only after specifying how q would act on the plant. The existing
code applies area according to the sign/threshold of q; it never applies q
as an acceleration or uses a continuous inverse actuator.

## 2. Original plant and the sign of control

For an unforced atmosphere write f(a,A)=-CD rho(a) A sqrt(mu a)/m.
The original plant is adot=f, lambdadot=n(a). Direct differentiation gives

    a0 lambda_ddot = a0 n'(a) f(a,A)
                   = [3 a0 mu CD rho(a)/(2 m a²)] A = g(a) A > 0.

Increasing area decreases a and increases mean motion. Thus positive q,
which asks a satellite to catch up to an overly distant forward neighbor,
correctly selects more drag. This is reasonable orbital phasing feedback.
The virtual along-track acceleration has the opposite sign from the
instantaneous tangential drag force. Confusing these would reverse stability.

Common low drag is already present. Actual incremental area above LOW is
nonnegative and bounded. Since sum q=0, a nonzero signed q vector cannot
generally be realized as delta_A=q/g around the LOW boundary. A continuous
surrogate is therefore a diagnostic idealization, not a hidden implemented law.

## 3. Explicit assumptions for a continuous tangent model

1. No external density pulse; identical satellites and the existing exponential atmosphere.
2. Perturb around a common decaying near-circular trajectory a_b(t), initially a0,
   with common nominal LOW area A_L. Keep the reference-radius spacing definition.
3. Retain first-order perturbations in delta_a and longitude. Let v=xdot.
4. Replace switching by an ideal instantaneous, signed incremental area delta_A=q/g(t).
   This is physically infeasible at LOW for negative q. An interior nominal area
   could allow small signed perturbations, but is not simulated here or assigned
   a new numerical value.
5. For numerical eigenvalues and the following Lyapunov proof, freeze the
   time-dependent coefficients at t=0. Do not infer a general time-varying or
   hybrid theorem from frozen eigenvalues.

To derive the tangent before freezing, c(t)=a0 n'(a_b(t)), v=c delta_a,
delta_adot=f_a delta_a+f_A delta_A. Therefore

    vdot = [cdot/c + f_a] v + c f_A delta_A
         = b(t) v + g(t) delta_A,
    b(t) = f(a_b,A_L) [-1/H - 2/a_b] > 0,
    g(t) = 3 a0 mu CD rho(a_b)/(2 m a_b²) > 0.

The term cdot/c=-5 f/(2 a_b) must be included; simply setting xddot to q
discards both altitude-density feedback and the moving reference.
A finite-difference test of the original plant's along-track acceleration
verifies b and g independently.

With the ideal signed input and frozen b the surrogate is

    xdot = v,
    vdot = -Kp L x - (Kd L - b I) v.

## 4. Equilibria, eigenvalues, and a limited Lyapunov result

The path eigenvalues are ell_k=2[1-cos(k pi/N)], k=0,...,N-1.
Each disagreement mode k>=1 satisfies

    zddot_k + (Kd ell_k-b) zdot_k + Kp ell_k z_k = 0.

All disagreement roots have strictly negative real part if and only if

    Kp > 0 and Kd ell_1 > b.

For the present baseline:

- b=1.787876897303224e-8 s^-1
- g=1.7392686347694972e-4 m^-1 s^-2
- ell_1=0.022338347549742954
- slowest pair: -1.025243372446028e-6 +/- 6.84307733163415e-6 i s^-1
- slow-mode exponential amplitude time: about 11.29 days

The common mode has roots 0 and b: no absolute-position/altitude regulation.
Thus the full system is not asymptotically stable to a fixed absolute orbit.
In disagreement coordinates (project with P=I-11^T/N), the equilibrium is
Px=0, Pv=0: perfect spacing and common mean motion. Uniform x is free.

For constant coefficients, on the disagreement subspace,

    V = (1/2) v^T v + (Kp/2) x^T L x
    Vdot = -v^T (Kd L - b I) v <= -(Kd ell_1-b)||v||².

V is positive definite on that subspace. If Kd ell_1>b, the largest invariant
set with Vdot=0 is v=0,Lx=0, hence x=0 in disagreement coordinates.
LaSalle gives asymptotic convergence; the modal roots give exponential decay.
This is a proof ONLY for the frozen linear signed-input surrogate.
The numerical eigenvalues are saved for every mode in linear_modes.csv/json.

If b were additionally neglected, Kp=1/T², Kd=2/T gives modal damping ratio
sqrt(ell_k), not critical damping for every mode. Only ell_k=1 is critical.
A 6 h design_time is not the train settling time. Longer paths have
ell_1 ~ pi²/N² and weaker damping of long-wavelength modes. Finite-N stability
does not establish uniform string stability or bounded disturbance gain as N grows.

## 5. What the actual hybrid system does not inherit

- No fixed absolute equilibrium: adot<0 even in LOW, while density stays positive.
- The all-LOW, zero-error, equal-a manifold follows common orbital decay.
- With equal a after the pulse, nonzero constant spacing errors can persist
  whenever all q_i<h_on. This is a deadband set, not proof of convergence to e=0.
  The rule thresholds q, which mixes position and rate, not |e| directly.
- For the actual input, P g(a)(A-A_L) is neither q nor a linear feedback.
  It can excite the common drag cost even though sum q=0.
- Sampling, two positive hysteresis thresholds, finite slew, dwell, bounded
  unilateral actuation, pulse forcing and time-varying density invalidate a
  direct application of the continuous Vdot calculation.
- The nonlinearity is not a symmetric sign(q) relay. Negative q causes LOW,
  never thrust or negative drag; h_off is positive, not -h_on.
- Dwell prevents arbitrarily fast repeated switching per satellite but does
  not by itself guarantee formation stability, string stability or low cost.
- Position/rate tolerances describe one recovery test; they are not mission
  requirements or a theorem of permanent recovery.
- This 24 h deterministic experiment cannot prove lifetime benefit/cost,
  robustness to noise, repeated disturbances, or infinite-horizon propagation failure.

## 6. What should be tested or proved next

Define acceptable gap/rate error from mission requirements, then minimize
area exposure and attitude transitions subject to it. Analyze hybrid
invariant error sets and return maps with the existing slew/dwell pulse size.
Check disturbance-to-error and disturbance-to-area gains separately; add
bounded observation error only with justified settings. Examine gain and
dwell sensitivity before a large N/disturbance sweep. Any claimed string
stability requires an explicit norm, input location and uniform-N bound.

## Literature correspondence (not a claim that this controller was copied)

- [Olfati-Saber and Murray, consensus problems](https://murray.cds.caltech.edu/Consensus_problems_in_networks_of_agents_with_switching_topology_and_time-delays):
  graph Laplacian, disagreement functions and algebraic connectivity provide
  the structural comparison; that work's protocols do not prove this hybrid system.
- [Herman et al., bidirectional platoon scaling](https://arxiv.org/abs/1410.3943):
  nearest-neighbor coupling and graph-distance-dependent amplification motivate
  separating finite-system stability from string stability. Its assumptions
  and scaling theorem are not directly transferred to this open satellite path.
- [Leonard, Formationkeeping of Spacecraft via Differential Drag](https://ntrs.nasa.gov/citations/19870004033):
  differential drag as a formation actuator and transformed relative-motion
  models are relevant context, not the source of the implemented controller.
- [NASA EO-1 differential-drag demonstration](https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20170011151.pdf):
  links differential drag with changes in semi-major axis and along-track motion.
  Operational performance is not inferred from this simplified model.

