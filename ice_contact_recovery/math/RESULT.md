# A finite-water reachability certificate for stream–cliff reconnection

1 November 2026. Standalone derivation and executable conditional forecast. No active Claude source or frozen campaign was changed.

## What this adds

The existing peak-only rescue inequality is incomplete for a finite flood. A peak can exceed the static activation discharge, while insufficient duration or water volume prevents the pond from reaching the abandoned toe. Here a sharp scalar reachability bound gives **both** minimum duration and minimum total water. It remains an optimistic necessary bound when the bypass incises but the toe does not lower.

This is proved mathematics for stated equations, not established historical originality. Scalar comparison, reservoir mass balance and optimal-control arguments are established tools. A potentially useful new application would be a quantitatively validated distinction between peak-feasible and actually reachable reattachment in an evolving ice/debris setting.

## Model and assumptions

Before first toe contact, let H be the pond free-surface elevation, z_t the fixed toe activation elevation, and A(H)>0 its area. Let R_0(H)>=0 be a static bypass discharge curve, nondecreasing and locally Lipschitz, and let I(t) be total inflow, bounded by 0<=I<=P. Water storage obeys

    A(H) H' = I(t) - R(t,H).

For the static model R=R_0. For the exclusion extension assume R(t,H)>=R_0(H) at every admitted H (for example, downward bypass incision with an otherwise unchanged monotone control). No toe discharge occurs before H=z_t. The initial state H_0<z_t and the complete initial pond storage are known. Volume means **total inflow**, not water above an ordinary-flow baseline. There is no uncounted precipitation, pore release, lateral inflow, meltwater, pond displacement or internal water source. Add all admitted water sources before using the certificate.

Any cliff-toe melting that lowers z_t during the pulse must be represented; a fixed-target bound cannot then be applied blindly. This model excludes unsteady inertial overshoot, waves, moving walls, pressure-driven flow and pressure-induced debris failure.

## Proposition and proof

Suppose P>R_0(z_t). Define

    G(h) = integral[H_0 to h] A(s)/(P-R_0(s)) ds,
    F(h) = P G(h).

Until first contact, P-R_0(H)>0. Along every admissible trajectory, even when R>=R_0,

    dG/dt = (I-R)/(P-R_0) <= 1,

and

    dF/dt = P(I-R)/(P-R_0) <= I.

The second inequality follows from

    I(P-R_0)-P(I-R) = P(R-R_0)+R_0(P-I) >= 0.

These inequalities also hold during draining intervals (H'<0); no assumption of monotone filling is needed. Integration up to the first-contact time gives

    t_contact >= t_min = G(z_t),
    V_contact >= V_min = P G(z_t).

For the static model, I=P until contact gives equality in both inequalities. Therefore, with a deadline T and total available water V, first contact is reachable **if and only if** T>=t_min and V>=V_min, when the forcing may be freely selected within the peak bound and stopped at contact. This is an idealized reachability test, not a recommendation to manipulate natural streams.

If P<=R_0(z_t), no finite-time upward crossing reaches z_t in the static model. At equality the target is approached asymptotically under the locally Lipschitz discharge law; above-static bypass discharge cannot improve reachability. A grain mobilization threshold or finite positive toe discharge requires a higher/modified target and its own discharge terms. Zero-depth first contact is not sustained cliff maintenance.

## Broad-crested control specialization

For R_0(H)=C[H-z_b]_+^p, p>0, and constant A,

    t_min = A integral[H_0 to z_t] dH/(P-C[H-z_b]_+^p),
    V_min = P t_min,
    P_critical = C[z_t-z_b]_+^p.

`recovery_budget.py` evaluates this one-dimensional quadrature; no CFD run is necessary. A measured/verified discharge curve can replace the ideal power law, but the software supplied here implements the power-law specialization only.

For g=z_t-z_b>0, smooth area A_t=A(z_t), and epsilon=P/P_critical-1 tending to zero from above,

    t_min = A_t/[p C g^(p-1)] log(1/epsilon) + O(1),
    V_min = A_t g/p log(1/epsilon) + O(1).

The coefficient is obtained by expanding the denominator near the target: P-R_0(H) = P_critical epsilon + p C g^(p-1)(z_t-H)+higher-order terms. Thus a peak-only classifier becomes arbitrarily overoptimistic close to its nominal threshold. This logarithm is a standard near-equilibrium slowing mechanism, not a new universal physical law.

## What would be a consequential environmental test

Freeze this bound using a reservoir storage curve and bypass discharge curve estimated independently from the CFD cases. Compare three predictors on held-out cases: (1) peak-only activation, (2) peak plus storage-only volume A(z_t-H_0), (3) the leakage-aware quadrature. Use rectangular pulses varying only duration and peak, then a separate pulse-order challenge. Do not fit C or A to each test outcome.

Score first-contact occurrence and time, with positive-contact and persistent ordinary-flow recovery treated as separate outcomes. Compare a frozen substrate with the same heat/flow forcing on the moving substrate. A larger unreachable domain caused by bypass incision, beyond storage leakage alone, would isolate the ice-specific contribution. Cases with toe lowering or additional meltwater must use a properly bounded extension or be declared outside scope.

No superiority over the existing simulator or real glacier models is established by the scalar examples. This result supplies a fast benchmark, an explicit baseline, and a theorem-backed rejection region, rather than a fitted prediction guarantee.

## Closest prior work checked

- [Kingslake et al. 2015, Modelling channelized surface drainage of supraglacial lakes](https://www.cambridge.org/core/journals/journal-of-glaciology/article/modelling-channelized-surface-drainage-of-supraglacial-lakes/E5154217B4D0307A242FE37C958F228E): primary full HTML inspected for lake-storage/incision competition. It already couples changing lake stage and channel incision; no broad first-coupling claim is tenable.
- [Schoof et al. 2023, The drainage of glacier and ice sheet surface lakes](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/drainage-of-glacier-and-ice-sheet-surface-lakes/6A699229F6F8C928BA52E83ABFEC16B7): primary full HTML discussion inspected. It already derives supply thresholds and different drainage regimes. The present target is reactivation of a previously abandoned competing toe under a finite peak and total-water budget, rather than lake drainage itself.
- The existing local `PULSE_PRIORITY_CHECK.md` and `NO_RESCUE_CERTIFICATE_AUDIT.md` already establish that timing, sediment-pulse hysteresis and thermal rescue barriers are not broadly new concepts. This package must be cited as a refinement of that work, not an independent priority claim.

Searches were targeted, not exhaustive. Historical novelty of this specific finite-budget application remains unestablished; the derivation alone does not justify a Nature-level claim.

## Run

    cd research/novelty_parallel_2026_11_01/math
    python3 -m unittest -v test_recovery_budget.py
    MPLCONFIGDIR=/private/tmp/recovery_budget_mpl python3 run_examples.py
    python3 recovery_budget.py --area 10 --coefficient .7 --bypass 0 --toe .2 --stage .05 --peak .1

The JSON and figure outputs are explicitly synthetic, generated from the declared scalar equations. They are not final experimental or environmental data.
