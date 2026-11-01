# Bare-corridor ordering theorem and a topology counterexample

31 October 2026. Conditional mathematics for an ideal parallel-route / plug-flow heat model. This does not assert an unobserved glacier discovery or historical priority.

## 1. Model and exact hypotheses

Three parallel routes connect the same upstream reservoir and downstream datum: a bare toe corridor of width W_a=(1-f)W>0; a covered toe lane of width W_c=fW (absent when f=0); and a bare bypass of width W_b>0. Each has a flat, evolving ice floor z_i(t), length L, and ideal rectangular control law

    Q_i = K_0 W_i [eta-z_i]_+^(3/2),
    sum_i Q_i=I(t)>0.

K_0 is identical per unit width. The covered lane's hydraulic floor includes its nonnegative debris thickness, not only the buried ice floor. If thickness evolves, its geometric evolution must be included explicitly below.

The bare routes receive identical inlet temperature T(t)>=0, identical transfer coefficient h(t)>=0, the same length L, and constant positive ice/water properties. At fixed Q_i, the declared plug-flow law removes heat

    P_i=rho_w c_p Q_i T [1-exp(-h W_i L/(rho_w c_p Q_i))].

The limiting value is P_i=0 at Q_i=0. Its planform-average ice-floor lowering rate is

    z_i'=-P_i/(rho_i L_f W_i L).

This is the fixed-throughflow plug-flow approximation: meltwater mass addition to the along-channel discharge is neglected. That approximation is small when the water Stefan number is small; it must be declared rather than called exact coupled mass conservation. The comparison theorem is exact for these stated equations. It is not yet proved for every full free-surface Stefan system.

The covered route may have h_c(t)<=h(t), but must have the same nonnegative inlet T, same L, and the same functional heat law. Its cover thickness is fixed or nondecreasing for the strengthened wetness result. No bank collapse, lane closure, lateral thermal conduction, sidewall melting, sediment relocation, or differential ice motion changes route elevations outside these equations.

## 2. Monotonicity lemma

Write q=Q_i/W_i, a=hL/(rho_w c_p), and

    m(q;h,T)=rho_w c_p T/(rho_i L_f L)
             q [1-exp(-a/q)].

For q>0,

    d/dq {q[1-exp(-a/q)]}
       =1-(1+a/q)exp(-a/q)>=0,

strictly positive when a>0. This follows from exp(x)>1+x for x>0. The continuous limit at q=0 is zero. The function is also nondecreasing in h and T. Therefore at equal h,T,L, a larger per-width discharge gives at least as much floor lowering.

An equal-per-width inflow alone is insufficient when h, T, L or properties differ. The relevant comparison is the whole melt function.

## 3. Bare-to-bypass ordering

Suppose z_a(0)<=z_b(0). For any common eta,

    q_a=K_0[eta-z_a]_+^(3/2)
       >=K_0[eta-z_b]_+^(3/2)=q_b.

The lemma gives m_a>=m_b, and hence the gap obeys

    (z_a-z_b)'=-m_a+m_b<=0

whenever z_a<=z_b. Thus z_a(t)<=z_b(t) is invariant and its nonpositive gap cannot increase under the model. Strictly lower bare corridors tend to deepen at least as quickly, reinforcing hydraulic access. A time-dependent common stage caused by the covered route does not invalidate the comparison: the two bare per-width fluxes are compared at that same stage at each instant.

This proves bare-to-bypass ordering. **By itself it does not prove the bare toe is wet:** a third, even lower route could carry all I.

Counterexample to the stronger statement without additional hypotheses: z_a=z_b=0, z_c=-10 m, with sufficiently small I>0. The covered route alone conducts and z_a<=z_b holds while Q_a=0. The common stage comparison is not a theorem that every route receives flow.

## 4. Strengthened no-disconnection result

Add z_c(0)>=z_a(0), h_c<=h, identical T,L and K_0, and nondecreasing cover thickness. If z_c>=z_a then q_c<=q_a. The melt lemma gives m_c<=m_a, so

    (z_c-z_a)'=m_a-m_c+(cover-thickness increase)>=0.

Consequently z_c>=z_a remains true. Together with z_b>=z_a, the bare toe is a lowest route at all times. Since I>0, its common stage must satisfy eta>z_a; otherwise every route would have Q_i=0, contradicting the water budget. Therefore

    Q_a(t)>0

throughout the admitted interval. The covered route changes eta and total partition, but cannot dry the connected bare corridor in this model.

If the toe width tends to zero, this statement gives no quantitative useful flow lower bound. For finite width, the partition can be bounded from the controls; wetness alone does not guarantee environmentally important heating or a surviving cliff.

## 5. Equal covered area and equal stock do not specify contact

Consider a rectangle W by L with a fixed debris height H>0 and coverage fraction 0<f<1. Compare two arrangements:

1. A longitudinal covered strip of width fW and length L. A bare corridor of width (1-f)W remains connected from entrance to exit.
2. A transverse covered strip of width W and length fL. It forms a complete insulating/raised cutset across the toe route.

Both have covered planform area fWL and the same debris bulk volume H fWL. Their mean cover fraction, mean layer thickness over the whole rectangle, and total stock are identical.

In the ideal route construction, arrangement 1 has a bare route at floor 0 and a bypass at floor 0. For every I>0 the bare route carries water by the theorem above, regardless of the raised covered strip.

For arrangement 2, let the controlling transverse sill be at H and the bypass at 0. Its common-stage law is

    I=K_toe[eta-H]_+^(3/2)+K_b[eta]_+^(3/2).

Whenever

    0<I<=K_b H^(3/2),

the stage eta=(I/K_b)^(2/3)<=H and the entire toe route is dry. Thus the same covered area and stock admit opposite hydraulic contact classes. A closure using only those means cannot uniquely determine contact, even without stochastic noise.

The result is an explicit **information-insufficiency counterexample**, not proof that the same dichotomy holds for every natural three-dimensional geometry. The transverse strip must actually be a hydraulic cutset, without side leakage, porous underflow, overtopping roughness or an alternative bare path. Debris of finite permeability may carry water below its crest. Ice and debris evolution may destroy the cutset. Those processes need independent treatment.

## 6. What would make this useful environmentally

The useful candidate closure is based on **connected bare access / minimum cross-stream barrier height**, alongside debris amount, rather than mean covered fraction alone. A computational test should hold total stock, area fraction, thermal properties and supplied heat/water budgets fixed, vary only spatial arrangement, and compare contact duration plus integrated ice loss against a resolved independent model.

A result could show that a mean-cover model confuses a protected bare corridor with a transverse disconnection barrier. That would be a concrete model limitation and a measurable forecast improvement if it survives uncertainty and environmental-scale conditions. The theorem alone does not establish that improvement. It also does not imply that connectivity is a previously unknown concept in hydrology or percolation.

## 7. Important counterexamples and transfer limits

- A warmer or longer bypass can melt below the bare toe and capture flow; equal thermal forcing/length is essential.
- A rougher bare toe with a smaller per-width conductance can carry less flow despite a lower floor; identical K_0 is essential.
- A covered lane starting below the bare route can monopolize the flow.
- Erosion of its cover or independent ice motion can destroy the covered-above-bare ordering.
- A zero-width bare path is not a finite-discharge corridor; f=1 is excluded.
- Finite-depth sidewalls, lateral routing, pressure coupling, turbulence and grain blockage are not resolved by independent parallel routes.
- The theorem concerns contact preservation, not a guaranteed positive lower bound on seasonal melt or glacier-wide heat consumption.
- A stationary wetted contact can still be almost fully insulated, and a disconnected toe can still melt atmospherically.
