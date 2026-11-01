# Hydraulic disconnection as a falsifiable ice-cliff survival mechanism

Independent derivation, 31 October 2026. This does not alter Claude's work. This is a candidate mechanism, not a claim of historical priority or an environmental discovery already demonstrated.

## 1. Why another shielding model cannot suffice

The already-proved scalar stock model has `s'=S0+D F(T,s)-C(s,Q)`, with `F_T>0`, `F_s<0`, `C_s>=0`. Its steady response cannot decrease with warming, and its cumulative melt cannot reverse under pointwise warmer forcing. Therefore changing an arbitrary shielding exponent is not a legitimate new mechanism. A state must change a physical flux that was constrained monotone in the old model.

The candidate is **hydraulic disconnection caused by toe aggradation**. The stream can bypass the cliff. Debris accumulating at the toe raises its bed, sends water down the bypass, and thereby reduces both toe transport capacity and warm-water contact. Existing sediment-routing and river-avulsion literature must be checked before any novelty claim; the possible contribution is a quantitatively predicted *ice-cliff exposure transition*, not a first theory of avulsion.

## 2. Conservation-derived geometry and hydraulics

Let s be retained solid-grain volume within a toe strip of width W and length L; phi is porosity. In the simplest uniform aggradation geometry,

    z_t(s)=z_t0+s/[(1-phi) W L] = z_t0 + a s.

This relation requires bounded lateral spreading and must be replaced by a resolved wedge/control-volume relation for a real cliff. Do not equate stock with bed elevation without that geometric assumption.

Represent a toe route and a bypass with a common head H, route bed/sill heights z_t(s), z_b, and measured or independently computed discharge curves:

    Q_t = K_t [H-z_t(s)]_+^p,
    Q_b = K_b [H-z_b]_+^p,
    I = Q_t+Q_b.

Here I is supplied streamflow and p>1 (e.g. p=3/2 for a declared broad-crested control). This is a quasi-steady network surrogate. It needs a hydraulic adjustment time much shorter than the stock evolution time and no critical backwater changes outside its domain. It conserves water exactly. It is not full CFD and does not resolve turbulent grain stress.

When both routes are active define A=p K_t(H-z_t)^(p-1), B=p K_b(H-z_b)^(p-1). Implicit differentiation gives

    dH/ds = a A/(A+B),
    dQ_t/ds = -a A B/(A+B) < 0,
    dQ_b/ds = +a A B/(A+B).

Thus declining toe discharge is derived from routing, not inserted as an arbitrary declining export closure.

The toe is completely dry once

    z_t(s) >= z_b + (I/K_b)^(1/p).

The resulting disconnection stock is

    s_dis = (1-phi) W L [z_b + (I/K_b)^(1/p) - z_t0].

A nonpositive value means the toe route is dry initially. This stock threshold is exact within the declared route law. It depends on discharge, bypass geometry, and the stock-to-elevation geometry, not on a fitted shielding coefficient.

## 3. Conservative toe and ice budgets

Use retained solid volume, not bulk volume:

    s' = S_ext + D F(T,s,Q_t) - C(s,Q_t).

F is removed ice volume per time; D is released solid volume per removed ice volume, with separate mantle and englacial terms if their geometry differs. External material includes no material already in D F. Exported solid volume integrates C and closes independently. An admissible export law must satisfy C(0,Q)=0, C(s,0)=0, and C_Q>=0. It should be obtained from resolved grain motion, experiments, or an independently bounded relation. The disappearance of C as Q_t->0 follows from hydraulic loss of contact, even if C increases with available stock at fixed Q_t.

Ice loss must distinguish atmospheric exposed melt, stream melt and covered melt. A defensible placeholder for theory only is

    F(T,s,Q_t)=F_atm(s) + F_stream(T,s,Q_t),
    F_stream <= rho_w cp Q_t T/(rho_i Lf).

The upper bound is supplied water enthalpy, not a heat-transfer prediction. Setting F_stream equal to the bound assumes perfect extraction and is not acceptable as an environmental point predictor. Every term needs its corresponding energy source. Actual stream heat extraction should come from Claude's verified flow–thermal–moving-interface calculation.

## 4. Exact tipping test, without inventing a bistable closure

Let G(s)=S_ext+D F(T,s,Q_t(s))-C(s,Q_t(s)). The total derivative contains a routing term:

    G_s = D F_s - C_s + (D F_Q - C_Q) Q_t,s.

Partial derivatives here hold Q_t fixed. Because Q_t,s<0, sufficiently strong transport sensitivity C_Q compared with melt-derived supply sensitivity D F_Q can make G_s positive. This is the physical feedback missing from monotone fixed-discharge stock models. It is not guaranteed: shielding F_s<0 and availability C_s>=0 oppose it.

An interior unstable threshold s_u exists only if measured/resolved G crosses from negative to positive, `G(s_u)=0`, `G_s(s_u)>0`. Then stocks below it are cleared toward a low-stock attractor or boundary, while stocks above it accumulate toward disconnection or a high-stock attractor. Bistability requires checking both endpoint dynamics and additional roots; a single positive derivative does not establish it.

A useful local dimensionless routing-feedback ratio at a root is

    R = [(C_Q-D F_Q)(-Q_t,s)] / [C_s-D F_s],

provided the denominator is positive. The root is unstable iff R>1. This is an exact rearrangement of the Jacobian, not by itself a new theorem or universal scaling. Its possible scientific use is to locate the boundary with independently computed derivatives and test it against coupled simulations.

## 5. A rigorous stock crossing certificate

If G(s,t)>=g_min>0 throughout an interval before s reaches a guaranteed disconnection stock s_dis,max, then

    t_disconnect <= (s_dis,max-s0)/g_min.

To remain disconnected under a bounded future inflow I<=I_max, it suffices that

    z_t(s) >= z_b + (I_max/K_b)^(1/p)

and stock does not decrease. This is a conditional persistence certificate. Atmospheric delivery can keep stock increasing after stream disconnection; if all sources also vanish, s remains constant in the no-export idealization. Rainfall, later stream overtopping, bank erosion, route migration, ice advection, and cliff collapse can invalidate persistence and must be represented or bounded.

## 6. Consequential prediction and decisive falsifier

Target output: **probability or bounded classification of cliff–stream contact persistence and integrated melt over a forcing interval**, not a grain feature map. Compare the same initial geometry and forcing with a strongest fixed-contact, mass-conserving shielding model and with a mobile-debris routing model. The new information would be a correctly predicted loss of contact and its seasonal melt consequence.

The decisive mechanistic falsifier is simple: if stock increases but the stream remains at the toe or export rises enough to clear it, the hypothesized feedback fails. It also fails as a major contribution if fixed-contact models already predict exposure/melt equally well, or if the effect is confined to physically implausible geometry/roughness.

Do not claim warming reversal unless a conservation-consistent pair of simulations produces lower integrated ice loss at higher supplied enthalpy and the difference exceeds closure/numerical uncertainty. A routing switch can reduce local cliff melt while increasing downstream melt; a glacier-wide claim requires accounting for where the diverted heat goes.

## 7. Observational/computational credibility

A computation-only study can verify existence, thresholds and environmental-scale parameter regimes for this mechanism. It cannot claim field prediction accuracy without independent observations. Fully resolved CFD/DEM is not necessary to prove a conditional threshold, but route discharge, grain export and heat extraction require independent validation or uncertainty bounds. The existing archive lacks co-located stream forcing and repeated toe stock, so it cannot validate this mechanism directly.

The useful contribution, if it survives testing and priority review, is a **routing-aware survival closure with an explicit physically measured switching variable**. It is more consequential than another static shielding coefficient because it changes whether the stream can deliver heat to the cliff at all.

## 8. Essential correction: moving ice bed creates accommodation

The fixed-ice-bed version above is insufficient for an environmental model. Let b(t)>=0 be vertical toe ice lowering. Then

    z_t=z_i0-b+a s,
    z_t'=-b'+a(S-C),     a=1/[(1-phi)WL].

Therefore aggradation capable of shutting off flow requires

    S-C > (1-phi) W L b'.

If the only supply is local englacial debris released by toe melting, define D as solid debris volume per pure-ice volume and A_toe=WL. Then S=D A_toe b', and

    z_t'=[D/(1-phi)-1]b' - a C.

For D<1-phi and C>=0 this is nonpositive. **Local melt release alone cannot raise the toe bed in this geometry.** Any model that predicts aggradation from low-concentration local debris while neglecting melt-created accommodation is physically wrong.

The required source is nonlocal mantle delivery or an independently represented collapse. For a retreating cliff lip at horizontal speed v_lip under a mantle of thickness d_top,

    S_mantle=(1-phi) d_top W v_lip chi,

where chi is delivered fraction. With no other sources, aggradation requires

    chi d_top v_lip/L > b' + C/[(1-phi)WL].

This is a measurable supply/accommodation condition. It shows why merely assuming “more melt releases more debris and blocks the stream” is not sufficient. Include a finite upper-mantle reservoir and its evolving source geometry; unlimited constant delivery is not conservative over long horizons.

## 9. Size-resolved stock and the coarse-tail question

Track solid stocks s_j for size/shape classes, released volume fractions f_j (sum=1), independent sources S_ext,j and exported volumes E_j:

    s_j' = f_j S_mantle + D_j A_toe b' + S_ext,j - C_j,
    E_j' = C_j,
    s=sum_j s_j.

All C_j must be limited by available class stock and independently justified grain force/transport capacity. A deep sediment-bed Shields relation is not automatically valid when grain diameter approaches flow depth. A force-based partial-submergence closure or resolved grain mechanics is needed.

Suppose a class subset A has no admissible entrainment under a declared forcing bound (for example a force-balance calculation proves U<U_c,j for every class j in A). Then its retained stock obeys exactly

    s_A(t)=s_A(0)+integral [f_A S_mantle + D_A A_toe b' + S_ext,A] dt.

Fine-grain export cannot remove that persistent coarse volume. A mean grain diameter can therefore miss the relevant *non-exportable tail*. This is an exact accounting consequence, not a novel sediment-armoring principle; priority must be checked.

If a resolved geometry analysis supplies a deterministic contact-disconnection stock s_dis(t), then sufficient accumulation requires s_A>=s_dis, but total volume alone does not guarantee a particular barrier shape without that geometry condition. The threshold must not be fitted to survival outcomes.

For a fixed bed geometry and a fixed critical retained coarse stock, cumulative mantle retreat volume required for threshold is proportional to 1/f_A. A small coarse tail can dominate lifetime even if median grain size remains fine. When the ice bed melts, use the full condition

    a s_A(t)-b(t) >= z_b+(I_max/K_b)^(1/p)-z_i0.

A coarse tail only creates rising topography if its delivery beats ice accommodation. This places a concrete environmental feasibility gate on the proposed mechanism.

## 10. Flood rescue window and budget-dependent outcomes

A high-flow event has three competing effects: it raises the head and may reconnect the toe, exports grains if their entrainment thresholds are crossed, and supplies more warm-water enthalpy that lowers ice and creates accommodation. These effects have different thresholds. Consequently a flood can rescue contact without removing all toe stock, or remove fine material while leaving an immobile coarse barrier.

For a pulse with bounded incoming discharge and independently bounded class capacities, a necessary clearing condition for reducing stock by Delta s is

    integral_pulse sum_j C_j dt >= Delta s + integral_pulse S dt.

For lasting reconnection, the post-pulse geometry must additionally satisfy

    z_t,post < z_b+(I_post/K_b)^(1/p).

Temporary overtopping during the pulse is not lasting rescue. A mathematically and environmentally useful target is the minimum pulse volume/duration/peak satisfying both the material budget and the post-pulse hydraulic criterion. If a class is immobile at every pulse velocity in the admitted domain, export-only rescue is impossible for that class, although ice-bed incision may still reconnect flow.

This optimization can be computed with deterministic uncertainty bounds. It is more useful than a generic demonstration of hysteresis: the output is a numerical intervention/forcing threshold with an explicit failure certificate. It remains conditional on independently verified capacity, geometry and heat exchange; no external agency suitability follows merely from implementing the optimizer.
