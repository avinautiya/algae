# Audit: finite-horizon no-initiation certificate for a coarse grain

31 October 2026. This is a conditional certificate for explicitly specified control, contact and heat-transfer closures, not a universal impossibility theorem for a real glacier stream.

## 1. Critical-depth force balance

Consider a stationary cube of side D on a horizontal plane with Coulomb sliding coefficient mu>0, constant drag coefficient Cd>0, grain density rho_s>rho_w>0, no lift, no pivoting/rolling, and quasi-steady drag. In an ideal rectangular critical-flow section, u²=g d. For 0<=d<=D, the submerged projected frontal area is D d and buoyancy volume is D²d. Sliding onset therefore satisfies

    0.5 rho_w Cd (g d) D d
       = mu g (rho_s D³-rho_w D²d).

Writing x=d/D and r=rho_s/rho_w gives

    0.5 Cd x²+mu x-mu r=0,
    x_c=(-mu+sqrt(mu²+2 Cd mu r))/Cd.

For numerical stability at small Cd use the equivalent positive expression

    x_c=2 mu r/[mu+sqrt(mu²+2 Cd mu r)].

If x_c<=1 this partial-submergence solution is valid. Otherwise full submergence gives

    x_c=2 mu(r-1)/Cd >1.

The branches meet at Cd=2 mu(r-1), x_c=1, so the piecewise law is continuous. Strict positive density contrast is required. The mu=0 limit is immediate mobility for any nonzero flow, not a useful no-rescue bound. Shape, orientation, friction and Cd need uncertainty bounds chosen to make mobilization as easy as physically admissible when constructing a conservative no-initiation test.

The force calculation is a sliding threshold only. If rolling, lift, collision, bed failure, undermining or mechanical collapse can move the grain at lower flow, this is NOT a bound against all motion. Either prove sliding is the easiest allowed mode or replace d_c with a lower bound across all admitted modes. Using an arbitrary high friction threshold would invalidate the claimed impossibility.

## 2. Discharge required at fixed route elevations

For two ideal broad-crested controls,

    Q_1=K_1[eta-z_1]_+^(3/2),
    Q_2=K_2[eta-z_2]_+^(3/2),
    Q=Q_1+Q_2,

and define gap g=z_1-z_2. For the toe route to reach depth d_c at its critical section requires head H_c=3d_c/2 above its crest, because critical depth d=2H/3 in the ideal loss-free broad-crest relation. The required total discharge is

    Q_move(g)=K_1 H_c^(3/2)+K_2[max(g+H_c,0)]^(3/2).

It is continuous and nondecreasing in g, strictly increasing when the bypass is active. For nonnegative total inflow, the common-stage equation has a unique stage whenever at least one route is conducting; evaluating it at eta=z_1+H_c gives exactly Q_move.

Consistency warning: a loss-free rectangular control has K=(2/3)^(3/2) W sqrt(g_accel). If K includes an empirical discharge coefficient or significant loss, d=2H/3 and u²=g_accel d are not automatically simultaneously valid. The force depth and the discharge curve must be calibrated/derived consistently. Avoid using g for both elevation gap and gravitational acceleration in code.

If the grain is not actually at the hydraulic control, local depth/velocity do not follow from crest head. A verified transfer from route discharge to the grain's local force is required. At obstacles comparable to depth, neglecting flow blockage and pressure forces can be particularly inaccurate.

## 3. Optimistic lower bound on the route gap before first motion

Let S_c>0 be retained coarse solid volume that cannot leave before the first coarse-grain motion event. Let A_s=(1-phi)WL>0 convert retained solid volume into uniform bulk thickness. At every time before that event, even if all fine material has already been removed,

    z_1(t)>=z_ice,1(0)-b_1(t)+S_c/A_s.

This requires the fixed toe-footprint, uniform-stock-to-elevation geometry. A coarse grain can remain in the control volume yet rearrange or roll into a lower position; in that case constant S_c does not guarantee constant route elevation. The certificate must exclude such rearrangement or use a geometry-derived lower bound, not volume alone.

Suppose the upper bound on toe ice lowering rate is

    b_1'(t)<=M_max
      =T_max/[rho_i L_f (1/h_max+S_c/(A_s k_d,max))].

This follows only from a fully covering uniform series-resistance layer with no lateral bypass of heat, no grain heat storage release beyond the prescribed budget, no radiation/atmospheric addition beyond what is included in T_max,h_max, and a temperate ice latent-heat budget. Use the largest admissible debris conductivity and heat-transfer coefficient to maximize melt. A current stored-hot-grain or conductive heat reservoir can invalidate this instantaneous rate cap even if a cumulative energy cap remains available.

A safer generalized form is simply b_1(t)<=B_max(t), independently established from all admitted heat sources and initial storage. The resistance formula is one route to B_max, not a universal fact about debris cover.

If the bypass crest cannot rise, z_2(t)<=z_2(0). Bypass lowering increases g=z_1-z_2 and makes toe mobilization harder. Ignoring it is therefore optimistic for rescue and valid for the lower bound

    g(t)>=g_lower(t)
      =z_ice,1(0)+S_c/A_s-z_2(0)-M_max t.

Over 0<=t<=tau,

    g(t)>=g_lower(tau).

If the bypass can aggrade, jam or receive debris, that can increase its elevation and favor toe routing. Its maximum possible rise must be subtracted from g_lower. Ignoring bypass rise is not conservative.

## 4. First-crossing proof

Assume all the closures above hold until the first coarse-motion time t_star, and inflow satisfies Q(t)<=Q_peak. Suppose

    Q_peak < Q_move(g_lower(tau)).

For contradiction, assume 0<t_star<=tau is a first motion event. Before it S_c is retained, so by continuity the gap bound holds at t_star. The function Q_move is nondecreasing; hence

    Q(t_star)<=Q_peak
       <Q_move(g_lower(tau))
       <=Q_move(g(t_star)).

The stage at t_star is therefore strictly below the stage needed to reach d_c. The assumed force law gives drag below sliding resistance, contradicting first motion. Thus no admitted coarse sliding can initiate during that finite interval.

The strict inequality matters. Equality is unresolved and cannot support a strict no-motion result. If Q_peak=0, the result is trivial when d_c>0. If g_lower+H_c<=0, the bypass is inactive even at onset; then Q_move reduces to the toe-only threshold and remains valid. Negative gaps are not invalid inputs. Large tau weakens the result as the optimistic ice lowering grows; failure of the inequality is not evidence that rescue occurs.

A pre-existing coarse motion at t=0, or an initially unstable contact, violates the first-crossing premise. A time-discontinuous inflow can be handled using the right-limit stage and the same bound; an impulsive force or unsteady inertia is outside quasi-steady drag.

## 5. What is and is not established

The theorem is a finite-horizon impossibility result inside a stated model: **no initial coarse sliding even after instantaneous removal of all fines and the maximum accepted toe melt**. It can supply a useful lower-bound rejection test for a forecast/optimization tool. Its environmental value would be demonstrated if it robustly rejects proposed flood-rescue histories that a mean-grain model predicts as successful, with the rejection confirmed by an independent resolved model.

It does not imply no reconnection: melting can reconnect flow even if coarse grains never move, and a grain can be bypassed geometrically. It does not imply no cliff survival, no downstream transport, or no glacier-wide warming response. Grain rolling, lift, rearrangement, ice fracture, bypass aggradation and three-dimensional lateral rerouting must be excluded or bounded before applying the result environmentally.

Historical originality has not been established by this audit. The mathematical proof combines monotone hydraulics, a force threshold and an optimistic geometric bound; the potential new scientific contribution is the environmental prediction and robust exclusion region, not the individual elementary inequalities.
