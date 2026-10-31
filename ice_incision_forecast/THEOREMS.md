# Sharp regional incision bounds under incomplete thermal information

Proof manuscript prepared 31 October 2026. The results below are proved for the stated reduced model. The precise combination of a meltwater-corrected map, sharp spatial bounds and forcing-uniform witnesses was derived for this project. Priority has not been certified; the source comparison is recorded separately. A proof establishes a mathematical implication, not that all its assumptions hold in a glacier.

## Model and notation

Consider a fixed reach x in [0,L] with prescribed positive inlet liquid discharge Q0 and inlet temperature excess theta0>=0 above the interface melting point. Let P=rho_w cp Q theta be advected sensible power and let a(x)>=0 be an effective heat conductance per unit length, in W/(m K). All water-to-ice heat becomes latent heat in temperate ice and its meltwater joins the liquid at zero temperature excess. There is no additional heating, conductive loss into cold ice, retained thermal storage or lateral water source. Properties are constant. Thus

```
dP/dx = -a(x) theta
dQ/dx =  a(x) theta/(rho_w Lf)
```

These equations include meltwater addition; they do not assume downstream Q is equal to inlet Q. They are a quasi-steady thermal model, not the Navier–Stokes/Stefan free-boundary system.

Write G(x)=integral_0^x a(s) ds, C0=rho_w cp Q0, P0=C0 theta0, and Ste=cp theta0/Lf. For a spatial partition let g_i=integral_cell_i a dx, with independently declared bounds l_i<=g_i<=u_i. A total conductance band and optional exact prefix conductance anchors define the admissible class. No probabilistic interpretation is imposed.

## Theorem 1 — nonlinear thermal map and observational equivalence

For theta0>0 and every finite G>=0, there is a unique s>=0 satisfying

```
G/C0 = s + Ste (s - 1 + exp(-s)).                         (1)
```

The solution obeys

```
G/[C0(1+Ste)] <= s <= G/C0,
P(G)=P0 exp(-s),
Q(G)=Q0 + [P0-P(G)]/(rho_w Lf),
theta(G)=theta0 exp(-s)/[1+Ste(1-exp(-s))].              (2)
```

P is strictly decreasing and strictly convex. In particular, any two admissible conductance fields with the same total G(L) have the same outlet power, discharge, temperature and integrated melt. This equivalence holds for every inlet Q0>0 and theta0>=0 in the model. If a(x,t)=gamma(t) a0(x), gamma>=0, it holds for every prescribed inlet history under the quasi-steady interpretation. Outlet histories alone cannot identify how conductance, and hence local melt, is distributed within an unobserved spatial block.

### Proof

The governing equations give d[P+rho_w Lf Q]/dx=0. Therefore Q=Q0+(P0-P)/(rho_w Lf). Where conductance is positive, divide by dG/dx=a to obtain

```
dP/dG = -P/[rho_w cp Q].                               (3)
```

Integrating with P=P0 exp(-s) gives (1). Its right-hand derivative with respect to s is 1+Ste(1-exp(-s)), which lies between 1 and 1+Ste for s>=0. Existence, uniqueness and the bracket follow. Substitution gives (2).

Let C=rho_w cp Q=C0+(cp/Lf)(P0-P)>0. Then P'=-P/C<0 and C'=-(cp/Lf)P'>=0, so

```
P''=P/C^2 + P C'/C^2 >0.                              (4)
```

All endpoint quantities depend on a only through G. Under a common nonnegative multiplier, replace G by gamma(t)G in each inlet-dependent map. Equal total base conductance therefore yields equal outlet histories. The case theta0=0 has P=0, Q=Q0 and zero melt, so the equivalence remains true. QED.

### Exact prefix information

For 0<theta_out<=theta0, put r=theta_out/theta0. Equations (2) imply

```
P_out/P0 = r(1+Ste)/(1+r Ste).                         (5)
```

This recovers s and hence cumulative G from an ideal bulk-temperature observation with the stated inlet information. A zero outlet excess supplies a limiting constraint, not a finite exact conductance. Point sensor temperatures are not automatically the required bulk observation.

## Theorem 2 — sharp regional bounds with explicit extremizers

Consider a target spatial interval inside one independent conductance block. Denote total conductance before, inside and after it by A,B,C. Sum the cell bounds in the three groups to obtain

```
A in [a_l,a_u], B in [b_l,b_u], C in [c_l,c_u],
g_l <= A+B+C <= g_u.                                  (6)
```

Assume this feasible set is nonempty. Heat delivered to the target is H(A,B)=P(A)-P(A+B), with the appropriate fixed upstream conductance offset included if the block follows an exact prefix anchor.

The **maximum** target heat is attained at

```
A_maxheat = max(a_l, g_l-b_u-c_u),
B_maxheat = min(b_u, g_u-A_maxheat-c_l).                (7)
```

The **minimum** is attained at

```
A_minheat = min(a_u, g_u-b_l-c_l),
B_minheat = max(b_l, g_l-A_minheat-c_u).                (8)
```

For either pair choose C=max(c_l,g_l-A-B). Every group total can be allocated within its individual cell bounds, yielding a complete attaining conductance profile. Hence these are sharp bounds over the stated admissible class, rather than conservative but unattained estimates.

### Proof

Strict convexity of P gives

```
partial_A H=P'(A)-P'(A+B)<=0,
partial_B H=-P'(A+B)>0.                               (9)
```

Feasible A ranges from max(a_l,g_l-b_u-c_u) to min(a_u,g_u-b_l-c_l). At a fixed A, feasible B lies between max(b_l,g_l-A-c_u) and min(b_u,g_u-A-c_l). Both endpoint functions are nonincreasing in A. Therefore the largest H uses the smallest feasible A and largest feasible B; the smallest H uses the largest feasible A and smallest feasible B. These give (7) and (8). The stated C is feasible by construction.

For a group with cell lower/upper sums l and u, any target v in [l,u] is realized by g_i=l_i+eta(u_i-l_i), eta=(v-l)/(u-l), with the fixed case u=l treated separately. Applying this in all three groups constructs the witnesses. The inequalities are non-strict when a group or heat input is zero. QED.

## Theorem 3 — forcing-uniform finite-horizon bounds

Suppose a(x,t)=gamma(t)a0(x), gamma(t)>=0, and all inlet scenarios use positive Q0(t), nonnegative theta0(t), and the same admissible base-profile class. Then the two spatial witnesses in Theorem 2 attain the minimum and maximum **time-integrated regional heat** over any finite horizon for which these integrals exist. They do not depend on the inlet forcing history or the common multiplier.

### Proof

At each time the corresponding P_t is decreasing and convex by Theorem 1. Replacing A,B by gamma(t)A,gamma(t)B preserves the monotonicities in (9). The feasible optimizing group totals (7)-(8) are independent of P_t, Q0, theta0 and gamma. The same two profiles therefore attain the pointwise-in-time extrema. Integrating inequalities over time preserves order, and integrating those witness profiles attains both bounds. QED.

For fixed bed-only melting over planform area S_R, regional ice-volume loss is E_R/(rho_i Lf) and mean vertical lowering is E_R/(rho_i Lf S_R). The resulting finite-horizon lowering bounds inherit sharpness. This conversion does not predict lateral wall erosion or resolve hydraulic/geometry feedback.

## Corollary — exact-prefix partitions and thermal observation requirements

Exact prefix conductances partition the admissible profile into independent blocks. If both ends of a target lie in one block, apply Theorem 2 with its fixed upstream offset. If they lie in distinct blocks, their cumulative conductances have independent feasible ranges; decreasing P gives the extremizers by pairing the smallest start prefix with largest end prefix for maximum heat, and the opposite for minimum. Witness allocations in distinct blocks are compatible.

If exact prefix conductances are known at **both target boundaries**, regional heat is uniquely determined in this model, although its within-region spatial distribution remains uncertain. Additional outlet observations do not replace those interior spatial constraints.

### Proof

Each exact prefix fixes the sum up to that boundary; subtract consecutive values to fix each intervening block total. Box constraints in different blocks then decouple. A cumulative prefix inside a block has minimum max(prefix lower sum, block lower total minus tail upper sum) and maximum min(prefix upper sum, block upper total minus tail lower sum). These extremes have explicit feasible allocations. Monotonicity of P proves the different-block result. When both endpoints are fixed, P at both endpoints is fixed. QED.

## Theorem 4 — optimal worst-case point forecast and threshold timing

Let [d_l(t),d_u(t)] be the sharp bounds on regional mean lowering supplied by Theorem 3. Among deterministic scalar forecasts based on the same information, the minimax forecast under absolute error is

```
d_star(t)=[d_l(t)+d_u(t)]/2,
R_star(t)=[d_u(t)-d_l(t)]/2.                            (10)
```

No scalar forecast can guarantee a smaller worst-case absolute error over this uncertainty class. If a positive declared incision threshold d_c is used, its earliest possible crossing is the first time d_u reaches d_c, and its latest crossing is the first time d_l reaches d_c. Each finite bound is attained by its corresponding spatial witness. Failure of d_l to cross within the forecast horizon means crossing is not guaranteed within that horizon, rather than that it never occurs.

### Proof

For any scalar forecast z, the two attaining witnesses imply worst-case loss at least max(|z-d_l|,|z-d_u|), which is at least (d_u-d_l)/2. The midpoint achieves that value for every admissible output between the bounds. This is the classical scalar interval-minimax argument applied to the sharp physical-model bounds, not a new general decision theorem.

All melt rates are nonnegative, so regional lowering is nondecreasing. By Theorem 3 the same spatial witness attains each envelope over the whole horizon. Before d_u first crosses, no profile can have crossed; its witness attains that earliest crossing. When d_l crosses, every profile has crossed, and its witness attains the latest crossing. QED.

This result concerns a declared incision observable. A hydraulic-failure or flood threshold requires separate physical justification and coupling. Scalar minimax forecasts for several regions need not correspond to a single conservative spatial profile; the tool retains a separate feasible nominal profile for its joint mass/energy ledger.

## What sharpness does and does not mean

Sharpness is with respect to the conductance uncertainty class in (6). Some extremizers may not be realizable by a particular Navier–Stokes geometry, turbulence model or thermal history. Additional physical constraints can shrink the class and its bounds. Coverage of a real forecast requires the true effective conductance and boundary conditions to lie in the assumed class; it is not established by the algebra alone.

The code includes meltwater mass addition but omits water heat storage, extra sources, cold-ice conduction and evolving hydraulics. Common-multiplier separability must be verified before using the forcing-uniform statement. The theorems do not prove a density-anomaly mechanism, a hydraulic instability or superiority over other models.
