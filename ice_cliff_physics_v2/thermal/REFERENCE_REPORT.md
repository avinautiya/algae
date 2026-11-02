# Independent closure check: the one-node model is not transiently equivalent to a debris slab

All cases are synthetic, with fixed shared water-exchange closure and no parameter fitting. The reference resolves only normal diffusion. These results are not environmental field accuracy or a claim of historical mathematical novelty.

## Main matched-parameter comparison

Physical inputs: d=0.05 m, k=1 W/(m K), rhoc=1.8e6 J/(m3 K), A=1 m2, Tin=0.5 C, F=400 W/K. Initial debris is at0 C. The lumped model uses exactly C=rhoc*A*d=90000 J/K and K=k*A/d=20 W/K. Diffusion time is4500 s. No effective capacity is tuned.

At one diffusion time:

| Nominal h d/k | Resolved interface heat (J) | Lumped interface heat (J) | Lumped relative error |
|---|---:|---:|---:|
| 0.1 | 2303.122 | 1605.944 | -30.27% |
| 1 | 14902.984 | 12531.860 | -15.91% |
| 10 | 31774.619 | 35425.409 | +11.49% |

At early times reference ice heat is close to zero; huge relative ratios there are not useful accuracy statements. The CSV additionally supplies an absolute difference normalized by steady flux times elapsed time, avoiding division by a nearly-zero reference.

## Independent analytic dry-storage limit

For an initially uniform hot slab, bottom temperature0 C and insulated top, the fraction of initial energy remaining is

    sum[n=0..infinity] 8/(pi²(2n+1)²) exp[-(2n+1)² pi² t/(4 tau)].

The dry pulse uses initially4 C debris with initial energy360000 J.

| t/tau | Analytic ice heat (J) | Resolved ice heat (J) | Lumped ice heat (J) |
|---|---:|---:|---:|
| 0.01 | 40621.650 | 40305.248 | 3582.060 |
| 0.1 | 128456.424 | 128356.361 | 34258.530 |
| 1 | 335253.484 | 335211.511 | 227563.401 |

## Refinement and verification

At fixed dt/tau=0.00025 and t/tau=0.2, doubling N=16,32,64,128 gives successive heat differences 2.30102, 0.576418, 0.144177 J; ratios 3.992, 3.998 agree with second-order space convergence.

At N=128 and t/tau=0.2, halving dt/tau from0.01 gives successive heat differences 26.8199, 13.6807, 6.91029, 3.47291 J; ratios 1.960, 1.980, 1.990 agree with first-order time convergence.

Twenty unit checks cover the exact wet-contact kernel, independent spatial-profile/RK4 quadrature, passivity, dry and stiff limits, exact steady Robin flux, analytic Dirichlet diffusion, analytic dry-storage release, time order and areal-coordinate equivalence. A separate5000-case seeded numerical-conditioning stress test has no failures; its maximum normalized one-node energy residual is4.08e-13. It tests arithmetic and declared equations, not material-property accuracy over that extreme parameter range.

## Scientific decision

The conservative one-node kernel is useful as an explicitly approximate closure and exact RC reference. Its perfect energy ledger is not evidence that its transient ice heat is accurate. For transient thermal-memory claims, use the resolved normal slab with independent mesh/time checks, or demonstrate an approximation error smaller than the effect being claimed. Small water-side h d/k alone does not establish a uniform-temperature approximation when the opposite side is held at melting temperature.

The slab still assumes homogeneous material, a fixed normal thickness/areal constitutive pair, temperate ice, quasi-steady water exchange and areally homogenized wet/dry forcing. It does not resolve lateral wet/dry gradients, debris rearrangement, liquid storage within debris, refreezing, ice conduction or a moving melting interface. None of these limitations disappears because a ledger closes.

Ten additional wet/dry-column partition checks verify conservative per-cell transfers and a coupled heat/remapping sequence. The original one-node/slab/partition suite totals30 tests; ten additional analytic Robin tests are documented in ROBIN_REFERENCE.md.
