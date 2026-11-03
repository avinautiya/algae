# Temperature-fit and basal-heat diagnosis

All calculations here use the already inspected Tower1 **TRAIN block as
development data**. There is no new validation, no Tower2 camera access and no
change to any previous failed field-comparison verdict.

## What the executable comparison establishes

The source-proxy capacity is1.855487MJ/m³/K; use the unchanged middle20mm gap
geometry,32debris/64ice cells,300s steps, known ice properties and a cold-content
complementarity interface. Search a fixed8×8 positive conductivity grid. All
models are passive conductors; positivity is physical admissibility, not evidence
that the conductivity bounds are independently measured at this site. The
two-layer boundary is fixed at half of the computational debris thickness; it
is a candidate model assumption, not observed stratigraphy. Every face integrates
the piecewise resistance exactly, even where a layer boundary crosses a cell.

| TRAIN-temperature-selected candidate | k upper/lower W/m/K | Raw TRAIN melt mmwe | Gain needed to match known TRAIN mass |
|---|---|---:|---:|
| Homogeneous, plain MSE |4.5/4.5|371.78|0.229|
| Homogeneous, depth-standardized loss |3.5/3.5|288.23|0.296|
| Two-layer, plain MSE |4.5/2.5|264.03|0.323|
| Two-layer, depth-standardized loss |4.5/2.5|264.03|0.323|
| Two-layer, basal T/gradient objective |3.5/3.5|288.23|0.296|

The already inspected TRAIN target is85.20165mmwe. It is read **after** all
temperature-based selections, solely to diagnose amplitude. Neither new layering
nor the tested goal-focused weighting eliminates that mismatch. The plain-MSE
two-layer candidate reduces its temperature MSE from3.4835 to2.4292C² and raw
melt from371.8 to264.0mmwe, but still needs a factor0.323 output correction.
Upper conductivity lies at the search boundary, so the finite grid does not
establish a globally unique optimum.

None of the64 tested pairs matches every observation within the **linearized
engineering box**`.4C + .02m*|observed local TRAIN gradient|`. That box is an
approximate depth-error propagation, not an exact joint uncertainty set or a
statistical interval. Failure of a finite grid does not prove that every possible
heterogeneous material model is impossible. It does rule out presenting these
specific fitted candidates as already physically validated.

Depth-specific derivative-OLS estimates using the same TRAIN thermistors are
1.3551e−6,5.2705e−7 and1.8510e−7m²/s atDT2,DT3 andDT4. Their conditional
conductivities are2.5143,.9779 and.3434W/m/K. Differentiation amplifies sensor
noise; these are diagnostics of the homogeneous closure, not independent
material measurements. The TRAIN mean temperatures fall strongly betweenDT3
andDT4, so a fixed half-depth two-layer partition is not guaranteed to represent
the relevant bottleneck. A temperature misfit is not proof of a particular
missing geological or interfacial process.

## Conservative layered equations

In each debris layer, `C T_t = (k(z) T_z)_z`. For an interval `[za,zb]`, its
conductance is `G=[∫za^zb dz/k(z)]^−1`. The cell-centered implicit equation is

```
(C Δz/dt)(Tj_new−Tj_old)
  = Gleft(Tleft_new−Tj_new) − Gright(Tj_new−Tright_new).
```

Summing cancels all internal face fluxes. Add the cold-ice domain and interface
latent mass to obtain the same global sensible-plus-latent ledger as the frozen
phase model. Tests cover exact unaligned two-layer steady series resistance and
the homogeneous limit against the frozen implementation. The unscaled raw
component is conservative. An empirical output gain is a different statistical
model and cannot inherit that ledger.

## Exact scale obstruction to absolute heat identification

When the ice interface stays at0C, replace all debris conductivities andC by
the same positive factor`s`. Every debris equation is multiplied bys, so its
temperature solution is unchanged. Known ice conductivity/capacity do not help
in this regime: initially temperate ice stays at0C and carries no sensible flux.
All debris fluxes and melt scale bys. This is an exact symmetry, not a numerical
optimization issue.

The selected two-layer TRAIN candidate has **zero cold-interface steps**. Across
the independently declared dry/proxy/saturated capacities, holding layer
diffusivities fixed produces temperature histories identical within3e−13C, while
raw melt is170.57/264.03/432.21mmwe. Even the dry member of this particular
family remains above the known85.20mmwe TRAIN target. That does not make the
source capacity wrong; geometry, missing processes and model discrepancy remain
possible, and the continuous model family has not been exhaustively excluded.

The homogeneous selected TRAIN candidate also has zero cold-interface steps.
Thus this particular TRAIN interval provides **no cold-ice admittance signal**
for identifying absolute heat scale, despite including known ice properties in
the software. Other colder records may break the symmetry, provided ice initial
state and deep-boundary uncertainty are also addressed.

## Local observability and the adjoint

Let `F_n=A_n(theta)T_{n+1}−B_n(theta)T_n−f_n(theta)=0`, wheretheta includes
material and geometry parameters. On a fixed active-set sequence, the tangent is

```
A_n S_{n+1} = B_n S_n + (dB_n)T_n + df_n − (dA_n)T_{n+1}.
```

For a basal-heat goalG, the discrete adjoint with terminal condition is

```
A_(N−1)^T lambda_N = dG/dT_N,
A_(n−1)^T lambda_n = dG/dT_n + B_n^T lambda_(n+1).
```

Then `dG/dtheta = partial_theta G − Σ lambda_(n+1)^T partial_theta F_n`, with
the initial-state contribution retained if it depends ontheta. Switching active
sets creates piecewise differentiability; a global smoothness claim is invalid.
This study uses finite differences, not an implemented exact adjoint.

For observation sensitivityJ and goal gradientb, a local temperature-only goal
estimate requires an observation combinationa satisfying`J^T a=b`. If a parameter
directionv has`Jv=0` but`b^T v!=0`, no observation weighting can recover the goal:
the data simply lack that information. Under a linearized error box
`|delta_y_i|<=epsilon_i`, the best linear combination has the familiar bound

```
min_a Σ epsilon_i*|a_i|, subject to J^T a=b.
```

If the constraint is infeasible, independent parameter bounds or another
measurement type are needed. This is standard linear observability/duality,
not a priority claim to a new theorem. In the warm scaling symmetry,
`v=(1,1,1,0)` for`logk_upper,logk_lower,logC,probe_depth`, while the log-melt
derivative alongv is1: absolute basal heat is not observable from temperatures.

The numerical Jacobian initially appeared full rank with smallest singular value
1.78e−4, but halving the finite-difference step reduces it to4.45e−5 and then
1.11e−5, precisely theh² truncation pattern. Exact simultaneous-scaling tests
give machine-precision invariance. Interpreting that tiny singular value as
physical absolute-scale information would be an error. Larger singular values
are stable under step refinement. The reported condition/amplification numbers
are engineering-scaled local diagnostics, not confidence intervals.

## Conservative integrated weak-form alternative

A useful next development is a positive thermal-network inverse on the actual
probe locations, rather than selecting conductivity by bulk temperature MSE.
For P1 finite elements on measured nodes, include the **consistent** mass matrix:

```
Mii = C(hleft+hright)/3,
Mi,i−1 = C hleft/6, Mi,i+1 = C hright/6,
Σj Mij [Tj(tb)−Tj(ta)]
  = Gleft ∫(Tleft−Ti)dt − Gright ∫(Ti−Tright)dt.
```

The measured Dirichlet-top temperature's mass-coupling term must not be dropped.
Infer positive conductances from TRAIN integrals only, report matrix rank and
errors in both regressors and endpoint temperatures. Overlapping windows are
correlated; do not invent independent-sample confidence intervals.

For an ice Dirichlet node at0C, the energy-consistent outward boundary reaction
includes last-element storage:
`Qice=Glast∫T5 dt − C*hlast*[T5(tb)−T5(ta)]/6`.
Using only its stiffness gradient drops that storage term. The sameC/G scaling
ambiguity still exists unlessC or another absolute energy scale is independently
constrained. A positive NNLS solution alone is neither identified physics nor
field validation. The parent is implementing this alternative separately; this
document does not claim its result.

## Decisions supported by this development study

- Do not launch a new melt test merely because layering or weights were added.
- Separate better TRAIN temperature fit from raw energy/ablation agreement.
- Do not interpret numerical Jacobian truncation as cold-interface observability.
- Retain source-capacity/geometry ambiguity and test omitted-process hypotheses
  with independent observations; do not assign a unique field cause from these data.
- Keep Tower2 outcomes unopened until a new candidate passes its identification
  and physical-comparison gates and predictions are frozen.

The contribution here is executable physical diagnosis and a demonstrated
identification obstruction. A new environmental mechanism or superior predictor
has not been established by this round.
