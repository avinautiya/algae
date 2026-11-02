# Version5 phase-mass repair: implemented and bounded verification

## Isolation and provenance

The exact source under audit is vendored as `swe2d_v4_pinned.py` (Git blob `d0d25b0c5f92aa4bb667db577e5f87f75b3ed086`, repository commit `05e7312e3bbf248acf9f646267a5b8ea06a8b32e`). It is unchanged. The new implementation is `swe2d_v5.py`, class `Sim`, a subclass using that same actual Numba hydraulic kernel. No Claude checkout, frozen output, job, or historical verdict was edited. This is an additive candidate implementation, not a production replacement.

## Repair

For an interior bare-ice cell melting vertically by `dz`, the new source applies

```
zb_new = zb_old - dz
h_new  = h_old + (rho_ice / rho_water) * dz
hu_new = hu_old
hv_new = hv_old
hT_new = hT_old
```

Ice is stationary and temperate. Meltwater enters at the melting point, with zero horizontal momentum and zero sensible heat relative to that point. The old conserved momenta and depth-weighted temperature are therefore retained; velocity and temperature dilute. The free surface drops only by `(1-rho_ice/rho_water)*dz` during an isolated source step, rather than by `dz`. Ghost cells and walls are excluded from melt, generated volume and region ledgers. No ice-floor clamp is reintroduced.

### Two time semantics are kept separate

1. **New `step_physical(dt, cool=True)`** advances the actual hydraulic/thermal kernel for `dt`, then uses the same heat integral and duration to lower local bare ice and add water. It updates `physical_time_s` and `t`. This is a first-order split shallow-water/thermal/Stefan floor model. It is not VOF or a full mechanical-energy/ALE solver.
2. **Existing `hydraulics(...)` followed by `lower_bed(dtm)`** remains the explicitly labelled quasi-steady, morphologically accelerated approximation. Hydraulic spin-up seconds are not environmental elapsed seconds. Each source records `heat_time_scaling=dtm/window_time` and the difference between its morphological latent budget and the actual last hydraulic window's heat removal. The software refuses to issue a combined physical-time thermal/latent certificate for this mixed-time path. Merely executing a public settling call also disables that certificate for subsequent history; it cannot silently be relabelled physical time.

The phase source corrects a mass omission in both paths, but only the new identical-time path can be checked as the executed physical split evolution. This distinction prevents a mathematically closed source step from certifying the old accelerated campaign.

### Native physical scope is bare connected ice only

Positive melt beneath an impermeable raised debris layer cannot simply be placed into overlying water without a drainage/storage model. Porous deposits also need a justified pore-volume/transfer closure. Accordingly `step_physical(cool=True)` refuses active covered ice **before advancing flow or removing heat**. A defensive source-level check also rejects positive covered phase transfer. Solid and porous refusal tests check that no half-applied state remains.

The old quasi-steady path retains its effective under-debris source transfer only as an explicit assumption; `covered_ice_melt_m3` and `under_debris_source_transport` expose it. No covered source is called physically verified. The `uniform` incision mode remains a nonlocal counterfactual in the legacy path and is refused as a native physical step.

## Ledgers and API

`source_log` and `last_phase_source` record the heat epoch, interval, ice and water volumes, latent heat, sensible-heat change, density-jump surface drop, source residuals, toe/bypass contributions and route-local source/throughflow diagnostics. A zero route flow gives an undefined ratio (`None`), not a zero or hidden infinity. The section diagnostic is not proof of uniform route discharge or negligible feedback everywhere.

`hydraulic_window_log` records inlet/outlet water and heat, the actual control-volume storage change, the ice heat sink and the numerical positivity correction. `combined_ledgers()` accounts separately for phase water, boundary flow, explicit state resets, and numerical added water. Its physical mass residual does **not** subtract away positivity-created mass; the correction is visible, rather than treated as a physical inflow.

Before any evolution, `reset_state` defines a new initial condition. After evolution, it records explicit mass/heat replacement and disables the physical-history certificate; it never erases completed history. A source can consume a mean-flux epoch only once. Resetting or modifying/removing deposits invalidates the available heat window. Direct external mutation of the inherited public arrays/forcing attributes remains the caller's responsibility; the existing v4 API is mutable and is not a full immutable restart system.

The original constructor and region APIs are retained: `region_heat()`, `region_vol`, `last_reach`, and the three-volume return from `lower_bed`. Diagnostics identify when their last hydraulic window is stale after geometry/source changes. At zero external inflow, normalized toe fraction is undefined rather than forced to zero. The new constructor rejects negative inflow/temperature, nonpositive essential coefficients, and unsupported CFL.

### Import and execute

Place **both** `swe2d_v5.py` and `swe2d_v4_pinned.py` alongside the existing solver. Relative and standalone imports are supported:

```python
from ice_cliff_survival.reactivation.swe2d_v5 import Sim
s = Sim(dx=0.25, Lpool=0.5, Lr=0.5, Ltail=0.25,
        Wl=0.25, Q=0.01, Tw=0.5, ramp_ratio=0)
for _ in range(10):
    result = s.step_physical(0.2)
print(s.combined_ledgers())
```

Requirements match the existing kernel: NumPy and Numba. The verification environment has NumPy2.0.2 and Numba0.60.0.

## Verification results

`VERIFICATION.json` contains the full unit-test log and actual-kernel closed/open results. Nineteen checks cover closed-cell phase mass and latent heat, density-jump surface motion, unchanged source-step momentum/enthalpy, ghost exclusion, duplicate-source rejection, explicit reset accounting, pseudo-time labelling, zero heat, repeated physical steps, invalid intervals/inputs, old region/counterfactual API, open boundary export, covered-source refusals, stale geometry, and the native scope guard.

In the closed actual-kernel example, mass residual is about `5.6e-14 kg` and thermal/latent residual about `5.5e-10 J`. Ten open-boundary physical steps give mass residual about `2.0e-13 kg` and thermal/latent residual about `9e-10 J`. These are small numerical examples, not environmental error bars.

`SHORT_TIME_REFINEMENT.json` retains a20-second numerical sensitivity with deliberately elevated prescribed heat transfer: `dt=1,0.5,0.25,0.125s` on the same coarse spatial mesh. Melt differences relative to the finest are approximately `1.30e-5,5.63e-6,2.48e-6,0` as fractions. This demonstrates a bounded splitting sensitivity in that case only. It does not establish spatial convergence, contact-event accuracy, environmental applicability or a universal order from four numbers.

## What remains unsupported

No full v4 campaign is relabelled v5. No field calibration or held-out field comparison has been rerun. No improvement over Gnielinski, constant-rate, other professional models or field observations is claimed. Long-duration contact-loss/recovery and mesh convergence are not established for v5. The legacy porous drag/insulation closures remain unverified. Cold ice, refreezing, mobile grains, under-debris drainage/pore storage, atmospheric exchange, mechanical energy, collapse and glacier mechanics remain outside this repair. The source-step tests cannot lift those limits.
