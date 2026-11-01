# Numerical uncertainty audit and gated finest-grid plan

Pinned report: `f3602e35e1595bcdbee4208a354a34944f1ae974`. This folder changes
no active solver or frozen result. No finest-grid hydrodynamic run was launched.

## Immediate blocker found while preparing refinement

The pinned `Sim.lower_bed` applies `maximum(zb-dz,pool_floor)` to **every** cell,
including non-ice cells. A zero-duration update raises initially lower outlet
tail cells even though no melt was supplied. It also reports melt volume from
the requested `dz`, rather than the actual floor-limited ice lowering. The
isolated preflight reproduces both issues at dx=0.1 without advancing hydraulics:

- 132 non-ice cells change during `lower_bed(0)`; maximum upward change1.055m.
- A deliberately cap-activating one-cell fixture reports0.00116667m³ while
  actual ice lowering is0.00016667m³.

These are operator defects, not estimates of the error in the completed
campaign. Actual campaign effects require a versioned comparison. The parent
agent is independently auditing and correcting that operator in a separate
folder. **The finest runner refuses to compute with the pinned bad operator.**
Refining a defective operator would not strengthen the environmental claim.

The parent's separate3.5d coarse paired development check completed: legacy
loss3.164769d versus ice-only, unbounded-lowering loss3.188465d, a0.023695d
shift. Legacy declared-minus-actual melt was0.016616m³, versus zero for that
corrected operator. These results are in `../bed_audit/paired_summary.json`.
They measure one operator sensitivity, not medium/fine convergence, and do not
transfer the old GCI to the repaired model. A new versioned grid sequence is
still needed. No finest-grid run was launched while Claude fixes the operator.

## What ±0.37 days means

Using the reported rounded event times3.165,2.721,2.478d at h=.1,.05,.025m:

| Quantity | Reproduction |
|---|---:|
| Observed order p |0.869603 |
| Richardson extrapolation |2.184224d |
| Fine-grid estimated error |0.293776d |
| GCI with safety factor1.25 |0.367220d |
| Conditional next-halving prediction |2.345007d |

This reproduces the report to its precision. It is a **conditional numerical
error estimate**, not a95% confidence interval, not observation uncertainty,
and not uncertainty from heat transfer, porous flow or real boulder geometry.
The estimate assumes the same well-defined observable and an asymptotic
single-power error expansion. Three monotone points cannot establish that
assumption. In this code, divider physical thickness changes with dx, so this
sequence also changes a geometric approximation; do not call it an isolated
PDE-discretization study without explaining the zero-thickness-wall limit.

The report's times are rounded to.001d. Perturbing each by±.0005d at all eight
corners gives p≈.8604–.8788 and GCI≈.3606–.3740d. This is a rounding-sensitivity
check, not a probability interval. Recover unrounded source data before quoting
more digits. `records/rounded_report_gci.json` contains the central calculation.

`uncertainty.py` refuses a GCI estimate for censored, zero-difference,
oscillatory or noncontracting sequences. It must not treat an event that never
occurs as time zero, infinity, or the final simulation time in this calculation.

## Minimal useful finer run after the operator is repaired and refrozen

1. **Early-loss check:** one solid-strip/no-removal run at dx=.0125 for3.5d,
   6h morphology steps. It checks the early contact-loss time and pre/post-loss
   flow curve. It cannot validate day12 removal or21d integrated melt.
2. **Actual recovery boundary:** choose an already frozen case near the
   recovery/nonrecovery boundary. Run from its original initial condition at
   dx=.0125 until at least the removal time plus the declared persistence
   window. A day12 removal with1d persistence therefore needs≥13d; reproducing
   the reported nine-day absence of recovery needs all21d.
3. Do not interpolate a coarse day12 state into the fine grid and call it an
   independent refinement of the complete history. A fine-grid restart must
   originate from the same fine-grid model, forcing and geometry.

Measured medium solid-strip/no-removal cost:2028s for21d. Assuming cubic
grid-cost scaling and approximately linear horizon cost gives:

| Finest horizon | Extrapolated one-worker wall cost |
|---|---:|
|3.5d |2704s ≈45min |
|13d |10043s ≈2.79h |
|21d |16224s ≈4.51h |

These are scheduling estimates, not completion promises. Startup, morphology,
quasi-steady convergence, toolchain and machine contention can invalidate them.
Running one worker avoids consuming Claude's compute allocation. This packet
does not launch that worker automatically.

## Executable plan and safety behavior

    python3 -m unittest -v test_uncertainty.py test_finest_check.py
    python3 uncertainty.py
    python3 finest_check.py --output /tmp/m3_finest_plan

The default command only produces a plan and fast operator preflight. Adding
`--execute` still refuses the pinned faulty source before hydraulics. A future
corrected/refrozen snapshot can be supplied explicitly with `--snapshot`;
its own frozen hashes must match its source and its true provenance must be
recorded. The report reference pin is not authentication of arbitrary supplied
directories.

After a valid preflight, execution uses a separate child process, separate
output directory, one numerical thread, wall-time limit and worker log.
It refuses overwriting an existing result and never writes under the snapshot.
The inherited stock `cases.run` has no resumable checkpoint API, so this wrapper
does **not** claim checkpoint support or continuous simulated-time progress.
Timeout preserves plan/log metadata only. Add a verified restart protocol before
using it in a restart-prone environment; do not claim restart equivalence from
this wrapper.

## Reporting correction

Replace “95% interval ±0.37d” wherever present with “conditional three-grid GCI
estimate≈0.37d, safety factor1.25; asymptotic regime not independently verified.”
State that it applies to the early-loss observable in the tested geometry,
not every removal event, all21d melt totals or the physical system. Following
the bed-update repair, rerun the convergence sequence under a new version before
carrying this number into the final comparison.
