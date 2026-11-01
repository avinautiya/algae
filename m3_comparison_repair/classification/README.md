# Corrected recovery labels and shared-history accounting

Pinned source: `f3602e35e1595bcdbee4208a354a34944f1ae974`. This folder is additive;
it edits no snapshot, active Claude file, process, checkpoint or frozen result.

## Executable correction

`recovery.py` separates maintained contact, genuine post-loss sustained
reconnection, finite-horizon failed recovery, initial dryness and unresolved
states. It also identifies known pre-intervention rewetting rather than
attributing it to a later removal. It preserves onset/loss observation brackets
and requires the declared persistence duration after the actual assessment
start. A terminal wet sample never automatically establishes one-day recovery.

`classify_run()` reads **actual** clearing time from `export[0].t` in seconds;
requested `removal_day` is reported separately. Requested clearing with no
actual export record is refused. A future flood study must pass its ordinary
flow return time as `assessment_start_day`: the present constant-inflow campaign
uses prescribed removal, not a flood hydrograph.

The frozen rules used here are theta=1e-3 and persistence=1 day. Interpolation is
piecewise linear and explicitly reported; persistence is confirmed only when
the conservative observed onset/end brackets leave at least that duration.
The bracket is temporal sampling uncertainty, not a confidence interval.
Sampled observations cannot exclude unobserved sub-step dry interruptions.

Run from this folder:

    python3 -m unittest -v test_recovery.py
    python3 recompute.py
    python3 recovery.py ACTUAL_RUN.json --output classified.json

Fourteen tests pass, including terminal-wet false positives, persistence after
actual removal, truncation, missing actual export records and causal grouping.

## A real bug in the pinned metric

`snapshot/metrics.py` uses:

    duration >= SUSTAIN_D or terminal_active

Therefore ANY final active wet interval is labelled sustained, even if it lasts
less than the declared one day. It also clamps an earlier episode start to the
removal timestamp without rechecking the post-removal duration. The new code
removes both shortcuts. This proves an implementation defect; it does not prove
the reported campaign's conclusions are changed by it. Raw records are needed
to determine which completed cases, if any, are affected.

## Correct CFD class accounting — REPORT-DERIVED, not raw reanalysis

The pinned Git tree contains no raw `results/eval/E*.json`. Consequently the
following table is transcribed from the committed report's event summaries.
It is not presented as a fresh independent calculation of CFD time series.

| Case | Report-supported operational class |
|---|---|
| E01 bare, no removal | Maintained contact |
| E02 corridor, no removal | Maintained contact |
| E03 corridor, day12 removal | Maintained contact; no recovery event |
| E04 strip, no removal | Lost contact; no recovery intervention |
| E05 strip, day1 removal | Maintained contact; no recovery event |
| E06 strip, day8 removal | Genuine sustained reconnection following loss |
| E07 strip, day12 removal | Failed recovery by observed day21 |
| E08 porous, no removal | Lost contact; no recovery intervention |
| E09 porous, day1 removal | Maintained contact; no recovery event |
| E10 porous, day12 removal | Failed recovery by observed day21 |

Among all SIX removal cases: **3 maintained, 1 genuine reconnection, 2 failed
recoveries**. Among the FOUR nonporous removal cases with reduced forecasts:
**2 maintained, 1 genuine reconnection, 1 failed recovery**. Thus four matched
binary post-removal outcomes are not four demonstrated reconnection events.
State the three-class outcomes, not a generic "4/4 recoveries" statistic.

## Repeated loss times are shared trajectories

`history_key()` hashes the physical initial state, constitutive parameters,
boundary data, forcing prefix and earlier interventions. It deliberately
excludes names, future clearing and grid resolution. The caller must pass
actual causal inputs; a hash is neither proof of statistical independence nor
proof that a solver reproduced those inputs correctly.

From the frozen physical definitions and report loss times:

| Distinct pre-loss history | Cases sharing it |
|---|---|
| Solid strip, no intervention before its loss | E04, E06, E07 |
| Porous strip, no intervention before its loss | E08, E10 |

The **five reported loss cases contain two distinct pre-loss physical-input
histories**. Connectivity-aware loss-time agreement in E04/E06/E07 is **one**
solid-strip history comparison. Later removal outcomes remain different
intervention tests; grouping loss history must not erase that distinction.
These are deterministic grouping counts, not an effective sample-size estimate.

## Additional frozen-series reproducibility problem

`freeze.py` stores reduced time series using `[::4]`. Non-bare frozen series end
at **day20**, although the declared horizon is **day21**; bare ends at day21.
Downsampling dropped the final endpoint. Never silently extrapolate them for
full-horizon maintained-contact or non-recovery conclusions.

The stored full-resolution scalar loss prediction for the solid-strip model
is 2.49420556 days. Interpolating the frozen DAILY samples instead gives
2.94424936 days, bracketed by days2 and3. That discrepancy is sampling, not a
new model error. `corrected_accounting.json` keeps the full-resolution frozen
scalar and the sampled reclassification separate. The sampled series flags
unconfirmed full-horizon labels as truncated; E06 already contains sufficient
post-removal persistence before day20.

## Minimal integration with Claude's current analysis

1. Reclassify raw completed CFD records with actual export timestamps; retain
   original classifications as versioned legacy output.
2. Report maintained, regained and failed contact separately. Preserve temporal
   censoring and actual onset brackets. Inspect pre/post clearing snapshots if
   attributing the jump specifically to clearing.
3. Export full-resolution reduced predictions, or include every event time AND
   the final endpoint. Do not overwrite the frozen predictions; create a
   versioned diagnostic export with provenance and declare the sample repair.
4. Use one shared-history loss-time comparison for E04/E06/E07, while retaining
   their different intervention outcomes. Future boundary/geometry tests supply
   genuinely new predictive evidence.
5. Publish raw evaluation records or provide them in the archive so that the
   corrected metrics can be independently reproduced.
