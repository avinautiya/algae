# JSON first-contact forecasting tool

This standalone interface wraps the **unmodified pinned research kernel** from
`ice_contact_recovery/math/recovery_budget.py`. Its exact source and SHA256 are
recorded in `PROVENANCE.json`. No active Claude file was changed.

The model is `A H' = I - C max(H-z_b,0)^p`, before first contact at a fixed toe
elevation. It determines whether **some freely selected inflow** bounded by the
supplied peak, total water and deadline can reach that elevation. It does NOT
predict success of an arbitrary hydrograph with those same budgets. It does
NOT predict sustained ordinary-flow recovery, finite positive toe delivery,
cliff lifetime, or validated real-world glacier melt.

Run from this folder:

    python3 forecast.py sample_input.json --output sample_output.json
    python3 -m unittest -v test_forecast.py

Python 3 and SciPy are required. The CLI returns exit 0 for valid model output,
including an excluded case. Invalid inputs or failed quadrature return exit 2,
with a structured JSON error on stderr and no stdout prediction. JSON never
uses nonstandard NaN/Infinity constants: infinite requirements appear as null
alongside explicit failure classifications. Input duplicate keys are rejected.
Output writes are atomic; the input cannot be overwritten as the output.

## Schema version 1

Required document fields are `schema_version` (integer 1), `evidence_label`
(nonempty string), `parameters`, and `budgets`. Optional fields are `notes`
(string) and `parameter_bounds`. Unknown fields are rejected.

Parameters:

| Field | Unit / range |
|---|---|
| `area_m2` | Positive reservoir plan area, m² |
| `coefficient` | Positive bypass coefficient, m³ s⁻¹ m⁻ᵖ |
| `exponent` | Positive dimensionless power p |
| `bypass_m` | Bypass crest elevation, m |
| `toe_m` | Fixed first-contact elevation, m |
| `initial_stage_m` | Known initial water level, m, same datum |
| `peak_m3_s` | Positive upper bound on TOTAL inflow, m³/s |

Budgets are finite, nonnegative `available_volume_m3` and `available_time_s`.
Volume is **total inflow**, not excess over a baseline. Initial stored water is
already represented by initial stage. Account for every admitted additional
source; this tool does not add meltwater or storage release for you.

`parameter_bounds` maps any subset of parameter names to `[lower, upper]`.
Both endpoints obey the parameter range and contain the nominal value. The
tool evaluates all distinct deterministic corners (at most 128 for seven
parameters). These are **sensitivity samples**, not confidence intervals,
probabilities, or guaranteed extrema over the interior. Correlations and
physical feasibility of cross-combinations must be justified independently.
Varying the exponent changes coefficient units: do not combine such corners
without a consistently specified physical discharge family.

## Classification and provenance

- `initially_in_contact`: initial level already reaches the target; no new
  water is needed for this zero-depth criterion.
- `excluded_static_model`: peak, deadline, or total volume fails a necessary
  condition for this static model.
- `feasible_static_model`: there exists an admissible selectable inflow
  reaching first contact. This is neither a probability nor sustained recovery.

Numerical margins are reported but no safety margin is inferred. Near equality
and the logarithmically slow activation threshold, assess parameter and
quadrature error before interpreting a categorical result. Integration warnings
are fatal rather than silently accepted.

Each result hashes the canonical input, wrapper and vendored kernel. The pinned
kernel hash is checked before returning a forecast. The sample output
is synthetic and reproducible. Data provenance, calibration, valid model
regimes and independent CFD/field validation remain the caller's responsibility.
