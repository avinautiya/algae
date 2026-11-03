# A separate source-bounded neutral-profile energy hypothesis

This is a new **unscored** hypothesis, not a revision of the failed four-timescale radiation-reservoir family. T1 and T2 camera targets are already inspected development data. No Tower3 or later Tower1 camera outcomes were opened. No new target comparison is reported here.

## One fixed proposed closure

Use the already-source-bound weather CSV incoming shortwave and longwave, with albedo 0.20 and emissivity 0.94 as explicit inherited assumptions, not measured site radiometric properties. The subsequently audited author MAIN_FRAME.m lines 434–436 establish that SAB1+SAB2+SAD1+SAD2 is total incoming shortwave, with albedo applied afterward; it must never be treated as already absorbed radiation. This specification retains the already-bound weather CSV SWIN rather than selecting among sources by target fit. Let the buried upper thermistor be a **surface-temperature proxy**, not a radiometric surface observation. Define

`Q* = 0.80 SW_in + 0.94 [LW_in − sigma (T_buried + 273.15)^4] − H_up − LE_up`.

All terms have units W m^-2; the neutral-profile feature exports define H and LE positive **upward**, so both are subtracted. Site-specific source MATLAB pressure in Pa is used in the source-agent density/humidity calculation. This avoids inventing sea-level pressure or confusing hPa with Pa. Same naive timestamp labels are bound by the author's weather/driver code and the matching-pressure signal; physical UTC remains unresolved. An extra midnight row in the MATLAB export is a coverage difference, not a demonstrated one-hour phase shift.

The profile feature closure is neutral, zero displacement, horizontally homogeneous and constant-flux, using four source instrument heights 0.5/1/2/2.7 m. Those assumptions are not validated by small Richardson number or high regression R². Vertical profiles do not supply an independent surface heat-flux measurement.

Use only the fixed tau=0 effective cold-deficit-repayment limit of `radiation_reservoir.step`, initialized by the same temperature-profile energy proxy. Negative net supply is retained as a deficit; it is not clipped into extra melt. Positive supply first repays that deficit, then becomes nonnegative deliverable energy G. Do not call this deficit measured cold-ice content. One local TRAIN gain A, constrained to [0,1], would give `m_mmwe = A sum(G)/334000`. The remaining `(1−A)G` is an explicit unresolved partition, not independently measured atmospheric/storage loss. No fitted time scale, stability parameter or extra target-calibrated radiation multiplier is allowed. This remains a conditional index, not a closed full surface-energy-balance model.

If data support the closure, score all four development periods with the old air and calendar rules and preserve all errors. No fresh validation or novelty claim would follow from these inspected periods.

## Actual input-support gate: NO-GO

`source_gate.py` evaluated native process feature coverage only, before any new predictions or target-error calculation. A supported time span must have two admissible native bracketing profiles, a gap no larger than 30 minutes, finite fluxes/pressure and positive fitted wind log-slope. Invalid rows remain invalid; they are not zeroed or bridged. The 30-minute criterion is explicit native-support accounting, not a physics theorem.

The four evaluation-period nominal support fractions are:

| Tower | Period 1 | Period 2 | Period 3 | Period 4 |
|---|---:|---:|---:|---:|
| T1, nominal computable neutral profiles | 92.0% | 93.9% | 97.3% | 97.4% |
| T2, nominal computable neutral profiles | 84.8% | 85.0% | 92.2% | 94.8% |
| T1, additionally abs(Ri) <= 0.03 | 30.3% | 32.1% | 35.8% | 33.8% |
| T2, additionally abs(Ri) <= 0.03 | 35.9% | 25.0% | 39.6% | 39.4% |

TRAIN nominal support is 96.0% and 90.2%, respectively. TRAIN near-neutral support is 24.5% and 32.7%. Requiring both temperature and wind gradient signs to be resolved under the supplied instrument-accuracy assumptions leaves only 4.6–13.7% of T1 evaluation time and 5.4–18.8% of T2 time. These tiers are diagnostics, not probabilities, confidence intervals or proofs of flux accuracy. The gradient-error bounds are conditional source-agent calculations; humidity uncertainty is not exhaustively resolved by these two sign tests.

The source features contain four explicit nonpositive wind-log-slope refusals at each tower, including TRAIN, and native gaps as long as 16.83/16.17 h. Therefore neither tower supports a complete, native full-period turbulent-flux energy comparison. Missing fluxes may be consequential precisely during non-neutral conditions. There is no supplied defensible flux bound that would make their omitted integral negligible. Normalizing only the surviving samples would change the estimand; zeroing them would invent measured zero flux. Neither is done.

No melt predictions or errors were calculated for this unscorable hypothesis. `NATIVE_SUPPORT_GATE.json` preserves all period/tier support fractions and source feature hashes. Five synthetic tests verify no gap bridging, no invalid-to-zero conversion and exact partial-span support. This gate does not establish that the true environmental heat budget is small or that an energy-balance model cannot work: it establishes that these particular neutral-profile exports do not independently close it.

## What would make a subsequent test legitimate

A separately frozen model would need a justified atmospheric exchange/stability closure valid over the missing regimes, actual surface-temperature and evaporation/wetness constraints, and explicit uncertainty or reconstruction rules. A fresh site/event with these inputs and independent ablation would then be needed for field performance. A generic stability function, storage term or conductivity inversion would also need contrast with established debris surface-energy models before any novelty claim. Do not grow another family merely to rescue these inspected camera errors.
