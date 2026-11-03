# Version2 fixed-parameter phase-model sensitivity

No camera observations, calibration artifact or evaluation errors were imported.
Selected diffusivity, capacity, geometry, initial profile, native forcing and
calendar endpoints are unchanged. Script:`model/phase_refinement.py`; record:
`model/frozen_piramide_v2/PHASE_NUMERICAL_REFINEMENT.json`. Three helper tests pass.

The actual published baseline is32debris/64ice cells,300s steps,2m ice. An
earlier task description's24ice-cell value was corrected against the frozen
source before running. Refinements use64/128cells at150s and75s; the deep-boundary
case uses64/256cells,4m ice,150s, preserving the refined ice-cell width.

| Method, proxy capacity | Original TRAIN melt mmwe | Original EVALUATION melt mmwe | Finest ratio change | Ice depth2→4m ratio change |
|---|---:|---:|---:|---:|
| Forward alpha |362.51798|264.66881|+0.005072%|−0.0002150%|
| Derivative OLS alpha |99.07683|74.89166|+0.0007788%|0%|

The ratio is EVALUATION/TRAIN raw mass. Its relative change is also the relative
change in a prediction when reapplying the same one-TRAIN-target calibration
rule; no target value is required to compute it. The report separately gives
the multiplier obtained when holding the original empirical gain itself fixed.
These differences are finite numerical/deep-boundary sensitivities, not rigorous
error bounds, confidence intervals or evidence of held-out field accuracy.

All complementarity products are zero and substep melt masses remain nonnegative.
Relative sensible-plus-latent ledger residuals are at most6.49e−13. The primary
refined interface reaches about−0.877C during cooling. Signed deep-boundary heat
is approximately−5354J/m² for the75s/2m case and−201J/m² for the150s/4m case;
despite this change, the melt-ratio depth sensitivity is tiny. The derivative
case stays at the melting interface and has zero deep-boundary heat.

## Fixed-geometry diagnostic

The original primary raw28-day mass implies0.68396m ice-equivalent retreat,
1.368times the0.5m debris computational depth,0.342times the2m ice domain and
34.2times the initial20mm sensor-to-ice gap. The derivative raw component implies
0.18971m retreat. No such motion was applied to the solver geometry. These
ratios flag the need to explain the chosen material frame and omitted ice
advection; they do not prove that debris thickness or the sensor gap actually
changed by the full retreat. Debris/sensors may move with the melting surface.

Any later empirical output multiplier changes an algebraic ice-equivalent
retreat diagnostic, **not** the solved geometry or its conserved energy ledger.
The optional script argument accepts an explicitly supplied multiplier without
opening a calibration artifact. No multiplier was supplied for this report;
the actual calibrated-retreat diagnostic remains pending the parent's separate
training-only calibration stage. Wet-debris pore phase heat and moving geometry
are still absent from version2.
