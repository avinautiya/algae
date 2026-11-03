# Review of the V2 coverage failure and proposed V3 endpoint amendment

## V2 result: unscorable

The parent reports that the authorized V2 scoring attempt found a nonfinite native camera distance at the frozen endpoint2022-12-30T03:00:00. The native operator requires finite exposed-stake distances at both endpoints. Therefore **V2 has no valid held-out ablation target, no computed evaluation error and no success/failure ranking against baselines**. This is a source-coverage failure, not evidence that the candidate is accurate or inaccurate. The prediction freeze had already been published; preserve it, its predictions and the failed scoring record. Do not replace the missing target with zero, interpolate it, use the ambiguously labelled supplied ablation column or infer it from the model.

The original V1 physical no-go also remains separate and preserved. V1 failed its no-refreezing conversion gate; V2 repaired that declared physical component but then encountered an unavailable endpoint.

This review has not read raw camera outcomes. It relies on the parent's explicit failure report. The independent result auditor is held until the parent announces that a new evaluation has actually been scored.

## Is V3 scientifically defensible?

**Yes, as a disclosed coverage-only protocol amendment, not as unchanged execution of the original freeze.** An availability check can select a usable source endpoint without selecting a favorable melt magnitude. The necessary protection is that the rule is fixed before checking availability and that no finite evaluation distance, ablation difference or model error informs the choice.

The proposed rule is deterministic: examine native endpoint labels within±12h of the original endpoint; choose the finite exposed-stake distance nearest in time; break equal-distance ties toward the earlier timestamp. If none exists, remain unscorable. The scan must inspect/release only timestamps and finite/nonfinite status, not magnitude, sign, inferred melt, model error or comparative skill. It must not reject a finite endpoint because its eventual melt increment is negative or unfavorable. Any source-quality/reset gate requires independent justification and must not become a second performance-based search.

No finite held-out ablation difference was available from the failed V2 attempt, according to the parent. This supports the claim that the amendment is based on availability rather than observed performance. It does not make the original endpoint valid or retroactively restore original-protocol success.

## Conditions before the new value is opened

1. Write a new V3 directory/specification and an immutable explanation of the V2 failure. Preserve old V1/V2 files and failed scoring logs. Record the original and replacement clocks, displacement in seconds, the exact rule, source manifest/CRC, availability-scan code and all eligible timestamp/finiteness outputs.
2. Keep all inferential choices fixed: temperature-fit diffusivity, conductivity/capacity scenarios, geometry, forcing holds, initial/deep ice assumptions, numerical method, training interval, empirical gain rule and values, primary candidate, baseline family, metric and tie tolerance. The already known TRAIN target is not a new selection opportunity.
3. Regenerate **all39** predictions to precisely the same replacement endpoint, retaining the18 forward scenarios,18 matched derivative-OLS scenarios and3 simple baselines. No favorable secondary scenario may replace the primary proxy/middle/left candidate.
4. Integrate each model to the **exact native camera clock**. If the replacement lies after the old horizon, extend the original state trajectory under the same forcing rule rather than resetting/recalibrating it. If it falls inside a native thermal or meteorological interval, use the declared partial-interval integration; do not silently snap to a different process timestamp. Refuse unsupported forcing gaps consistently.
5. Keep the TRAIN gain unchanged when its dose/target support is unchanged. If code recomputes it, verify equality against the previously frozen gain. The only prediction change should arise from the declared endpoint/support adjustment, not a new fit or parameter selection.
6. Bind new normalized physical predictions, final calibrated predictions, camera/operator metadata, specification and code into new seals/hashes. Publish and independently verify the new final prediction freeze **before** releasing the selected finite endpoint's outcome magnitude.
7. Score all models with the same native917kg/m³ exposed-stake-difference operator and exact new support. Unknown camera precision/reset semantics remain unknown. Negative finite increments are reported or trigger a separately justified predeclared validity failure; they do not authorize picking another endpoint.

## What the amended result could establish

The paper may report a transparently amended, one-block conditional hindcast against its fixed baseline family. It must give the original unscorable endpoint and the availability-only amendment before presenting any eventual ranking. It must not say that the original V2 test passed or that its endpoint was selected blindly without qualification.

Availability can be weather- or instrument-dependent. Conditioning on an available endpoint does not guarantee representative missing-data performance. No statistical confidence interval, independent-frame bootstrap, multi-glacier superiority or new environmental mechanism follows from this amendment. The same fixed-geometry, assumed initial ice state, finite deep boundary, omitted wet-debris phase heat and empirical-output-gain limitations remain.

## Audit plan

After the parent announces a successful amended score, run the additive cold-result auditor against V3's actual normalized predictions, physical freeze, calibrated prediction seal, amended analysis specification and result artifact. Verify all39 models, the unchanged primary, the common observation operator, dose-ratio arithmetic, unscaled physical ledgers, cumulative errors, all scenario results and conditional shared-middle measurement-error thresholds. Report the outcome exactly, whether favorable, tied, unfavorable or still unscorable. Do not infer a missing observation from predictions.
