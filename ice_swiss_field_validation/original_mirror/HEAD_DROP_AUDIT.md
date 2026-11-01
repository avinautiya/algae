# Independent head-expression audit

## Finding

A source-level sign inconsistency is confirmed in `Hydraulic_thermodynamics_parameters.ipynb`, zero-based cell131, function `friction53`. The original author-mirror notebook was retrieved by a bounded range request, verified against original ZIP CRC32 1778057482 and uncompressed size 1,922,595 bytes. Its SHA256 is `33750db4007b72418345f2ee3e93fb87165ac97eeb4a24e5d98fdd4dc149bf45` and it is byte-identical to the previously retrieved GitHub notebook at commit `9116f98dc7fef42be82dc4de50316f4a676ca8e9`. No authors' file or product was modified.

The code computes:

```julia
slope = (z5 + h5 - z3 + h3)/(dist["P3"] - dist["P5"] ± 2)
```

Cell133 reads `z3,z5` from original channel-floor elevation CSVs. Cell137 defines `h3,h5` as positive water depths. Thus the piezometric head drop is `(z5+h5)-(z3+h3)`, whereas the source expression adds downstream depth. The excess head is exactly `2*h3`. Equal positive depths should cancel from a head difference; the archived expression fails this limiting case. The paper's Eq5 defines hydraulic slope as the difference between upstream and downstream bottom-elevation-plus-stage heads, consistent with the subtraction, not this implementation. [Primary paper](https://tc.copernicus.org/articles/15/5133/2021/tc-15-5133-2021.html#section3).

This establishes an algebraic inconsistency in the audited source. It does not establish the full provenance of every number in the published article, invalidate all results, or prove a new model's field superiority.

## Nominal independent recomputation

`audit_head_drop.py` retains the seven dates in the archived exported friction CSV, original direct floor elevations, the notebook's first-same-day depth convention, discharge products and saved rounded tracer velocities. It evaluates both width-fixed and stage-fixed cross-section branches. It intentionally does not reconstruct the authors' stochastic particles; original exported Monte Carlo means differ from evaluating nonlinear expressions at nominal inputs. No stochastic bit-for-bit reproduction is claimed.

| Case | Archived-expression head m | Piezometric head m | Nominal width-fixed f, archived expression | Nominal width-fixed f, piezometric expression |
|---|---:|---:|---:|---:|
| 11 July AM | 0.95 | 0.49 | 0.4125 | 0.2127 |
| 11 July PM | 0.95 | 0.49 | 0.3444 | 0.1776 |
| 16 July | 2.26 | 1.56 | 0.1735 | 0.1198 |
| 25 July | 1.65 | 0.55 | 0.1782 | 0.0594 |
| 30 July early | 1.23 | 0.23 | 0.4946 | 0.0925 |
| 30 July late | 1.23 | 0.23 | 0.1987 | 0.0372 |
| 31 July | 1.04 | 0.28 | 0.3584 | 0.0965 |

All inputs and outputs are in `HEAD_DROP_NOMINAL_CASES.csv` and JSON. The original 16 July exported Monte Carlo mean is 0.171, close to but different from the nominal archived-expression value 0.1735. This supports the source-path interpretation without pretending exact particle reproduction. The corresponding nominal subtraction result 0.1198 is a sensitivity calculation, not a validated replacement coefficient.

Both width-fixed and stage-fixed Darcy factors scale linearly with head drop when the other same-case inputs are fixed. Manning roughness scales with the square root of that factor. The notebook's selected branch computes `ffw_53` from fixed width; the alternative computes `ffh_53` from prescribed mean stage. The same signed head expression feeds both.

## Downstream product scope

- Cell144 computes `ffw_53`, `nw_53`, `ffh_53`, and `nh_53` using `friction53`; all inherit the signed-head expression.
- Cell176 exports `Friction_factor_P5_P3.csv` and `Manning_roughness_P5_P3.csv` from the width-fixed arrays (excluding the first 10 July observation). This source dependency is explicit.
- Cells228/229 use `mean.(ffh_53)` for discrete experiment-time Gnielinski values `nu3_g`; these are affected through the stage-fixed friction branch.
- The continuous daily/hourly Gnielinski arrays in those same cells explicitly set `fr=0.3`; they do not inherit this source expression directly. Do not automatically replace or reinterpret their archived CSV values.
- Standard Dittus–Boelter predictions and direct measured temperatures/elevations do not acquire this sign error merely by appearing in the same archive.

## Limits before a scientific correction or accuracy claim

The true frictional energy loss also depends on velocity-head changes and other open-channel assumptions; substituting a piezometric difference does not resolve those approximations. The notebook uses same-day rather than simultaneous depths and floors, and the two 11 July P5 depths are deduced rather than directly observed. On 31 July downstream depth is from the pressure sensor. The first-same-day selection also uses the earlier P3 depth for both 11 July experiments; this audit preserves that behavior rather than quietly improving it.

Near the low-head late-July cases, measurement uncertainty is a material fraction of the corrected nominal head. For independent illustrative errors of 0.1 m per floor and 0.05 m per depth, root-sum-square head uncertainty is about 0.158 m; correlations are unknown, and the source endpoint uncertainties are not a universal probabilistic model. A conservative sum bound is 0.30 m, exceeding the nominal 30 July head of 0.23 m. Thus those nominal inferred friction values should not be asserted as precise field truth.

Six independent limiting/property tests pass: depth cancellation, equal-head zero, common-datum invariance, endpoint sign reversal, negative-depth rejection and friction proportionality in both cross-section branches. `plot_head_drop.py` generates PNG/PDF comparisons from saved values, without fitting. The source-level discovery is meaningful for auditing hydraulic calibration, but there is no demonstrated improvement in held-out melt or discharge prediction here. Preserve published baselines as published and version any future correction separately.
