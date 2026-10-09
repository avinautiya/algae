# Claims and limitations audit (paper-ready; each claim points to an experiment and result file)

## Claims supported by the evidence

| # | Claim | Evidence |
|---|---|---|
| C1 | Representing glacier algae (any of the tested optics) improves plot-scale broadband albedo prediction over a no-algae model with identical ice/dust treatment, at S6 2017 (previously seen sites; development-reused data, not an untouched campaign). MAE improvement 0.031 (95 % day-block interval 0.018–0.042, 9 blocks). | `records/h1_albedo/h1_contrasts.csv`, P6-H1-1 |
| C2 | At PROMICE KAN_L/KAN_M, satellite-retrieved albedo with Tier A optics predicts overpass absorbed SW with MAE 33.3 W m⁻² (bias +2.1) on 126 station-days / 11 station-year blocks. Every other model is worse. | `records/h2h4/h2h4_summary.csv`, P6-H2H4-1 |
| C3 | Algae vs no algae at stations: absorbed-SW MAE improves by 8.8 W m⁻² (3.3–14.7). This is below the pre-set 10 W m⁻² minimum meaningful improvement. | `records/h2h4/h2h4_contrasts.csv` |
| C4 | An uncertainty-aware, provenance-complete satellite-interpretation and coupling software stack exists, with round-trip-validated outputs and double-counting refusals. Software claim only. | `phase4/tests/test_scene_product.py`, `model_integration/tests` |

## Claims NOT supported (state them as negative results)

| # | Statement | Evidence |
|---|---|---|
| N1 | Pigment/cell-informed (TD-DFT) optics improve independent albedo over established empirical optics: **no** (H1 −0.0019 MAE). | P6-H1-1 |
| N2 | They improve station absorbed SW: **no** (H2 M1 − M3 = −2.3 W m⁻², interval spans 0; M3 biased bright). | P6-H2H4-1 |
| N3 | The algae interpretation adds albedo information beyond a simple four-band empirical conversion: **no** (H4 +0.004, spans 0). The conversion is a stated fallback, not a published method. | P6-H2H4-1 |
| N4 | Improved melt or ablation prediction: **untested** (H3 prerequisites unmet / rule-dependent). | P5-H3-1, P6-SW-1 |
| N5 | Physics forward reflectance matches field HCRF: **no**. Four-band RMSE 0.20–0.35 is worse than a training-mean baseline; the visible is too bright by 0.2–0.3. Tested geometry adjustments do not explain the full reference-model bias. | P6-FWD-1, P7-FWD-1 |
| N6 | Station albedo from biology under the radiometer at new sites: **insufficient evidence** (2 test station-days). | P6-H5-1 |
| N7 | Quantum-chemical absorption changes glacier albedo/absorbed SW by a practically important amount after calibration: **unresolved** (corrected from "no", P8-INF-1). Plug-in contrasts are ≤ 0.008 albedo / ≤ 6 W m⁻². With calibration draws propagated, the B3LYP-D contrasts span −0.020 to +0.063 albedo, with P(\|Δ\| ≥ 0.01) up to about 0.5 at high loading. Measured-MAC uncertainty is omitted in the 24-draw version. | `records/molecular_contribution/uncertainty/posterior_contrast_summary.csv` |
| N8 | The forward-model visible/NIR deficit is explained by the tested pigment spectra: **no**. About 0.44 (day-block 95 % 0.16–0.67; denominator stable) of the observed broadband albedo–abundance slope is unexplained abundance-associated darkening relative to the specified reference model. It is largest at B8 (0.92), where the modelled pigments barely absorb. This is association, not a causal attribution. Development data (S6 2017). | `records/molecular_contribution/unexplained_fraction_day_bootstrap.csv` |

## Supported mechanistic statement

- **M1.** TD-B3LYP (6-31G*, PCM) reproduces the two-band structure of the isolated glacier-algal chromophore within about 0.1 eV without calibration; CAM-B3LYP is 0.5–0.8 eV blue. Agreement is with data later used for calibration, not independent validation. (`records/molecular_contribution/raw_vs_measured.csv`)

## Interpretation rules (P8-INF-1)

- **Practical thresholds** (0.01 albedo, 10 W m⁻²) are protocol importance thresholds. They are not measurement uncertainty or detectability; plot-albedo measurement uncertainty was not quantified separately.
- **Scope of "cannot resolve":** such statements apply to the evaluated observations and assumptions (four-band S2 / broadband albedo, posterior-mean optics, stated thresholds). Lack of evidence is reported as insufficient evidence, not impossibility.
- **Data labels:** S6 2017 is development-reused data; Williamson 2020 MAC/HPLC are calibration targets; PROMBIO has 2 test station-days.

## Words and claims to avoid

"first", "novel", "agency-ready", "operational", "NASA-endorsed", "improved melt prediction", "better climate prediction", "validated biological detection" (the scene classes are model-consistent interpretations, not validated detections), species maps, and causal effect rankings.

## Standing limitations

- **Posterior-mean optics** (D1) in H1–H5: no joint calibration draws propagated there. They are propagated only in the controlled substitution experiment (24 draws per treatment), where they widen the contrasts substantially.
- **k-prior overlap with S6 plots** (D3, not model-symmetric).
- **Campaign-informed ice structure and dust prior** for H1.
- **M2 is a simple fallback.**
- **Footprint mismatch** between radiometer and pixel.
- **Emulator SZA domain 44–56°.**
- **Low-sun SW replacement:** energy ≤ 0.5 %, but it decides H3 prerequisite verdicts.
- **PROMBIO** usable for radiation at only 3 development and 2 test station-days.

## Preferred predictive baseline

For broadband albedo and absorbed shortwave: **established empirical algae optics (BioSNICAR Tier A)** with the same ice/dust treatment. The molecular (TD-DFT) work keeps mechanistic value (pigment-level absorption, Fe-complex increment) but shows no operational predictive benefit in these tests.
