# Claims and limitations audit (paper-ready; each claim points to an experiment and result file)

## Claims supported by the evidence

| # | Claim | Evidence |
|---|---|---|
| C1 | Representing glacier algae (any of the tested optics) improves plot-scale broadband albedo prediction over a no-algae model with identical ice/dust treatment, at S6 2017 (previously seen sites). MAE improvement 0.031 (95 % day-block interval 0.018–0.042, 9 blocks). | `records/h1_albedo/h1_contrasts.csv`, P6-H1-1 |
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
| N5 | Physics forward reflectance matches field HCRF: **no**. Four-band RMSE 0.20–0.35 is worse than a training-mean baseline; the visible is too bright by 0.2–0.3. | P6-FWD-1 |
| N6 | Station albedo from biology under the radiometer at new sites: **insufficient evidence** (2 test station-days). | P6-H5-1 |
| N7 | Quantum-chemical absorption changes glacier albedo/absorbed SW after calibration: **no** (≤ 0.008 / ≤ 6 W m⁻² vs measured MAC or across functionals). | `records/molecular_contribution/expA_pair_differences.csv` |
| N8 | The forward-model visible/NIR deficit is a pigment-spectrum problem: **no**. About 40 % of the observed albedo–abundance slope is non-pigment darkening (B8 −0.087 per decade unexplained). | `records/molecular_contribution/albedo_slope_day_bootstrap.csv` |

## Supported mechanistic statement

- **M1.** TD-B3LYP (6-31G*, PCM) reproduces the two-band structure of the isolated glacier-algal chromophore within about 0.1 eV without calibration; CAM-B3LYP is 0.5–0.8 eV blue. Agreement is with data later used for calibration, not independent validation. (`records/molecular_contribution/raw_vs_measured.csv`)

## Words and claims to avoid

"first", "novel", "agency-ready", "operational", "NASA-endorsed", "improved melt prediction", "better climate prediction", "validated biological detection" (the scene classes are model-consistent interpretations, not validated detections), species maps, and causal effect rankings.

## Standing limitations

- **Posterior-mean optics** (D1): no joint calibration draws propagated.
- **k-prior overlap with S6 plots** (D3, not model-symmetric).
- **Campaign-informed ice structure and dust prior** for H1.
- **M2 is a simple fallback.**
- **Footprint mismatch** between radiometer and pixel.
- **Emulator SZA domain 44–56°.**
- **Low-sun SW replacement:** energy ≤ 0.5 %, but it decides H3 prerequisite verdicts.
- **PROMBIO** usable for radiation at only 3 development and 2 test station-days.

## Preferred predictive baseline

For broadband albedo and absorbed shortwave: **established empirical algae optics (BioSNICAR Tier A)** with the same ice/dust treatment. The molecular (TD-DFT) work keeps mechanistic value (pigment-level absorption, Fe-complex increment) but shows no operational predictive benefit in these tests.
