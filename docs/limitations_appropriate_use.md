# Limitations and appropriate use

**Status:** research-grade. Nothing in this repository is agency-ready or operationally validated.

## Appropriate uses
- **Hypothesis testing and sensitivity studies** of algal optics in offline models, using `model_integration` Modes A/B with a declared host contract.
- **Screening Sentinel-2 scenes** for where algal darkening is *supported vs ambiguous*, read with the scene report and its exclusions.
- **Comparing optical representations** under identical forcing (`phase4/seb.py` paired runs). The outputs are modelled increments conditional on the SEB.

## Inappropriate uses
- **Melt, runoff or SMB attribution.** The reference SEB over-predicts observed ablation at KAN_M with measured albedo (1.5–1.7×), H3 is UNTESTED, and χ = 0.3 is a fitted scenario.
- **Species composition maps:** f_n is not identifiable from 4 bands.
- **Adding a Δα to an observed albedo:** double counting, which Mode B refuses.
- **Applying the S6-calibrated σ, τ and radius prior elsewhere** without stating it is a transfer.
- **Treating 95 % intervals as calibrated.** They are conservative at plot scale and include τ. Calibration uncertainty is not sampled (D1).
- **Using the k prior as independent information for S6 2017 plots** (D3).

## Known failure modes
- **Systematic broadband HCRF over-prediction** at plot scale (+0.1 to +0.3). The diagnosis is in `docs/forward_diagnostics.md`.
- **Training-site band climatology** predicts plot reflectance better than any physics model; the physics advantage is in abundance retrieval only, and weak.
- **Sen2Cor SCL** mislabels dark ice. Dark/unclassified pixels are kept and flagged, so check the flag share in each report.
- **Out-of-domain illumination** (no emulator SZA node) is refused, not extrapolated.
