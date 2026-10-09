# Evidence table (2026-10-09; update with every result)

Legend: ✔ yes · ◐ partial or weak · ✘ no · — not applicable.

| Component | Implemented | Numerically verified | Empirically validated | Transferred | Integrated into a named host | Status / evidence |
|---|---|---|---|---|---|---|
| TD-DFT pigment spectrum (Phase 1) | ✔ | ◐ | ✘ | — | — | Level 2 CAM/B3LYP TDA jobs still running; convergence checks in `phase1/` |
| Spectral calibration to in vivo MAC | ✔ | ✔ (quadrature refinement converged) | ◐ | — | — | AR(1) vs iid structural difference is large; ρ adequacy UNRESOLVED |
| Cell/packaging optics | ✔ | ✔ (unit tests) | ✘ (packaging and scattering not tested separately) | — | — | — |
| Emulator / surrogate | ✔ | ◐ | — | — | — | Direct-BioSNICAR benchmark running; qualification gate refuses unbenchmarked use |
| Abundance retrieval (plot scale) | ✔ | ✔ | ◐ | ◐ (S6 → S Greenland) | — | `records/heldout_v2`: log-score gain only; RMSE worse; no gain over prior on reverse fold (P5-CORR-1, D1–D3) |
| Forward plot reflectance | ✔ | ✔ | ✘ | ✘ | — | Physics over-predicts broadband HCRF; climatology wins. Diagnosis: `docs/forward_diagnostics.md` |
| **H1** independent albedo (pigment/cell vs established optics) | ✔ (`phase4/h1_albedo.py`) | ✔ | pending | ✘ | — | Queued under budget; decision rule frozen |
| **H2** absorbed SW at stations | ✘ | — | ✘ | ✘ | — | Not yet tested: needs KAN_L/KAN_M hourly radiation and station pixels |
| **H3** observed ablation | ◐ (SEB) | ✔ (closure, convergence, gaps) | ✘ | ✘ | reference SEB only | Prerequisites met only at KAN_L 2016 and 2022, and 2016 has no S2 L2A, so ≤ 1 block; **UNTESTED** (P5-H3-1) |
| **H4** vs empirical satellite albedo | ✘ | — | ✘ | ✘ | — | Not yet tested |
| Reference point SEB (`phase4/seb.py`) | ✔ | ✔ | ◐ (1.5–1.7× over-prediction) | ✘ | is the host | χ = 0.3 is fitted, not validated |
| Host coupling Modes A/B | ✔ (`model_integration/host_modes.py`) | ✔ (unit tests) | ✘ | — | reference SEB (offline) | No E3SM/MAR integration |
| Scene interpretation product | ✔ (`phase4/scene_product.py`) | ✔ (synthetic tests, round trip) | ✘ | ✘ | — | Research-grade; pixel scale not validated |
| Uncertainty: joint calibration draws | ✘ | — | — | — | — | Single `posterior_mean` draw (D1); follow-up |
| Compute budget | ✔ | ✔ (fixture tests) | — | — | — | `docs/compute_reliability.md` |

No claim of "first", "agency-ready", "better climate prediction" or "improved glacier melt" is supported by this table.
