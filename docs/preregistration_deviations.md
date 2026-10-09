# Deviations between the held-out pre-registration and its implementation

The pre-registration (`docs/preregistration_heldout.md`, frozen before any score existed) is not edited.
Deviations found afterwards are recorded here, and the frozen results are interpreted with them in view.

## D1 — molecular uncertainty was NOT integrated over joint calibration draws (P5-PREREG-1)

- **Pre-registered text:** "Optics: molecular parameters from joint posterior draws".
- **Implementation:**
  - `phase4/emulator.py:211–212` builds the pigment MAC as `cal.mac_D` (tier D) or `cal.mac_C` (tier C).
  - `Calibration.mac_D/mac_C` (`phase2/tddft_calibration.py`), called without arguments, evaluate the posterior MEAN of (ΔE, w, f, φ). That is one fixed optical model.
  - The joint draws exist only in Phase 3 (`cal_draw`), which the held-out experiment does not use.
- **Consequence:**
  - The held-out predictive distributions of `tddft_C`, `tddft_D` and `tddft_D_iid` contain no molecular-calibration uncertainty. Their spread comes from the retrieval posterior (abundance, community fraction, radius, dust, k) and from the fitted discrepancy τ.
  - The fitted τ (0.72 dex for tddft_D on the primary fold) absorbs whatever the omitted calibration uncertainty would have contributed, together with all other structural error. It cannot be read as calibration uncertainty.
  - The frozen results (`records/heldout_v2/`) must be described as **"posterior-mean optics plus fitted discrepancy"**, not "integrated over joint molecular uncertainty".
- **Not done retrospectively:** the frozen experiment is not re-run under a changed specification. A follow-up that carries joint draws through the whole chain is specified in `docs/glacier_model_validation_protocol.md` (§ Uncertainty).

## D2 — the primary-fold block interval is degenerate

- The rule asked for a day-block bootstrap interval. The primary test site has 2 sampling days, so the bootstrap can only return the two day means; the reported interval (0.117–0.168) is their range.
- The pre-registration anticipated the interval being too wide, not degenerate.
- The frozen result is therefore reported as: "direction consistent; generalisation uncertainty at the primary site not estimable".

## D3: the HCRF/albedo anisotropy prior k was derived partly from scored plots (found 2026-10-09, after scoring)
- **What the code does.** `empirical_data.anisotropy_prior()` builds the k prior N(0.90, 0.175) from all 51 spectra in biosnicar-py's `ARF_master.csv` (`data/empirical/biosnicar_field_ARF.csv`).
  - Those ARFs are exactly HCRF / hemispherical albedo of S6 2017 plots. Checked numerically: ARF = HCRF(`cook2020_archive_hcrf.csv`) / albedo(`Albedo_master.csv`) to 4 decimals for shared plots.
  - 29 of the 51 are counted plots in the `s6_2017` held-out set.
- **Consequence for `records/heldout_v2`:**
  - In the reverse fold (test = S6 2017), the prior on the nuisance k includes the anisotropy of 29 of the 47 scored plots. This is pooled (one mean and SD over 51 spectra) and applies to all models equally, so it is a mild, model-symmetric leakage.
  - It does not favour any optics variant. It can make every physics model's reverse-fold score slightly optimistic relative to a true transfer.
  - The primary fold (test = S Greenland 2021) is unaffected.
- **Consequence for H1.** H1 predicts albedo directly and does not use k, so H1 predictions are unaffected. But the statement in the H1 protocol that the albedo spectra "were never used in any fit" is wrong: they entered the k prior. See the erratum in `docs/glacier_model_validation_protocol.md`.
- **Frozen results are not re-scored.** A leakage-free k prior (leave-one-day-out ARF) is a recorded follow-up.

### D3 correction (2026-10-09, later the same day)
- **Wrong statement.** D3 called the k-prior leakage "model-symmetric ... does not favour any optics variant". That is wrong.
- **Why.** k multiplies each model's own forward reflectance, and the physics models' brightness biases differ (+0.10 to +0.26 broadband HCRF at S6). The benefit of a test-informed k prior therefore differs by model and need not cancel in a contrast.
- **Analysis and planned sensitivity:** `docs/data_use_register.md`.
