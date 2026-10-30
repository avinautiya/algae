# Limitations (milestone 1)

1. **Contribution is not unique.** Separating algae, mineral and cryoconite darkening on bare ice with a Bayesian inversion and a weathering-crust column is already published and released (Dark Zone paper; Zenodo archive). Only the conditional hypothesis about variable absorber chemistry remains, and it may yield no gain.
2. **No independent truth for mineral, cryoconite or dark-particle loadings** in any dataset in hand. Fit quality cannot establish correct attribution; absence of truth does not show an inversion is wrong.
3. **Data scope is Greenland only.** No Himalayan spectra or samples. Nothing here supports Himalayan transfer.
4. **Two Greenland campaigns are reused and ancestrally entangled.** Optics and fitted parameters derive from the same groups (Chevrollier 2023; S21 supplies both optics and retrieval targets). The "transfer check" is not an untouched test. The only untouched candidates (PROMBIO 2024, `heldout_v2`) are reserved.
5. **HCRF is not albedo.** S Greenland and Dark Zone spectra are directional; the free anisotropy factor k absorbs, but cannot validate, the difference. The Dark Zone method assumes Lambertian; so do we in the transfer check, with k free.
6. **Error model is partly assumed.** sigma_obs = 0.015 is from documentation not re-verified against the source paper; sigma_total and AR(1) inflation are estimated from development days. No model-discrepancy process is fitted.
7. **Baseline E is MAP, not NUTS**, on an emulator whose outputs we could not compare with TensorFlow. Its reported error definition is unknown, so agreement with 1.8 % is indicative only. The scaler pickle was not loaded; constants were read from raw bytes and column order assumed from the paper.
8. **Dark Zone paper read partially**: supplementary tables S1-S4 and figures S1-S8 were not read; the SNICAR-fx release in the paper (Zenodo 18609613) was not downloaded and is not shown to equal commit c1780ab.
9. **snowlaps is a snow tool** with snow-algae optics; none of its numbers transfer to glacier-ice algae.
10. **B is unspecified numerically** until two measured end-member spectra are located. C (PINN) is not justified or run.
11. **Melt, absorbed energy and albedo-prediction endpoints cannot be tested** with available data.
12. **Computation is small by design**: single-thread, budgeted; no chemistry, training or satellite download was run.
