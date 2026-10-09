# Repair ledger

Reviewed public commit: `814295ab2809069a1c892bf7331e37ebbab30878`. Each entry records: source function; the verified defect or limitation; the repair; the test; affected outputs; the before/after result; and the remaining limitation.

Statuses:
- **FIXED_AND_VERIFIED**: repaired, with a regression test passing and production-relevant behaviour verified.
- **IMPLEMENTED_AWAITING_PRODUCTION_RUN**: code repaired and tested; the production result does not exist yet.
- **DATA_LIMITATION**: cannot be resolved with the available data or resources; stated as such.
- **UNRESOLVED_DEFECT**: known defect, not yet repaired.

Test commands are given per entry. Outcomes are the actual pytest results on this container (4 CPU, 16 GB, PySCF 2.14.0).

---

## 1. Execution environment and job durability

### P1-ENV-1: container termination vs numerical failure
- **Source:** environment / `scripts/phase1_production.sh`.
- **Finding:** the processes were killed by provider/container recycling, not by memory exhaustion.
  - The memory log (`phase1/results/memlog.txt`) shows ≥ 9 GB available at every termination.
  - The CAM-B3LYP Level 2 RSS was 6.4–7.2 GB on 16 GB.
  - Terminations happened both while the session was idle (within about 5 min) and while it was active (once, after about 2.4 h).
- **Storage:** `/home/user/algae` and the scratchpad survive restarts (persistent); processes do not (ephemeral).
- **Repair:** durable job execution.
  - Job definitions live in `phase1/jobs.py`.
  - Each job writes its outputs to `phase1/results/v2/<job>/` on persistent disk.
  - Excited-state solves keep chunk-level checkpoints in `<outdir>/td_ckpt/`.
  - Job state is kept in `results/v2/jobs_state.json`, written atomically.
  - A dead PID is detected and the job relaunched.
  - Restart: `python phase1/jobs.py run`.
- **Resource statement:** a durable 4+ core host (or a GPU host verified with gpu4pyscf for PCM response) is not available in this session. GPU support is NOT assumed: it has not been verified or benchmarked here.
  - All jobs run here, and each restart costs at most one checkpoint chunk.
  - Full-TD CAM-B3LYP at Level 2 (> 1 day of Davidson work) is NOT scheduled, because it cannot finish reliably here.
  - It needs a host with ≥ 4 cores and ≥ 12 GB RAM that stays up for about 2 days. Command: `python phase1/run_phase1.py --level level2 --skip-opt --start-xyz phase1/results/level2/opt_B3LYP_final.xyz --tddft-functionals CAM-B3LYP --nstates 30 --td-conv-tol 1e-5 --outdir phase1/results/v2/L2_CAM_FULL30`.
- **Status:** FIXED_AND_VERIFIED for job durability (tests below). DATA_LIMITATION (compute resource) for full-TD CAM-B3LYP.

### P1-JOB-1: blocking queue with false completion
- **Source:** `scripts/phase1_production.sh::run`.
- **Defect:**
  - Any nonempty `summary.json` was treated as completion.
  - Failures did not stop the queue, and `QUEUE DONE` was printed regardless.
  - Independent jobs (Level 1, Fe) were blocked behind CAM-B3LYP.
- **Repair:** `phase1/jobs.py`, a dependency-aware registry and runner.
  - A job is complete only when `completion.validate_run` accepts its output.
  - Failed jobs are not retried within a run, and stop after 3 attempts unless `--retry-failed` is given.
  - Dependencies on failed jobs never run.
  - CPU threads and committed memory are budgeted, including PySCF processes the runner did not start.
  - The runner exits non-zero if any job failed.
  - The old script is superseded and kept only for history.
- **Tests:** `phase1/tests/test_jobs.py`: queue failure, a nonempty-summary-only job, dependency blocking, and relaunch of an interrupted job (2 passed).
- **Affected outputs:** all Phase 1 production runs.
- **Before/after:** before, a failed or partial run could print QUEUE DONE; after, it is recorded as `failed` with reasons and exit code 1.
- **Status:** FIXED_AND_VERIFIED.

## 2. TD/TDA checkpoint and completion logic

### P1-CKPT-1: stage checkpoint validity
- **Source:** `phase1/qc.py::_staged_tda` (removed). Replaced by `phase1/tdcheckpoint.py`.
- **Defects:**
  - Resume compatibility was checked only by `nstates`.
  - A stage was marked done without checking that all roots converged.
  - Saved MO-basis vectors were reused after an SCF restart without checking the orbitals.
  - Saves happened only at stage boundaries.
- **Repair:**
  - **Operator fingerprint:** geometry (Bohr), atoms, charge/spin, basis hash and Cartesian flag, functional, RKS/UKS, density-fitting basis, grid level/prune, solvent model/eps/Lebedev order, equilibrium flag, TD/TDA, singlet, frozen orbitals, root count, PySCF/numpy versions. The tolerance is excluded, so changing it allows reuse.
  - **Orbitals:** stored with the vectors. On resume they must be the same canonical set (orbital energies within 1e-6 Eh, mixing only among orbitals within 1e-5 Eh) with unitary occupied/virtual overlaps. Vectors are then rotated, `X' = U_ooᵀ X U_vv`; otherwise the checkpoint is invalidated.
  - **Restart vectors and stage credit kept separate:** latest restart vectors; the last stage at which *every* root converged, with its energies and oscillator strengths (`stage_history`); and the final result.
  - **Writes:** atomic (tmp + fsync + rename) with `latest.prev.npz`; a lock prevents concurrent writers; schema, shape, finiteness and root count are validated. A corrupt latest file falls back to prev; a corrupt pair is quarantined and the solve restarts.
  - **Chunked solves:** PySCF 2.14's `lr_eigh`/`real_eig` expose no callback, so the solver runs in chunks of `--td-chunk` iterations (`max_cycle` + `x0`). This is chunk-level checkpointing, not per-iteration.
- **Defect found during testing:**
  - A restart from Ritz vectors alone (the old code's behaviour) converged root 4 of a test molecule to the wrong eigenvalue: 0.3805 instead of 0.3756 Eh, with small residuals.
  - Repair: restarts use the Ritz vectors augmented with the solver's standard initial guess, orthonormalised. If the solver rejects that set, it falls back to the plain Ritz vectors, then to a fresh guess, and logs the fallback.
  - Consequence: the legacy CAM-B3LYP TDA result (whose second stage restarted from Ritz vectors only) is NOT used as a result, only as an unvalidated starting guess.
- **Tests:** `phase1/tests/test_tdcheckpoint.py` (16 passed). They cover:
  - uninterrupted vs plain energies AND oscillator strengths (TDA and RPA);
  - interruption before stage completion, then resume;
  - between stages and a changed final tolerance;
  - `max_cycle` without convergence (not success);
  - changed geometry, functional, basis or solvent with identical `nstates` (not reused);
  - orbital sign changes (transformed; energies and f equal to the reference);
  - non-canonical rotation and occupied–virtual mixing (invalidated);
  - corrupt latest (recovered from prev); truncated pair (quarantined, clean restart);
  - an already-completed restart (verification pass); a concurrent writer (refused).
- **Scope:** these verify restart mechanics on formaldehyde/PCM, not target-pigment accuracy.
- **Status:** FIXED_AND_VERIFIED.

### P1-CKPT-2: unconverged results exported as production
- **Source:** `phase1/qc.py::run_tddft`, `phase1/run_phase1.py::main`.
- **Defect:** unconverged roots only produced a warning; energies and oscillator strengths were exported and consumed downstream.
- **Repair:**
  - `run_tddft` returns `status` and per-root `converged`.
  - `completion.td_problems` checks root count, convergence, finite positive energies, finite non-negative f, and ascending order.
  - `completion.write_marker` writes `COMPLETE.json` (production or provisional_geometry) only when there are no problems; otherwise `DIAGNOSTIC.json`.
  - Artifacts are checksummed; `validate_run` rejects tampered or partial outputs.
  - `run_phase1` exit codes: 0 production, 3 provisional geometry, 4 diagnostic only.
- **Tests:** `phase1/tests/test_completion.py` (5 passed): failure detection, provisional vs production geometry, artifact tampering, unconverged run → diagnostic only, summary without marker → rejected.
- **Status:** FIXED_AND_VERIFIED.

### P1-LEGACY-1: production reference (Level 2 B3LYP full TD, 30 roots)
- **Audit:** `jobs.py::legacy_audit`.
  - The run log contains no runtime "roots not converged" warning (residual tolerance 1e-6).
  - 30 roots, 2.949–6.260 eV, all f ≥ 0.
  - Geometry not converged: only 1 of 5 criteria met.
- **Marker:** `phase1/results/level2/COMPLETE.json`, status `provisional_geometry`.
- **Remaining:** the geometry sensitivity is pending (jobs L2_OPT → L2_B3LYP_TDA15_RELAXED vs L2_B3LYP_TDA15).
- **Status:** IMPLEMENTED_AWAITING_PRODUCTION_RUN.

## 4. Geometry and parsing

### P1-GEOM-1: unconverged geometry justified by the energy change
- **Source:** `phase1/results/level2/stop_criterion.json`, `phase1/README.md`, `records/MASTER_SUMMARY.md`.
- **Defect:** a 5e-6 Eh ground-state energy change was used to justify the geometry. That does not bound excitation-energy error.
- **Repair:** the claim is withdrawn in `stop_criterion.json`, which now records criteria met/not met and `geometry_status: not_converged`. The README and summary corrections are pending (see P1-DOC-1).
- **Measurement:** jobs L2_OPT and L2_B3LYP_TDA15_RELAXED.
- **Status:** IMPLEMENTED_AWAITING_PRODUCTION_RUN.

### P1-GEOM-2: optimiser history and finality
- **Source:** `phase1/qc.py::optimize_geometry`.
- **Defect:** the trajectory was truncated on entry, and `_final.xyz` was written even when the optimisation had not converged.
- **Repair:**
  - Trajectory and per-evaluation criteria (`_steps.jsonl`: E, ΔE, gradient RMS/max in Eh/Bohr, displacement RMS/max in Å, segment) are appended.
  - `_opt_record.json` holds thresholds with units and the convergence flag.
  - `_final.xyz` is written only when converged; otherwise `_unconverged_endpoint.xyz`.
  - The job runner resumes optimisations from `_last.xyz`.
- **Test:** `test_xyz_geometry.py::test_optimizer_final_only_when_converged_and_history_appended` (passed).
- **Status:** FIXED_AND_VERIFIED.

### P1-XYZ-1: XYZ parsing
- **Source:** `phase1/qc.py::xyz_file_to_atom_block`.
- **Defect:** blank lines were removed before parsing, so a valid blank comment line shifted every frame.
- **Repair:** `read_xyz_frames` parses frame by frame (count, one comment line, atom lines), validates coordinates, and recovers the last valid frame on truncation or malformation (with a warning), or raises in strict mode.
- **Tests:** 7 tests in `test_xyz_geometry.py` (passed): blank comment, multiple frames, malformed coordinates, truncated last frame, truncated only frame, malformed after valid, strict mode.
- **Status:** FIXED_AND_VERIFIED.

### P1-COV-1: root coverage judged at the wrong wavelength
- **Source:** `phase1/run_phase1.py` (coverage flag at `--lam-min` 300 nm) vs `phase2/tddft_calibration.py` (fits 260–750 nm; normalises the HPLC shape and the Fe increment over 265–600 nm).
- **Defect:**
  - Coverage was checked at 300 nm, but the calibration consumes 260–750 nm.
  - The check itself (highest root + 2·FWHM reaches λ_min) is not a convergence test.
- **Repair:**
  - `spectra.root_count_sensitivity` and `tddft_calibration.root_count_sensitivity` measure how much the broadened MAC in the downstream window changes when the top 5 roots are removed.
  - Phase 1 records it in `summary.json` at the run FWHM and at 0.6 eV.
  - The calibration records it at the posterior-mean (ΔE, w) and warns above 1 %.
  - The old flag is kept for information only, and `--lam-min` now defaults to 250 nm.
- **Tests:**
  - `phase1/tests/test_coverage.py` (2 passed): windows agree with the calibration constants; a root inside the window is flagged although the old criterion passes it.
  - `phase2/tests/test_phase2.py::test_calibration_root_count_sensitivity_flags_truncated_window` (passed).
  - `test_completion.py` checks the record (passed).
- **Measured on the production spectrum** (`results/level2`, 30 roots, ΔE +0.06 eV, w 0.62 eV), removing the top k roots:

  | k removed | max change in 260–750 nm (rel. to max) | 265–600 nm integral | at 260 nm |
  |---|---|---|---|
  | 5 | < 1e-5 | < 1e-7 | 1e-5 |
  | 10 | 0.4 % | 6e-5 | 1.2 % |
  | 15 | 0.6 % | 1.5e-4 | 1.9 % |
  | 20 | 10 % | 2.6 % | 28 % |

- **Remaining limitation:** roots above the 30 computed are not computed. Their effect is bounded only by this trend (evidence, not proof). The 15-root TDA runs are adequate for band positions above about 300 nm, but not for the MAC at 260 nm (about 2 %).
- **Status:** FIXED_AND_VERIFIED.

### P1-DOC-1: unsupported statements in documentation
- **Source:** `phase1/README.md`, `records/MASTER_SUMMARY.md`.
- **Defects:**
  - Geometry was justified by a 5e-6 Eh energy change.
  - "TDA within about 0.1 eV", "slightly less reliable" and "1e-3 suffices" were stated as facts for this molecule.
  - "No charge-transfer error" was inferred from the small fitted shift.
  - The runs were described as executed by the superseded queue script.
- **Repair:**
  - These statements are replaced by what was measured, or marked pending with the job that measures them.
  - The unaudited chemical-model choices (carboxylate, tautomers, conformers, basis set, explicit solvent) are now listed as limitations.
- **Status:** FIXED_AND_VERIFIED (text). The measurements it points to are tracked under P1-CMP-1 and P1-GEOM-1.

### P1-CMP-1: approximation checks with state matching by overlaps
- **Source:** new `phase1/compare_states.py`.
- **Implementation:**
  - Comparisons:
    - A: B3LYP full vs TDA.
    - B: B3LYP TDA vs CAM-B3LYP TDA.
    - R: 15 vs 25 roots.
    - G: start vs relaxed geometry.
    - T: per-run tolerance stages, 1e-3 → 1e-4 → 1e-5.
  - States are matched by AO-basis transition-density overlaps with Hungarian assignment; overlaps below 0.5 are reported as unmatched.
  - The legacy full-TD run kept no vectors, so it falls back to a cosine of the dominant-transition labels, which is coarser and labelled as such.
  - Acceptance (R, T): bright states (f ≥ 0.05) change by < 0.01 eV and the peak by < 2 nm. A, B and G are reported as measurements.
- **Comparison C** (CAM full vs TDA): not run. The CAM full-TD calculation is a compute DATA_LIMITATION (P1-ENV-1).
- **Tests:** `phase1/tests/test_compare_states.py` (3 passed, real formaldehyde runs):
  - self-overlap is the identity;
  - TDA ≥ RPA energies for matched states;
  - a permuted state order is recovered;
  - the label fallback works;
  - the tolerance-stage table is consistent.
- **Status:** IMPLEMENTED_AWAITING_PRODUCTION_RUN (needs L2_B3LYP_TDA15, L2_CAM_TDA15, L2_B3LYP_TDA25, L2_B3LYP_TDA15_RELAXED). Command: `python phase1/compare_states.py`.

## 8. Cell and ice physics

### P2-FILM-1: `film_only` also placed mineral dust in the algal film
- **Source:** `phase2/biosnicar_bridge.py::_layer_concs` (used by `run`/`run_multi`; reached via `phase4/emulator.py` and the `dust_film` scenario of `phase4/bias_study.py`).
- **Defect:**
  - With a split crust and `film_only=True`, every impurity, dust included, was concentrated into the 2 mm film.
  - Dust concentrations are bulk values for the surface-ice sample, so they belong uniformly in the 2 cm crust.
- **Repair:** `_layer_concs(..., in_film=None)`: cell counts (unit 1) follow `film_only`, other impurities stay uniform unless a caller places them explicitly.
- **Test:** `phase2/tests/test_phase2.py::test_film_only_places_algae_not_dust_in_film` (passed). It checks layer placement, column-amount conservation, and that BioSNICAR dust-only albedo is independent of `film_only`.
- **Before/after** (BioSNICAR, rds 1500 µm, SZA 47°; `records/repair_film_dust_sensitivity.json`):

  | dust | dusty-ice broadband albedo, before → after | algal Δα at 4×10⁴ cells mL⁻¹, before → after |
  |---|---|---|
  | 2×10⁵ ppb | 0.6453 → 0.6462 | 0.1444 → 0.1446 |
  | 5×10⁵ ppb | 0.5819 → 0.5874 | 0.1260 → 0.1276 |
  | 10⁶ ppb | 0.4936 → 0.5114 | 0.1006 → 0.1059 (+5 %) |

- **Affected outputs:** the `dust_film` row of `records/bias_study_williamson2020/` (to be regenerated with the corrected bridge under the held-out validation, task 6). Scenarios without dust, or without a film, are unchanged.
- **Status:** FIXED_AND_VERIFIED. Downstream regeneration: IMPLEMENTED_AWAITING_PRODUCTION_RUN.

## 11. Satellite retrieval

### P4-S2-1: BOA offset rule, grid georeferencing, illumination
- **Source:** `phase4/s2_io.py::_scale_offset`, `_read_scene`; `phase4/run_phase4.py::load_grid`.
- **Defects:**
  1. **BOA offset.** An explicitly recorded offset of 0 and a missing offset were treated the same: −0.1 was imposed whenever the baseline was ≥ 04.00.
     - Items whose DNs or metadata already account for the offset (`earthsearch:boa_offset_applied`) were not distinguished.
     - Local files were hard-coded as baseline "0", so baseline ≥ 04.00 local data would have received no offset at all.
  2. **Grid.** The window was snapped to the source pixel grid (`round_offsets`), but the output transform was built from the unsnapped bounds. Georeferencing was therefore wrong by up to half a source pixel.
     - Bounds in another CRS were not reprojected.
     - SCL (20 m) was resampled independently of the reflectance bands.
  3. **Illumination.**
     - The tile-mean sun elevation was used without checking it.
     - A missing value silently became 45°.
     - The local path defaulted to SZA 45°.
- **Repair:**
  - Explicit offset rules: recorded values are never overwritten; contradictory metadata raises; an unknown baseline raises. The local path needs `--processing-baseline` and `--sza`.
  - All bands and SCL are warped (GDAL WarpedVRT) onto one output grid whose transform is returned. Bounds must be whole pixels, and `aligned_bounds` snaps them outward onto the tile grid. A `dst_crs` different from the tile's is reprojected.
  - The solar zenith is computed at the AOI centre (NOAA approximation) and checked against the item's tile mean. The read raises above 1°; otherwise the AOI value is used.
  - Per-read provenance (scaling source, sun, grid, valid fraction) is stored in `scene.item["_read"]`.
- **Tests:** `phase4/tests/test_phase4.py` (4 passed):
  - `test_s2_scaling_and_slope`: every offset case, including "never −0.2".
  - `test_s2_reader_grid_alignment_and_mask`: synthetic COGs with a tile origin not on the 20 m grid; exact 2×2 averages; the SCL cloud columns line up exactly.
  - `test_s2_reader_reprojects_and_checks_sun`: reprojection to UTM 23N; a wrong acquisition time raises; a missing sun elevation raises.
  - `test_solar_zenith_matches_sentinel2_metadata`: the NOAA formula is within 0.3° of the Element 84 value at the tile centre.
- **Before/after on the scenes used:**
  - All scenes used (2017: baseline 00.01; 2019: 02.13) have recorded offset 0 and baseline < 04.00, so the offset repair changes no existing output (checked against live STAC metadata).
  - Reader comparison, 6 km S6 AOI at 20 m (`records/repair_s2_reader_comparison.json`):
    - old georeferencing error 0.53 px (10.7 m) north–south;
    - band means change by ≤ 3×10⁻⁴ reflectance;
    - SZA changes by 0.23–0.34° (tile mean → AOI), which does not move the integer SZA node used by the emulator for any of the three scenes tested.
  - The 10 m satellite validation window was already centred on the 10 m grid, so it is unchanged.
- **Remaining limitations:**
  - Atmospheric-correction residuals over bright ice (Sen2Cor), adjacency, and BRDF/topographic illumination are not corrected; they are absorbed by the per-band σ and the k nuisance and need the sensitivity analysis under task 6.
  - The emulator rounds SZA to whole degrees (≤ 0.5°).
- **Status:** FIXED_AND_VERIFIED.

### P2-PACK-1: is in-vivo packaging applied twice? (verified: no)
- **Source:** `phase2/cell_optics.py::CellModel.optics`, `phase2/empirical_data.py::pigment_macs_480`; the callers are `run_phase2`, `audit_optics` and `phase3/forward_model`.
- **Finding:**
  - The pigment spectra used are the *in vivo-form, unpackaged* specific absorption spectra (Williamson et al. 2020, after Dauchet et al. 2015, whose model applies cell packaging on top of them). The phenolic MAC is from solution extracts.
  - The packaging factor Q* is applied in exactly one place (`CellModel.optics`) per call. Tier A (Chevrollier 2023) uses measured cell optics with no further packaging.
  - No code path applies Q* to an already packaged spectrum.
- **Status:** FIXED_AND_VERIFIED (no defect).

### P2-VAC-1: `vacuole_fraction < 1`: compartments shaded independently
- **Source:** `phase2/cell_optics.py::CellModel.optics`.
- **Defect:** the vacuole pigment and the whole-cell absorbers (water, chloroplast pigments) each got their own self-shading factor, as if the other compartment were absent. The summed absorption could exceed the geometric cross-section: up to 1.27 × S/4 in the visible for f = 0.5.
- **Repair:** `pigment_packaging.joint_absorption` uses μ-random rays through the cell and a concentric vacuole of the same shape. Along each ray, τ = a_cell·l + a_vac·l_vac; the absorbed fraction 1 − e^−τ is shared in proportion to each compartment's optical depth. This is exact in ray optics for homogeneous compartments.
- **Test:** `test_vacuole_joint_absorption_bounded_and_consistent` (passed). It covers the f = 1 limit (equal to Q*·aV within 1 %), the dilute limit, the opaque bound ≤ S/4, and that the old treatment exceeds that bound.
- **Before/after** (cylinder 10 × 20 µm, calibrated tier C phenolics + chlorophylls/carotenoids; `records/repair_vacuole_joint_absorption.json`), mean 400–700 nm absorption:

  | vacuole fraction | change |
  |---|---|
  | f = 1 | unchanged (production setting) |
  | f = 0.5 | −29 % |
  | f = 0.25 | −23 % |

- **Affected outputs:** none in production (vacuole_fraction = 1 everywhere). Only `run_phase2 --vacuole-fraction` sensitivity runs.
- **Remaining limitation:** a concentric, same-shape vacuole is assumed; refraction at compartment boundaries is ignored.
- **Status:** FIXED_AND_VERIFIED.

## 9. Inference

### P4-INV-1: float32 likelihood
- **Source:** `phase4/inversion.py::GridPosterior`.
- **Finding:**
  - The quadratic form rr − 2G + a was evaluated in float32. Each term is about 10⁴ for bright ice at σ ≈ 0.01, so the O(1) result loses about 3 digits.
  - Measured on a production emulator (21×11×28×6, refined to 0.1 dex; 400 synthetic pixels, σ = 0.01–0.015; `records/repair_inversion_float64.json`), float32 vs float64:

    | quantity | max change |
    |---|---|
    | posterior mean log B | 4.5e-5 dex |
    | SD | 1.9e-4 |
    | q025/q975 | ≤ 1.1e-4 |
    | radius mean | 0.14 µm |
    | RF_algae | 0.003 W m⁻² |
    | log evidence | 1.2e-3 |

  - So the approximation was harmless at this σ, but not guaranteed for smaller σ.
- **Repair:** float64 throughout (runtime 8.9 → 12.6 s for 400 pixels).
- **Test:** `test_loglik_equals_multivariate_normal_and_is_float64_accurate`. It checks against scipy's multivariate normal with Σ = diag σ² + s_k² F Fᵀ and mean m_k F to 1e-8 (passed). This also verifies the k-covariance algebra.
- **Status:** FIXED_AND_VERIFIED.

### P4-INV-2: undefined derived quantities averaged as zero
- **Source:** `GridPosterior.run` (derived: pigment, bba, rf).
- **Defect:** NaN nodes were replaced by 0 and still carried posterior weight, which biased the means towards 0.
- **Repair:** renormalise over finite nodes, report `<q>_nonfinite_mass`, and return NaN if that mass exceeds 1e-3.
- **Test:** `test_derived_quantities_not_zero_filled` (passed).
- **Affected outputs:** maps of derived quantities wherever the emulator had NaN nodes. The production emulators checked have finite derived arrays (the float32/float64 comparison shows max bba change 6e-6).
- **Status:** FIXED_AND_VERIFIED.

### P4-INV-3: χ² at the MAP as a goodness-of-fit test
- **Source:** `GridPosterior.run` (`chi2`), `run_phase4.py`, `satellite_validation.py` ("0 % of pixels fail the χ² test").
- **Defect:**
  - χ² was evaluated at the MAP with k at its conditional mode and compared with χ²(4).
  - With 4 bands and 3–4 fitted states plus k, the residual has ≤ 0 degrees of freedom, so the test has almost no power.
  - "0 % fail" is therefore not evidence that the model fits.
- **Repair:**
  - New outputs: `mahal_map` (the k-marginal Mahalanobis distance) and `ppp`, a posterior predictive p-value Σ_z w(z) P(χ²₄ ≥ D(z)). The ppp is conservative and its low power is stated where it is printed.
  - `chi2` is kept for comparison.
  - The MASTER_SUMMARY statement is to be corrected in the final summary rewrite (task 7).
- **Test:** `ppp` is checked to lie in [0, 1] in `test_derived_quantities_not_zero_filled`. The power limitation is documented, not testable.
- **Status:** FIXED_AND_VERIFIED (code); the summary text is pending (UNRESOLVED until task 7).

### P4-INV-4: MCMC cross-check priors and convergence diagnostic
- **Source:** `inversion.mcmc_pixel`, `priors.mcmc_prior_params`.
- **Defects:**
  1. With the dust axis active, the MCMC had no dust prior: it was uniform in log10(dust + 100), the interpolation coordinate. The grid uses the measured log-normal.
  2. With prior `scale` ≠ 1, the MCMC f_n prior was not widened as it is on the grid.
  3. Split R-hat was computed across walkers of ONE ensemble. Walkers are coupled by the moves, so this is not a test of independent chains.
  4. Walkers were started from grid-posterior draws, so the "cross-check" was not independent of the grid.
- **Repair:**
  - ln dust is sampled with the grid's log-normal prior, and ln r likewise; both are mapped to the emulator coordinates explicitly.
  - `mcmc_prior_params` reproduces the grid priors, including the scale and dust prior.
  - `n_ensembles` (default 4) independent ensembles, each started from prior draws.
  - `rhat` is now across ensembles; `rhat_walkers` is kept for reference.
- **Tests (passed):**
  - `test_mcmc_prior_params_match_grid_priors`: moments agree at scale 1 and 1.5.
  - `test_mcmc_uses_dust_prior_and_independent_ensembles`: with uninformative data, the ln-dust posterior equals the truncated prior, and R-hat < 1.05 across 3 ensembles.
  - `test_split_rhat_detects_disagreeing_chains`.
- **Status:** FIXED_AND_VERIFIED. The slow end-to-end test `test_grid_vs_mcmc_and_calibration` (BioSNICAR emulator, now with 4 prior-started ensembles) was running at the time of writing; its result is recorded below when complete.

### P4-INV-5: spatial pooling reused a pixel's own data in its prior
- **Source:** `run_phase4.py` (`--spatial-pooling`, off by default and off in every recorded run).
- **Defect:** the smoothed first-pass map that sets each pixel's prior included that pixel (double counting).
- **Repair:** leave-one-out Gaussian smoothing (the kernel centre weight is removed).
- **Test:** none automated (an option not used in any output). Status: IMPLEMENTED_AWAITING_PRODUCTION_RUN (no production run uses it).

## 6. Provenance

### P0-PROV-1: cache reuse by row count, path or manual version
- **Sources:**
  - `phase3/run_phase3.py::run_or_load` (reused a CSV whenever the row count matched).
  - `phase2/tddft_calibration.py::cached_calibration` (keyed by source label and sticks; ignored settings such as seed and chain length, the molar mass, the width and the measured data).
  - `phase4/emulator.py::build_emulator` (tag held the Phase 1 *path* and a hand-bumped PHYSICS_VERSION; non-atomic `np.savez`).
- **Repair:** `phase2/provenance.py`.
  - Canonical content fingerprints (arrays by bytes, dataclasses, frames).
  - Code fingerprints of the physics modules, Phase 1 content fingerprints, the empirical-data fingerprint, the BioSNICAR git revision, and an environment record.
  - Atomic writes, and a `<file>.meta.json` sidecar with fingerprint + result checksum.
  - `run_or_load` reuses only on a matching sidecar.
  - `cached_calibration` keys on content + kwargs + data + code.
  - The emulator tag is a content fingerprint; `Emulator.save` is atomic; a corrupt cache (BadZipFile/EOF/OSError) triggers a rebuild.
  - emcee move RNGs are seeded (`sampler.random_state`) in the calibration and in `mcmc_pixel`.
- **Tests:** `phase2/tests/test_provenance.py` (6 passed):
  - same row count with a different design → recomputed (the old code reused it);
  - a settings change → recomputed;
  - a result edited after writing → recomputed;
  - different kwargs and molar mass → separate calibration entries;
  - the Phase 1 fingerprint follows content;
  - calibration chains are bit-identical for the same seed whatever the global NumPy state;
  - the emulator save leaves no temporary file.
- **Consequence:** existing caches (emulators, run tables) carry no new-style fingerprint and will be rebuilt by the next production run. This is intended: the physics code also changed in this repair.
- **Environment pin:** `python phase2/provenance.py` prints versions and code hashes. It is written into every sidecar.
- **Status:** FIXED_AND_VERIFIED.

## 12. Daily energy

### P4-DAY-1: daily-mean forcing formula, hard-coded overpass hour, completion, seasonal label
- **Source:** `phase4/multi_scene.py::summarise`, `main`.
- **Defects:**
  1. RF_daily = RF_overpass × mean₂₄(SW_meas) / SW_meas(15 UTC). But RF_overpass was computed with the emulator's *clear-sky model* SW, so the measured overpass SW in the denominator mixed two irradiances. The albedo reduction is RF / SW_model.
  2. The overpass was hard-coded as hour 15 (the actual times are 15:04–15:14 UTC), and the hourly-label convention was ignored.
  3. A scene counted as complete if `summary.json` existed.
  4. The mean over 6 clear-sky dates was labelled a "season" value.
  5. "Melt" was all extra energy into melt, without an SEB.
- **Repair:**
  - Δα = RF_overpass / SW_model(overpass), where SW_model is recorded by run_phase4 (`sw_down_model_w_m2`) or recomputed exactly as the emulator does. Then RF_daily = Δα × mean₂₄(SW_meas).
  - The measured overpass SW is interpolated at the STAC datetime (hour-centre convention) and used only as a clearness diagnostic.
  - run_phase4 writes `COMPLETE.json` last, with output checksums; multi_scene aggregates only marked runs. Legacy runs are accepted only with `--aggregate-only --accept-legacy` and are labelled per scene.
  - The output is labelled as the sampled-dates mean, and the melt column is renamed `melt_potential_*`.
- **Tests:** `test_multi_scene_daily_forcing_uses_model_sw_and_overpass_time` and `test_multi_scene_completion_requires_marker` (2 passed).
- **Before/after** (same six 2019 runs re-aggregated; `records/repair_multi_scene_daily/` vs `records/multi_scene_2019_tddft/`):

  | quantity | before | after |
  |---|---|---|
  | sampled-dates mean daily algal RF (ours) | 11.17 W m⁻² | 10.38 W m⁻² (−7 %) |
  | range | 6.23–13.61 | 6.26–12.54 |
  | Tier A | 14.58 | 13.49 |
  | potential melt (ours) | 0.289 cm w.e. d⁻¹ | 0.268 cm w.e. d⁻¹ |

  - The measured/model clearness at overpass was 0.90–0.98.
  - The 2019 KAN_M SW is the tilt-uncorrected `dsr`. That affects the daily mean directly: a DATA_LIMITATION to be bounded in the SEB task.
- **Remaining:**
  - The diurnal SZA dependence of Δα is approximated by its overpass value. It is quantified as < 10 % in literature_comparison and is to be replaced by the hourly SEB integration (task 7).
  - Actual melt needs the SEB (task 7).
- **Status:** FIXED_AND_VERIFIED (formula, completion, labels). Actual melt: IMPLEMENTED_AWAITING_PRODUCTION_RUN, pending task 7.

## 5. Fe(III) investigation

### P1-FE-1: Fe increment analysis inputs, naming and scope
- **Sources:** `phase1/fe_increment.py`, `phase1/fe_complex.py::spin_check`, `phase1/run_phase1.py` (xTB geometry status).
- **Defects:**
  - `eps_of` interpolated each run's own pre-broadened spectrum (written from 300 nm, while the analysis grid starts at 280 nm, so values below 300 nm were NaN). It did not check that ligand and complex shared functional, basis, solvent, TDA and width, and it accepted any folder with a `summary.json`.
  - `lmct_nm` named a state "LMCT" without any character analysis.
  - `spin_check` (vertical single points at the sextet geometry) was described as a ground-state check.
  - xTB-only geometries were not marked as exploratory.
  - The concentration assumption of the vis_ratio metric was not stated.
- **Repair:**
  - Spectra are rebuilt from the sticks with one FWHM on one shared 280–800 nm grid, with root-coverage diagnostics.
  - `completion.validate_run` is required (provisional or exploratory runs only with `--allow-provisional`, reported per row), and the settings must match the ligand run (otherwise it raises).
  - `lmct_nm` is renamed `strongest_vis_state_nm/_f`.
  - `spin_check` is documented as vertical only.
  - `--xtb-geometry` runs get geometry status `exploratory_xtb`, which is never production.
  - The docstring states the concentration caveat (the measured dA is a lower bound if complexation was incomplete; `shape_r` is concentration-free), and that nothing here feeds tier D, which uses the measured increment only.
- **Test:** `phase1/tests/test_fe_increment.py` (1 passed). It checks the 280 nm coverage, a basis mismatch → raise, a missing marker → raise, an exploratory run rejected by default and flagged when allowed, and the neutral column names.
- **Affected outputs:** none produced yet; the Fe TD-DFT jobs (FE_CAT, FE_TROP_XTB) are queued.
- **Status:** FIXED_AND_VERIFIED (code). Results: IMPLEMENTED_AWAITING_PRODUCTION_RUN.

### P1-ENV-2: OOM kill of the CAM-B3LYP TDA job by concurrent analysis work (incident)
- **What happened:** the runner admitted two 7 GB Level 2 jobs (memory budget: total − 1.5 GB), which left about 1.2 GB for everything else. My concurrent analysis processes (calibration comparison, held-out run, SEB validation) pushed the container cgroup over its limit.
  - The kernel killed `L2_CAM_TDA15` (RSS 7.06 GB; `dmesg`: memory cgroup out of memory).
  - The job had not reached its first TD checkpoint (it was still in SCF/density-fitting setup). About 50 min of work was lost; no checkpoint or result was corrupted.
- **Repair:**
  - `jobs.MEM_RESERVE_MB = 3000` keeps 3 GB free for analysis work. Two 7 GB L2 jobs therefore no longer run together; they run sequentially.
  - The runner was restarted with `--retry-failed`. The running jobs, which run in their own sessions, were not touched, and `L2_CAM_TDA15` was requeued; it waits for memory.
  - Analysis jobs now run at `nice 10` with bounded chunk memory.
- **Test:** `test_jobs.py` still passes (2 passed). The reserve is a configuration value.
- **Status:** FIXED_AND_VERIFIED (cause and policy). The CAM result: IMPLEMENTED_AWAITING_PRODUCTION_RUN.

## 14/12. Offline surface energy balance (SEB) integration

### P4-SEB-1: actual vs potential melt (new capability)
- **Source:** new `phase4/seb.py`. Forcing: PROMICE KAN_M hourly data (GEUS Dataverse doi:10.22008/FK2/IW73UU, file `KAN_M_hour.csv`, sha256 bb34825a…, CC BY 4.0), June–August 2016–2019, subset `data/empirical/promice_KAN_M_hour_JJA_2016_2019_seb.csv`.
- **Model:**
  - Hourly point SEB: measured SW↓, LW↓, T_a, q_a, U and p; bulk turbulent fluxes with a Richardson-number stability correction, z0 = 1 mm.
  - A 0.10 m ice slab carries cold content. Melt occurs only at 0 °C with positive net flux.
- **Tests:** `phase4/tests/test_seb.py` (4 passed):
  - energy closure at the melting point;
  - zero net flux and no melt for a cold surface;
  - a zero-thickness slab equals the instantaneous balance, and the slab delays melt;
  - paired runs give exactly 0 without algae, and 0 < actual ≤ potential with algae.
- **Validation against independent observations** (`records/seb_2019/seb_validation_KAN_M.json`; bare-ice days with complete data):
  - Turbulent fluxes vs GEUS's own estimates (2017): bias +1.5 / +1.3 W m⁻²; r = 0.78 / 0.84.
  - Modelled surface temperature vs the temperature from measured LW↑ (0.1 m slab): bias −0.05 to +0.04 K, RMSE 0.76–0.87 K. Without the slab: −0.4 to −0.8 K and 2.2–3.3 K, which is why the slab is used.
  - Melt vs the measured ice-surface lowering (ρ_ice 900): modelled/observed totals are 0.61/0.38 (2016), 0.44/0.27 (2017), 0.60/0.36 (2018) and 1.33/0.43 m w.e. (2019); daily r = 0.49–0.80.
  - **Closure check:** melt energy from fully MEASURED fluxes (net SW and LW radiometers + GEUS turbulent fluxes) also exceeds the observed lowering, by 1.45–1.58× (2016–2018) and 3.0× (2019; uncorrected `dsr`, no tilt correction). The SEB reproduces the measured-flux melt to within about 10 %. The gap therefore lies in the observations (radiometer tilt, ablation sensor, density) rather than in the SEB physics, but it is not resolved: absolute melt is uncertain by a factor of 1.5–3 (**DATA_LIMITATION**).
- **Algal coupling without double counting:** paired runs on the retrieved S6 surface (2019 sampled dates, interpolated, 8 Jul–29 Aug).
  - "On" uses the retrieved bare-ice BBA, which contains the algae; "off" adds back the retrieved algal albedo reduction. Meteorology is identical.
  - Results (`records/seb_2019/paired_algae_2019.json`):

    | optics | actual algal melt (m w.e.) | potential (m w.e.) | actual / potential |
    |---|---|---|---|
    | pigment-aware (tddft_D) | 0.108–0.116 | 0.126 | 0.86–0.92 |
    | Tier A | 0.143–0.152 | 0.163 | 0.88–0.93 |

    The ranges cover z0 1e-4–1e-2 m and slab 0.05–0.3 m.
  - The ratio is the robust result: the potential-melt shortcut overstates algal melt by 8–14 %.
  - The absolute values inherit the 2019 SW limitation and the closure gap above.
- **Remaining limitations:**
  - KAN_M (1270 m) meteorology is used for the S6 (~1000 m) surface.
  - The albedo between sampled dates is interpolated.
  - Cloudy days are represented by the measured meteorology, but the algal Δα comes from clear-sky retrievals.
- **Status:** FIXED_AND_VERIFIED (SEB code and tests). The actual-melt numbers are reported with the DATA_LIMITATION above.

## 10. Held-out predictive evaluation

### P4-HO-1: the "out-of-sample" field validation was not held out by site
- **Source:** `phase4/field_validation.py::run`.
- **Defects:**
  - Hyper-parameters (σ, radius prior) were chosen by leave-one-out over BOTH sites pooled.
  - τ was fitted leave-one-out on pooled samples.
  - Zero counts were dropped.
  - Gaussian observation errors were used.
  - The baseline regression was trained on the pooled data.
  - The model comparison was by empirical-Bayes evidence, which is not a held-out score.
- **Repair:** new `phase4/heldout.py`.
  - Two folds: S6 → S Greenland (primary) and the reverse.
  - Everything is fitted on the training site only: σ and the radius prior by training-site evidence; τ by maximum likelihood of the training observations using inner leave-one-out posteriors.
  - Test predictive = posterior marginal ⊗ N(0, τ²).
  - Poisson count likelihood with zeros kept (S6, known counted volume; zeros at V = 0.016 mL); log-normal observation error with an assumed SD of 0.10 dex for S Greenland (count numbers not published).
  - Four frozen optical models (Tier A; measured in-vivo MAC tier C; TD-DFT tier C; TD-DFT tier D).
  - Baselines on the same training data and scoring: climatology, the literature prior, band ratio and ridge (λ by training LOO), plus the published Cook et al. (2020) values on S6 as a reference.
  - Forward test: log p(4-band HCRF | measured abundance), and the broadband HCRF error against the measured spectrum (350–2500 nm).
  - Paired bootstrap CIs.
- **Tests:** `phase4/tests/test_heldout.py` (5 passed):
  - zero counts and the counted volume in the likelihood;
  - the score is a normalised density score;
  - the Gaussian baseline recovers the truth, and dropping zeros biases it (the old behaviour);
  - test-site evidence cannot change the training-site settings;
  - the vectorised tables equal the reference `GridPosterior.run` (evidence, marginal, forward density).
- **Status:** IMPLEMENTED_AWAITING_PRODUCTION_RUN (running: `phase4/results/heldout_v1`). Results are recorded below when complete.

### P4-SAT-2: satellite validation: pixel identity, leakage, E[log B] vs log E[B], sub-pixel mixing
- **Source:** `phase4/satellite_validation.py` (`invert_scene`, `within_pixel_sd`, site metrics, `main`); `phase4/field_validation.py::run`.
- **Defects:**
  1. Pixel ids were window-relative "row_col", and the two primary scenes have different windows. `within_pixel_sd` therefore pooled samples from different places as if they shared a pixel. In the recorded tier-D run: 9 merged ids vs 14 true scene-pixels. The SD was 0.835 dex with 11 "dof", against 0.841 dex with 6 dof when corrected.
  2. The "field calibration" applied to the pixels (σ, radius prior, τ) was fitted on all S6 counts, including the 20 samples then used to score the pixels (label leakage into coverage).
  3. The site-level comparison used the mean of log counts (geometric mean; zeros dropped). A pixel sees area-mixed reflectance.
- **Repair:**
  - Pixel id = scene + absolute pixel-centre coordinates. `within_pixel_sd` refuses ids without them and reports a 95 % χ² interval.
  - `FV.run(exclude=…)` keeps the satellite-compared samples out of the hyper-parameter and τ fits.
  - The site metrics now report both the geometric mean and log10 of the arithmetic mean (zeros included), with biases against each.
  - New `subpixel_mixing` maps a within-pixel abundance spread to the value a pixel retrieval should return.
- **Tests:** `test_within_pixel_sd_does_not_merge_pixels_across_scenes` and `test_subpixel_mixing_limits` (2 passed).
- **Measured** (`records/repair_satellite_scale.json`; tier D emulator, 21 Jul 2017; patches log10 B ~ N(3.6, sd); linear areal mixing):

  | within-pixel SD (dex) | retrieved | geometric mean | arithmetic mean |
  |---|---|---|---|
  | 0.40 | 3.76 | 3.59 | 3.78 |
  | 0.83 (measured) | 4.09 | 3.59 | 4.34 |
  | 1.20 | 4.31 | 3.58 | 4.80 |

  **Consequence:** at the measured spread, an unbiased pixel retrieval should read about +0.5 dex above the mean log count. The recorded site-level "agreement within 0.09–0.27 dex" against the geometric mean is therefore not evidence of an unbiased retrieval: relative to the mixing expectation it indicates a low bias of roughly 0.2–0.4 dex, or a failure of the linear-mixing assumption. That claim must be withdrawn from the summary (task 7).
- **Affected outputs:** `records/satellite_validation_tddft_tierD/`, `records/satellite_validation_williamson2020/` (legacy; to be re-run with the repaired code, which needs network reads of the 2017 scenes).
- **Status:** FIXED_AND_VERIFIED (code and analysis). Re-run: IMPLEMENTED_AWAITING_PRODUCTION_RUN.

## 7. Calibration statistics

### P2-CAL-1: stage-2 Gaussian approximation, φ clipping, covariance fallback, correlated residuals
- **Source:** `phase2/tddft_calibration.py::_fit_magnitude`, `calibrate`.
- **Defects:**
  1. Stage-2 draws of (ln f, φ) came from a Laplace approximation; draws outside [0, 1] were clipped, piling mass on the bounds.
  2. A non-positive-definite Hessian silently fell back to cov = diag(1e-4).
  3. Optimizer success was never checked.
  4. The discrepancy SD ln s was fixed at its MLE, not propagated: its posterior SD was reported as 0.005.
  5. The residuals at 2 nm spacing were treated as independent, but they are strongly autocorrelated. Lag-1 autocorrelation is 0.99; the effective sample sizes are about 13 of 168 (HPLC shape) and about 5 of 246 (extract MAC), from `records/repair_calibration_residual_acf.json`. All posterior widths were therefore far too narrow.
- **Repair:**
  - Stage 2 is now exact grid quadrature over (ln f, φ ∈ [0, 1], ln s ∈ [−8, 3]) with flat priors. Edge-mass diagnostics are recorded and ln s is propagated.
  - `_fit_magnitude` raises instead of falling back; the MAP optimizer success is recorded.
  - Residual model `residuals='ar1'` (default): AR(1) correlation between neighbouring wavelengths, ρ by maximum likelihood for each data set. `'iid'` is kept as a structural scenario. The residual model is an emulator setting (`EmulatorConfig.cal_residuals`), so both enter the held-out comparison.
  - `stride` thins the data for sensitivity tests.
- **Tests:**
  - `test_provenance.py::test_calibration_chains_reproducible` (passed) exercises the new stage 2.
  - The full production calibrations below were executed.
  - `test_fit_magnitude_has_no_silent_covariance_fallback` (passed): a failed optimizer raises, and its status is reported.
- **Before/after on the production spectrum** (`records/repair_calibration_stage2.json`; mean ± SD):

  | variant | ΔE (eV) | f | φ |
  |---|---|---|---|
  | old (Laplace + clip, iid) | 0.061 ± 0.004 | 39.5 ± 2.2 | 0.49 ± 0.04 |
  | grid, iid | 0.061 ± 0.004 | 39.3 ± 2.3 | 0.50 ± 0.04 |
  | grid, iid, every 5th point | 0.061 ± 0.004 | 38.6 ± 5.1 | 0.52 ± 0.10 |
  | grid, iid, every 10th point | 0.061 ± 0.004 | 37.3 ± 6.7 | 0.57 ± 0.14 |
  | grid, AR(1) (ρ_shape 0.95, ρ_mac 0.995) | 0.041 ± 0.008 | 69 ± 27 | 0.11 ± 0.02 |
  | grid, AR(1), every 5th point (ρ_mac 0.97) | 0.041 ± 0.008 | 61 ± 41 | 0.16 ± 0.06 |

  - The iid posterior is not stable under thinning; the AR(1) posterior is, approximately.
  - The two residual models disagree on the split between f and φ. The Fe-complexed share of visible absorption is 0.53 (iid) vs 0.19 (AR(1)), while the visible tier-D MAC differs by only +10 % (51k vs 57k m² kg⁻¹; extract 68k).
  - The intensity factor f and φ are therefore **not uniquely attributable** with these data. That is a stated result, not something to be tuned away.
- **Remaining:**
  - ρ_mac sits at the grid edge (0.995); the discrepancy is almost random-walk-like.
  - The extract's two 2 nm channels are the only constraint on level vs shape.
  - The choice between residual models is tested on held-out field data (P4-HO-1: tddft_D vs tddft_D_iid), not by preference.
- **Status:** FIXED_AND_VERIFIED (statistics). Structural ambiguity: DATA_LIMITATION.

### P4-INV-6: grid refinement targets
- **Source:** `run_phase4.py --log-b-step` (default 0.1 dex), the grid posterior.
- **Measurement** (`records/repair_grid_refinement.json`; 150 synthetic pixels, tier-D emulator; reference step 0.025 dex), max |difference|:

  | step | log B mean | 95 % bounds | RF_algae | radius mean |
  |---|---|---|---|---|
  | 0.1 | 0.012 dex | 0.07 dex | 0.80 W m⁻² | 33 µm |
  | 0.05 | 0.003 dex | 0.021 dex | 0.007 W m⁻² | 7 µm |

- **Targets set here:** posterior mean < 0.01 dex and interval bounds < 0.03 dex. Only 0.05 meets them, so the map default is changed from 0.1 to 0.05. The field validation and the held-out experiment already use 0.05.
- **Affected outputs:** the legacy maps (0.1 dex) have interval-bound discretisation errors of up to 0.07 dex. They are to be regenerated.
- **Status:** FIXED_AND_VERIFIED (default). Map regeneration: IMPLEMENTED_AWAITING_PRODUCTION_RUN.
