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

### P4-SEB-1: modelled vs potential melt (new capability) — SUPERSEDED by P4-SEB-2 (numbers below withdrawn)
- **Source:** new `phase4/seb.py`. Forcing: PROMICE KAN_M hourly data (GEUS Dataverse doi:10.22008/FK2/IW73UU, file `KAN_M_hour.csv`, sha256 bb34825a…, CC BY 4.0), June–August 2016–2019, subset `data/empirical/promice_KAN_M_hour_JJA_2016_2019_seb.csv`.
- **Model:**
  - Hourly point SEB: measured SW↓, LW↓, T_a, q_a, U and p; bulk turbulent fluxes with a Richardson-number stability correction, z0 = 1 mm.
  - A 0.10 m ice slab carries cold content. Melt occurs only at 0 °C with positive net flux.
- **Tests:** `phase4/tests/test_seb.py` (4 passed):
  - energy closure at the melting point;
  - zero net flux and no melt for a cold surface;
  - a zero-thickness slab equals the instantaneous balance, and the slab delays melt;
  - paired runs give exactly 0 without algae, and 0 < modelled increment ≤ potential with algae.
- **Validation against independent observations** (`records/seb_2019/seb_validation_KAN_M.json`; bare-ice days with complete data):
  - Turbulent fluxes vs GEUS's own estimates (2017): bias +1.5 / +1.3 W m⁻²; r = 0.78 / 0.84.
  - Modelled surface temperature vs the temperature from measured LW↑ (0.1 m slab): bias −0.05 to +0.04 K, RMSE 0.76–0.87 K. Without the slab: −0.4 to −0.8 K and 2.2–3.3 K, which is why the slab is used.
  - Melt vs the measured ice-surface lowering (ρ_ice 900): modelled/observed totals are 0.61/0.38 (2016), 0.44/0.27 (2017), 0.60/0.36 (2018) and 1.33/0.43 m w.e. (2019); daily r = 0.49–0.80.
  - **Closure check:** melt energy from fully MEASURED fluxes (net SW and LW radiometers + GEUS turbulent fluxes) also exceeds the observed lowering, by 1.45–1.58× (2016–2018) and 3.0× (2019; uncorrected `dsr`, no tilt correction). The SEB reproduces the measured-flux melt to within about 10 %. The gap therefore lies in the observations (radiometer tilt, ablation sensor, density) rather than in the SEB physics, but it is not resolved: absolute melt is uncertain by a factor of 1.5–3 (**DATA_LIMITATION**).
- **Algal coupling without double counting:** paired runs on the retrieved S6 surface (2019 sampled dates, interpolated, 8 Jul–29 Aug).
  - "On" uses the retrieved bare-ice BBA, which contains the algae; "off" adds back the retrieved algal albedo reduction. Meteorology is identical.
  - Results (`records/seb_2019/paired_algae_2019.json`):

    | optics | modelled algal melt increment (m w.e.) [withdrawn] | potential (m w.e.) | increment / potential [withdrawn] |
    |---|---|---|---|
    | pigment-aware (tddft_D) | 0.108–0.116 | 0.126 | 0.86–0.92 |
    | Tier A | 0.143–0.152 | 0.163 | 0.88–0.93 |

    The ranges cover z0 1e-4–1e-2 m and slab 0.05–0.3 m.
  - [Withdrawn wording: an earlier version called this ratio robust and called the 8–14 % gap a correction. It was neither. The calculation had the longwave error, silent gap filling and shared-surface Tier A described in P4-SEB-2, and the ratio depends strongly on the untested subsurface-shortwave assumption.]
  - The absolute values inherit the 2019 SW limitation and the closure gap above.
- **Remaining limitations:**
  - KAN_M (1270 m) meteorology is used for the S6 (~1000 m) surface.
  - The albedo between sampled dates is interpolated.
  - Cloudy days are represented by the measured meteorology, but the algal Δα comes from clear-sky retrievals.
- **Status:** superseded by P4-SEB-2.

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

### P1-PROT-1: protonation-state sensitivity (previously unaudited)
- **Finding:** the Level 2 pigment was computed only as the neutral acid. The carboxylic acid (pKa ≈ 3–4) is probably deprotonated at vacuolar pH. Tautomers, conformers, a larger basis set and explicit solvent are also unaudited.
- **Action:**
  - New molecule `level2_carboxylate` (C18H17O12, charge −1). Its start geometry is the neutral L2 endpoint with the carboxylic H removed (`phase1/results/v2/inputs/level2_carboxylate_start.xyz`; atom count and charge checked).
  - New job `L2_COO_OPT_TDA15`: B3LYP/PCM optimisation (60 steps) followed by TDA, 15 roots, tolerance 1e-5.
  - Comparison P in `compare_states.py`: neutral vs carboxylate, states matched by overlaps.
- **Queue note:** the running runner process loaded the job table before this job existed. It is picked up when the runner is next started (`python phase1/jobs.py run`).
- **Tautomers, conformers, basis set, explicit solvent:** not scheduled. They need several Level 2 optimisations (days of CPU on this container). This is a **DATA_LIMITATION (compute)**, listed as such in the summary.
- **Status:** IMPLEMENTED_AWAITING_PRODUCTION_RUN.

## 13. Uncertainty propagation

### P3-UNC-1: independent calibration marginals, iid bootstrap of an LHS, structural scenarios
- **Source:** `phase3/param_space.py::default_parameters`, `phase3/stats_tools.py::describe`, `phase3/run_phase3.py`.
- **Defects:**
  - The four calibration parameters (ΔE, FWHM, f, φ) were sampled as independent marginals, which dropped their strong posterior correlation (f–φ).
  - The CIs of Monte Carlo statistics came from an iid bootstrap of a Latin-hypercube sample, which is not an iid sample.
  - The calibration's residual model (a structural choice) was not represented.
- **Repair:**
  - A single input `cal_draw ~ U(0, 1)` indexes a joint posterior row (ΔE, w, f, φ together). The Sobol inputs stay independent, and the molecular index is that of the joint draw. `--independent-calibration` keeps the old behaviour for sensitivity runs.
  - The MC design is R = 10 independent LHS replicates by default; CIs of the mean and median use the replicate spread (Student t), labelled `ci_method`.
  - `--cal-residuals {ar1, iid}` runs the structural scenario as a separate set, not mixed into the probabilistic spread.
- **Tests:** `phase3/tests/test_phase3.py::test_joint_calibration_draw_keeps_correlations_and_replicate_ci` (passed). It checks the parameter sets, that a draw maps to one joint row (correlation −0.9 preserved), and that the replicate CI covers the truth.
- **Double counting:** the field-fitted model error τ is applied only to retrievals (Phase 4), not to the Phase 3 forward forcing distribution, so it is not counted twice. Abundance enters Phase 3 as the measured S6 distribution, not as retrieval output.
- **Status:** FIXED_AND_VERIFIED (code). Phase 3 production re-run: IMPLEMENTED_AWAITING_PRODUCTION_RUN.

### P2-HARM-1: harmonisation of model settings across phases (verified)
- **Checked:**
  - Phase 2 (`run_phase2`), Phase 3 (`forward_model`) and Phase 4 (`emulator`) use the same MAC window (350–800 nm), UV rule ("hold"), cell asymmetry mode (Mie), cell refractive indices, water absorption, chlorophyll/carotenoid set, and the BioSNICAR revision fe74eeef.
  - Phase 3 and Phase 4 use bubbly ice with the same density defaults (450/690 kg m⁻³).
- **Stated differences:**
  - Phase 3 uses one population-mean cell, with size and concentration distributions as parameters; Phase 4 uses the two measured species mixed by f_n.
  - The calibration is fitted on 260–750 nm, while the optics use 350–800 nm: the 750–800 nm part is an extrapolation of the calibrated band model (MAC there is < 1 % of its visible mean).
- **Status:** FIXED_AND_VERIFIED (consistency check; no defect).
- **Joint vs independent calibration draws** (`records/repair_joint_calibration_draws.json`; visible 400–700 nm tier-D MAC from 3000 draws):
  - AR(1) calibration (f–φ correlation −0.24): joint SD 21.7k vs independent 22.5k m² kg⁻¹, on a mean of 57.7k. The two are almost the same; the visible MAC uncertainty is about ±38 %.
  - iid calibration (f–φ correlation −0.81): joint SD 1.73k vs independent 3.63k, on a mean of 52.0k. Independent marginals doubled the spread.
  - The joint draws therefore matter under the iid scenario. Under AR(1), the posterior width itself dominates.


### P4-SEB-2: SEB corrections after review of 453e751; investigation of the 2019 ablation discrepancy; Tier A on its own albedos
- **Source:** `phase4/seb.py`, `phase4/multi_scene.py`.
- **Defects found in P4-SEB-1:**
  1. **Longwave.** All of LW↓ was absorbed. The correct term is ε(LW↓ − σT⁴); the error was +(1 − ε)LW↓ ≈ +5–6 W m⁻².
  2. **Time axis.** Gaps up to 6 h were interpolated, including shortwave. Remaining missing SW was set to 0. `paired_algae` dropped incomplete rows, so the time axis collapsed across gaps and the slab was stepped across them.
  3. **Albedo.** It was clipped to [0.05, 0.95] and missing daily albedo was filled with 0.6, both silently.
  4. **Convergence and closure.** Neither timestep convergence nor energy closure had been demonstrated.
  5. **Tier A pairing.** The "on" albedo for Tier A was the pigment-aware retrieval's albedo (a shared surface). Tier A was therefore an attribution of a different Δα to the same surface, not a Tier A forecast.
- **Repair:**
  1. LW_net = ε(LW↓ − σT_s⁴).
  2. Gaps and timestamps:
     - The forcing is placed on its own complete hourly axis (missing hours are kept as rows).
     - Only whole gaps of ≤ 2 h are interpolated, and never for shortwave; any missing shortwave makes the hour missing.
     - The slab restarts after a gap.
     - Counts of filled hours, missing hours and restarts are reported.
  3. Albedo must lie in [0, 1] wherever forcing exists (otherwise it raises). Days without measured albedo are excluded explicitly.
  4. Explicit sub-stepping (`n_sub`) and an energy budget (Σ net flux·dt = melt energy + Δ slab heat content + restart adjustments).
  5. Each optical model uses its own retrieved bare-ice BBA ("on") and its own Δα ("off" = on + Δα). `multi_scene` now aggregates `tierA_bba_mean`.
- **Tests:** `phase4/tests/test_seb.py` (7 passed):
  - LW closure, including the size of the former error;
  - cold-surface zero flux;
  - energy closure to 1e-6 and sub-step convergence (< 0.2 % at 16 sub-steps);
  - slab limits;
  - actual time axis with gaps; missing SW not zeroed; only whole short gaps filled;
  - albedo bounds;
  - paired runs give 0 for equal albedo and are ≤ potential.
- **Real-data checks** (2019, 8 Jul–29 Aug):
  - Energy closure residual relative to the energy input: 1e-15.
  - The algal increment changes by 0.4 % from 1 to 16 sub-steps (pigment-aware 0.1127 → 0.1122 m w.e.; Tier A 0.1485 → 0.1479), so hourly stepping is adequate.
- **Validation against each ablation record** (`records/seb_2019/seb_validation_by_sensor.json`). Daily sums on complete bare-ice days; model total / observed total (n days):

  | year | z_ice_surf (combined) | z_pt_cor (pressure transducer) | z_stake_cor (sonic on stake) |
  |---|---|---|---|
  | 2016 | 1.66 (17) | 1.54 (17) | 1.56 (11) |
  | 2017 | 1.67 (25) | 1.72 (25) | no data |
  | 2018 | 1.59 (38) | 1.49 (40) | insufficient (2) |
  | 2019 | 3.03 (50) | 3.07 (50) | **1.23** (50) |

  Surface-temperature bias is −0.13 to −0.02 K.
- **2019 investigation** (PROMICE readme `data/promice_aux/AWS_data_readme.pdf`, sha256 dea361c4…; sensor/QC extract `data/promice_aux/KAN_M_hour_MJJAS_2016_2019_qc.csv`):
  - **Radiation.** `tilt_y` is missing for all of 2019, so no `dsr_cor`/`usr_cor` exist. The 2019 albedo is uncorrected usr/dsr (r = 0.9999 with usr/dsr). The 2019 radiation is therefore not tilt-corrected (tilt_x −2.3 ± 1.9°).
  - **Surface height.** Per the readme, `z_ice_surf` follows the pressure transducer (`z_pt_cor`) during the ablation season, is adjusted manually at maintenance jumps, and is not a direct observation across gaps.
  - **Sensors disagree in 2019.** Over June–August the pressure transducer gives 0.57 m of lowering, z_ice_surf 0.75 m and the independent stake sonic ranger 2.12 m. The SEB matches the stake record (1.23) much better than the transducer (3.07). The 3× discrepancy is therefore specific to the transducer-based record, not a general SEB failure; which 2019 sensor is right is not determined here.
  - **Omitted process tested:** shortwave deposited below the surface (penetration into bubbly ice / weathering crust), which melts ice internally without lowering the surface (`sw_subsurface_frac` χ; `records/seb_2019/seb_subsurface_sw_test.json`).
    - χ = 0.3 brings the 2016, 2017 and 2018 ratios to 1.06, 0.94 and 0.98 (all records consistent).
    - With the same χ, 2019 lies between the two sensors (1.94 against the transducer, 0.79 against the stake).
    - χ is estimated from these same 2016–2018 ablation data. It is consistent with, but does not prove, subsurface absorption: a radiometer bias of similar size, or ablation-sensor footprint/representativeness, cannot be separated with these data. The measured-flux closure (P4-SEB-1) shows the same 1.5× gap.
  - **Footprint:**
    - The downward radiometers (2.5–3 m height) see a footprint of tens of m², the pressure transducer a single borehole, the stake sonic another point.
    - A comparison of station albedo with the Sentinel-2 pixel containing KAN_M is not done here (it needs the 2019 scene reads): **UNRESOLVED**.
- **Modelled algal melt increment, 2019 sampled period** (`records/seb_2019/paired_algae_2019_v2.json`; each model on its own albedos; ranges over z0 1e-4–1e-2 m and slab 0.05–0.3 m):

  | optics | χ | increment (m w.e.) | increment / potential | surface melt with algae (m w.e.) |
  |---|---|---|---|---|
  | pigment-aware (tddft_D, legacy calibration) | 0 | 0.108–0.116 | 0.85–0.92 | 1.02–1.27 |
  | pigment-aware (tddft_D, legacy calibration) | 0.3 | 0.071–0.076 | 0.57–0.60 | 0.58–0.83 |
  | Tier A (own albedos) | 0 | 0.143–0.151 | 0.88–0.93 | 1.11–1.36 |
  | Tier A (own albedos) | 0.3 | 0.096–0.101 | 0.59–0.62 | 0.64–0.89 |

  - These are **conditional modelled estimates of surface melt**. They depend on the SEB, KAN_M meteorology used for the S6 surface, uncorrected 2019 radiation, clear-sky retrievals interpolated between dates, and χ.
  - The subsurface-energy assumption alone changes the increment/potential ratio from about 0.9 to about 0.6. With χ > 0, part of the withheld energy may still melt ice internally; that internal melt is not included in the "surface melt" increment.
  - **No independent observation validates these increments.** They are not "actual algal melt", and the potential-melt gap is not a validated "correction".
  - The pigment-aware column uses the legacy (pre-review) calibration through the legacy 2019 retrievals. It is to be regenerated with the frozen primary calibration.
- **Status:** SEB code FIXED_AND_VERIFIED (tests, closure, convergence). Melt increments are conditional modelled estimates (DATA_LIMITATION: no independent algal-melt observation; 2019 ablation sensors disagree by 2.5×). Footprint comparison: UNRESOLVED.

### P4-HO-2: held-out scoring repairs after review of 453e751
- **Source:** `phase4/heldout.py`.
- **Defects:**
  1. `summary_stats` used the module default observation SD (0.10 dex) instead of the requested `--obs-sd-sgris`, so sensitivity runs changed the log score but not the coverage, CRPS or point metrics.
  2. Metrics referred to different targets and sets: the log score covered all samples via the count likelihood; coverage used the observation predictive; CRPS used the latent predictive; RMSE used the latent median on positives.
  3. Ridge standardisation was fitted on the whole training site and reused inside the inner LOO, so preprocessing leaked into the λ choice.
  4. The paired bootstrap resampled samples as if independent, although samples from one sampling day are correlated.
  5. The broadband comparison set measured HCRF over 350–2500 nm against k × BBA over 300–2500 nm, and the columns were labelled as albedo-like.
- **Repair:**
  1. The requested SD reaches every metric through `observation_sd`.
  2. All metrics use the predictive distribution of the observation. The log score is reported on all samples and on the positive subset that the other metrics use.
  3. Standardisation is fitted inside each inner fold and on the full training set for the final fit.
  4. `block_bootstrap` resamples sampling days and flags fewer than 5 blocks.
  5. The measured broadband HCRF is on the same 300–2500 nm support (300–350 nm held at the 350 nm value; its irradiance share is reported per sample). Columns are renamed `bb_hcrf_*` and `logp_hcrf_bands`, and the docstring states HCRF vs albedo and the role of k.
  - Primary model and primary contrast are taken from the pre-registration; contrasts are labelled `preregistered`.
  - The expensive tables are fingerprinted by their own code (`table_code_fingerprint`), not by the scoring code.
- **Tests:** `phase4/tests/test_heldout.py` (8 passed). New ones: the requested observation SD changes coverage, CRPS and log score; every ridge fit sees features standardised on its own rows; the block bootstrap groups by day and flags 2 blocks.
- **Status:** FIXED_AND_VERIFIED (code). Production run: IMPLEMENTED_AWAITING_PRODUCTION_RUN. The superseded run (pre-review code, old calibration) was stopped before any score was computed; it touched no chemistry job.

### P2-CAL-2: stage-2 sampling without clipping; ρ propagated; refinement verified
- **Source:** `phase2/tddft_calibration.py::_stage2_grid`, `calibrate`.
- **Defects:**
  - φ draws were jittered ±half a cell and then clipped to [0, 1], putting point masses at the bounds. End cells had the weight of full cells.
  - ρ was fixed at its maximum-likelihood value, which sat on the grid edge (0.995), so the residual-correlation uncertainty was not propagated.
- **Repair:**
  - Midpoint-rule cells, with the end cells halved (node mass = density × cell width); draws are uniform within the chosen cell, so they stay in [0, 1] without clipping.
  - Stage 1 samples ρ_shape jointly with (ΔE, w, ls), with a flat prior on [0, 0.999].
  - Stage 2 marginalises ρ_mac over the grid nodes (extended to 0.998 and 0.999) that carry non-negligible marginal likelihood.
- **Tests:**
  - `test_stage2_no_boundary_pileup_and_quadrature_refinement` (no draws exactly at 0 or 1; the refined grid agrees).
  - `test_provenance.py`: all 8 passed after the edit (335 s).
- **Refinement on the production spectrum** (`records/repair_stage2_refinement.json`; ΔE 0.041, w 0.61; grid ×1 vs ×2):

  | quantity | grid ×1 | grid ×2 |
  |---|---|---|
  | ln f mean (posterior SD 0.76) | 4.121 | 4.130 |
  | φ mean | 0.1035 | 0.1034 |
  | φ 95 % interval | 0.0637–0.1520 | 0.0639–0.1506 |

  The quadrature is converged for the quantities used downstream.
- **Finding:** the ρ_mac posterior mass sits at 0.998–0.999, the top of the grid. The residual process is effectively a random walk: AR(1) is at its limit and the discrepancy is dominated by a smooth level/shape mismatch. The f and φ posteriors are conditional on this discrepancy model; a smoother discrepancy model (e.g. a Gaussian process with fitted length scale) is not implemented (**UNRESOLVED**, stated as a model-structure limitation).
- **Status:** FIXED_AND_VERIFIED (sampling, propagation, refinement). Discrepancy-model adequacy: UNRESOLVED.

### P1-ENV-3: second OOM kill and a full disk (incident, 12:18–12:37 UTC)
- **What happened:**
  - At about 12:18 UTC the kernel OOM-killed `L2_B3LYP_TDA15` (7.0 GB) and `L2_OPT` (5.5 GB). Two of my analysis processes were running at the same time: the calibration refinement grids, peaking at about 1.5 GB, and the superseded held-out run. The 3 GB runner reserve was not enough with two Level 2 jobs plus that load.
  - `L2_B3LYP_TDA15` had not reached its first TD checkpoint after 2 h: under CPU contention one 5-iteration chunk took longer than that. All of its work was lost.
  - At 12:36 UTC `L1_FULL` failed with "No space left on device". PySCF writes about 3.8 GB of density-fitting tensors per Level 2 run to `/tmp`, and the killed runs had left five orphaned files (19 GB).
- **Repair:**
  - The orphans (verified not open by any process) were deleted; 19 GB is free again.
  - Each job now gets its own scratch directory (`TMPDIR`/`PYSCF_TMPDIR` = `<job>/scratch`), which is emptied before every (re)launch.
  - TD solves checkpoint every 2 Davidson iterations (`TD_CHUNK = 2`).
  - My analysis processes run under an address-space limit (`prlimit --as` 3.5 GB on the held-out run) at `nice` 15.
  - The runner was restarted. The running CAM-TDA15 and FE_CAT jobs were not touched; TDA15 and L2_OPT are marked interrupted and requeued, and L1_FULL is retried (attempt 2).
- **Test:** `phase1/tests/test_jobs.py::test_scratch_is_per_job_and_cleaned_before_relaunch` (3 passed in total).
- **Status:** FIXED_AND_VERIFIED (cause and policy). The lost production runs: IMPLEMENTED_AWAITING_PRODUCTION_RUN.

## Held-out result (pre-registered evaluation; `records/heldout_v2/`)

Run: `phase4/heldout.py` at commit d57264e and later, with the frozen primary (`docs/preregistration_heldout.md`) and observation SD 0.10 dex for S Greenland. Tables were computed 12:26–14:29 UTC and scoring finished after the container restart at 15:33 UTC.

### Pre-registered contrast: `tddft_D` − `tierA_empirical`, mean log predictive score of the counts
| fold | difference | 95 % day-block bootstrap | n samples / days |
|---|---|---|---|
| primary (S6 2017 → S Greenland 2021) | **+0.153** | 0.117 to 0.168 | 18 / **2** |
| secondary (S Greenland → S6) | **+0.537** | 0.360 to 0.725 | 46 / 9 |

**Against the pre-registered rule:**
- The rule is *formally* met: the primary-fold interval excludes 0 and the secondary fold has the same sign.
- **The primary-fold interval is not a credible uncertainty.** With only 2 sampling days, the day-block bootstrap can only return the two day means (0.117 and 0.168) and their average, so the "interval" is simply their range. On a per-sample basis, `tddft_D` scores better than Tier A on 17 of 18 S Greenland samples. The secondary fold (9 days) has a usable interval.
- **Classification:** the pre-registered improvement in held-out abundance prediction over Tier A empirical optics is supported in direction at both sites, and robustly in the secondary fold. The primary fold's site-level uncertainty cannot be quantified with 2 sampling days. This is a modest improvement in log score (0.15–0.54 nats per sample), not a demonstration of better melt prediction.

### Descriptive results (no winner selected; not confirmatory)
- **Primary fold, log score (higher is better):**

  | model | log score |
  |---|---|
  | measured_mac_C | −0.81 |
  | ridge | −0.83 |
  | tddft_D_iid | −0.88 |
  | tddft_D | −0.96 |
  | tierA | −1.11 |
  | tddft_C | −1.23 |
  | climatology | −1.34 |
  | literature prior | −1.43 |
  | band ratio | −1.65 |

  RMSE is 0.43–0.45 dex for measured MAC, tddft_D and tddft_D_iid, 0.39 for Tier A and 0.50 for ridge.
- **Secondary fold:** tddft_C (−6.11), measured_mac_C (−6.23), tddft_D_iid (−6.27), literature prior (−6.30) and tddft_D (−6.38) are close together; Tier A scores −6.91 and ridge −7.83.
- **The primary is not the best-scoring physics variant in either fold.** Measured in vivo MACs and the independent-residual calibration score slightly better. Per the pre-registration, these are reported, not adopted.
- **Coverage:** 95 % coverage is 1.00 for the physics models in the primary fold, so the predictive intervals are wide; the training-site discrepancy τ is 0.72 dex for tddft_D and 1.06 for Tier A. They are conservative rather than calibrated.
- **Forward test (4-band HCRF given the measured abundance):**
  - Primary fold: tddft_D vs Tier A −0.18 (block range −1.37 to 0.28, i.e. no difference). Secondary: +2.30 (0.87 to 3.83).
  - The training-site statistical baselines (band climatology, regression on log B) predict the bands far better than any physics model (by about 6 nats).
  - The physics models' advantage is in inverting reflectance for abundance, not in forward spectral fidelity.
- **Broadband HCRF:** all physics models over-predict the measured broadband HCRF. The bias is +0.23 (tddft_D), +0.19 (Tier A) and +0.29 (tddft_C) in S Greenland, and +0.10 to +0.26 at S6. This is a systematic broadband/NIR brightness (or k) mismatch that the retrieval absorbs through k and σ. **UNRESOLVED.**
- **Scope:** no statement about glacier-melt prediction follows from this evaluation.

### P5-CORR-1: corrected interpretation of the frozen held-out result (`records/heldout_v2/`, now read-only with `FROZEN_PROVENANCE.json`)
The previous classification, "improvement supported in direction at both sites", overstated the evidence. Corrected:
- **Log score vs RMSE.** On the primary fold, `tddft_D` has a better log predictive score than Tier A (+0.153 per sample) but a WORSE abundance RMSE: 0.447 vs 0.394 dex. CRPS is −0.012 better for the primary, with a 2-day block range spanning 0. The log-score gain comes from a better-calibrated predictive spread (smaller τ: 0.72 vs 1.06 dex), not from more accurate point retrievals.
- **Literature prior.** On the reverse fold, `tddft_D` does NOT beat the literature prior (N(3.56, 0.78) from S6 2016) in log score: −6.376 vs −6.302, difference −0.074, 9-day block interval −0.385 to +0.050. On that fold a fixed prior from independent S6 data is as good as the retrieval.
- **Interval.** Two primary sampling days cannot support a generalisation confidence interval. See `docs/preregistration_deviations.md` D2.
- **Optics.** The experiment used posterior-mean optics, not joint calibration draws (`docs/preregistration_deviations.md` D1).
- **Scope.** The experiment tests abundance retrieval and forward reflectance at field-plot scale. It does not test albedo, absorbed shortwave or melt prediction.
- **Status:** classification corrected to **"formal rule met; evidence of improved abundance prediction is weak (log-score only, not RMSE; not better than an independent prior on the reverse fold; primary uncertainty not estimable)"**.

### P5-RES-1: shared cgroup-aware resource budget and resumable scoring
- **Problem.** Admission used host `MemTotal` (16 GB) instead of the container's memory cgroup limit (14.35 GB). Analysis processes were not budgeted at all, and that caused the OOM kills P1-ENV-2/3.
- **Fix.**
  - `common/resources.py` reads the cgroup v1/v2 limit and unreclaimable usage, and keeps a file-locked reservation ledger (`phase1/results/resources/budget.json`). It admits on memory, threads and disk.
  - `common/run_budgeted.py` provides budgeted launches with thread-pool caps, an `RLIMIT_AS` backstop and per-run scratch.
  - The `phase1/jobs.py` runner now admits every job through the same ledger.
  - `phase4/heldout.py::cached_baselines` stages baseline fits per fold, so scoring resumes at fold boundaries.
- **Tests:**
  - `python -m pytest common/tests/test_resources.py`: 7 passed.
  - `phase1/tests/test_jobs.py`: 4 passed.
  - `phase4/tests/test_heldout.py::test_baseline_scoring_resumes_from_fold_stage`: passed.
  - All use fixtures; none pushes the container to its limit.
- **Recovery guarantees** per phase are in `docs/compute_reliability.md`.
- **Limitation.** Thread reservations made before the runner restart (on-going jobs registered after the fact) can temporarily exceed the core count. New work is held back until that clears; nothing is killed.

### P5-DIAG-0: deviation D3 found while preparing the forward diagnostics
- **Finding.** The k prior (`ARF_master`) equals HCRF / albedo of 51 S6 2017 plots, 29 of them counted held-out plots. Verified to 4 decimals.
- **Recorded as:** `docs/preregistration_deviations.md` D3; a protocol erratum (appended, frozen text unchanged); the data manifest; and the `h1_albedo.py` docstring.
- **Impact.** H1 does not use k. The reverse fold of `records/heldout_v2` is mildly optimistic for all physics models alike. Frozen results are not re-scored.

### P5-DIAG-1: controlled forward-optics diagnostic (`phase4/forward_diagnostics.py`)
- **Design.** Direct BioSNICAR on the S6 2017 plots with measured counts. Nuisances are fixed a priori at the literature-prior centre and varied one at a time across groups A–D.
  - Exact split of the HCRF error into a geometry term (k_prior − k_obs)·A_model and an albedo term k_obs·(A_model − A_obs).
  - Residual spectra.
  - A prior ensemble (40 common random draws) gives bias, z, PIT and coverage.
- **Not representable**, listed in the output: water films, roughness, snow patches, dust profiles, non-diffuse cloud spectra, separate packaging/Fe ablations.
- **Run.** Queued under the budget (smoke test first). Results go to `docs/forward_diagnostics.md` when complete.

### P5-COUPLE-1: host coupling (`model_integration/`)
- **Reconciled from** `codex/glacier-model-integration` (0fc0096). The CI workflow file was deliberately excluded.
- **Fix.** The adapter's film+dust refusal was stale: on this branch dust stays uniform. The README was updated.
- **New `host_modes.py`:**
  - `HostContract` requires the host's algae treatment.
  - Mode A (forward, prescribed abundance) refuses observed-albedo hosts, and explicit-algae hosts unless their scheme is switched off and `replace` is used.
  - Mode B (state estimation) keeps the observation, gives a signed counterfactual, and flags every output "not an independent validation".
  - `reference_seb_paired` runs `phase4/seb.py` with distinct output labels.
- **Tests:** `python -m unittest discover -s model_integration/tests`: 23 passed.

### P5-SCENE-1: scene interpretation product (`phase4/scene_product.py`, schema `algae-scene-product/1.0.0`)
- **Classes:** cloud/shadow, water, snow, non-ice, model-inconsistent, prior-dominated, supported, no meaningful darkening, ambiguous.
  - "Supported" = P(Δα ≥ 0.01) ≥ 0.9 AND Bayes factor ≥ 10 against no-algae. Thresholds fixed before any scene output.
  - Dark/unclassified SCL pixels are kept and flagged. No species fractions are reported.
- **Exports:** COG ×2, NetCDF4 CF-1.8, provenance sidecar and report.
- **Validation every run:** checksums, COG layout, georeference, exact value round trip, coordinates, masks and ranges.
- **Tests:** `phase4/tests/test_scene_product.py`: 4 passed (synthetic emulator, tamper detection, immutability, out-of-domain SZA refusal).
- **Scene run:** queued under the budget (S2A_22WEV_20190723, 3 km).
- **Docs:** `docs/scene_product.md`, `docs/schema/scene_product_v1.json`, `docs/limitations_appropriate_use.md`, `docs/evidence_table.md`. Dependencies are in `requirements.txt` (netCDF4 added).

### P5-SURR-1: surrogate qualification binding and domain rejection (`phase4/surrogate_api.py`)
- **Binding.** Benchmark records are bound to `benchmark_fingerprint()` (surrogate code, tolerances, design version) in addition to the emulator content tags. Stale or unbound records are rejected (tested).
- **Domain.** Out-of-domain requests raise `DomainError` by default; `strict=False` clips and flags. Zero algae is not a domain point.
- **Tolerances.** Δα tolerance 0.005 added, and the tolerances documented as error-budget based (half of 0.01). Set before any benchmark result.
- **Contract.** Units, concentration convention, illumination, 300–2500 nm support and ranges are in the module docstring.
- **Tests:** `phase4/tests/test_surrogate.py`: 3 passed.
- **Re-benchmark.** The benchmark started earlier (old code) cannot qualify by design. A re-benchmark on the cached emulators is queued after it.
- **Not yet benchmarked:** SZA > 60° (hourly coupling at low sun needs nodes at 70 and 80°). Such requests are refused.
