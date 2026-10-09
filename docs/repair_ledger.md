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
