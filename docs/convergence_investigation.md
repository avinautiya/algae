# Bounded convergence investigation of the restartable TD-DFT solves (P9-CONV)

Started 2026-10-20 14:45 UTC. The production L1_FULL run was **paused** (runner and solver stopped). All checkpoints, SCF chk files and logs are preserved and copied to `phase1/results/ckpt_backup_20261020_L1_before_diag/`; diagnostics work on copies only (`phase1/convergence_diag.py`, outputs in `phase1/results/convergence_diag/`). No tolerance, root count, functional, basis, solvent or TD/TDA setting has been changed in production. The unconverged L1 spectra remain excluded from every scientific conclusion.

## 1. Numerical configuration actually in force

Source: `phase1/results/convergence_diag/config_L1_B3LYP.json` and `config_L1_CAM-B3LYP.json` (`python3 phase1/convergence_diag.py config`). The solver was inspected in the installed package, not from online documentation.

| Item | Value |
|---|---|
| Code commit at start of investigation | a302c4f |
| PySCF / NumPy / Python | 2.14.0 / 2.5.3 / 3.13.16 |
| Response type | **Full TD-DFT (RPA)**, singlets, non-equilibrium linear-response IEF-PCM (ε 78.4; optical ε 1.78), Lebedev 29 |
| Object | `PCMTDDFT → TDDFT → TDHF` |
| Eigensolver | `pyscf.tdscf._lr_eig.real_eig` (non-Hermitian, real), file sha256 96186601… |
| Functionals | B3LYP (`B3LYPG`); CAM-B3LYP (`CAMB3LYP`) |
| Basis | 6-31G* with Cartesian d (256 AO); DF def2-universal-jkfit; DFT grid level 3 |
| Molecule | Level 1 chromophore, 24 atoms, 114 electrons; 57 occupied, 199 virtual |
| Response space | A_size = 57 × 199 = 11 343 per half (vector length 22 686) |
| Requested roots | 30 |
| Convergence criterion | per-root residual norm ≤ tol, vectors normalised to X·X − Y·Y = 1; staged 10⁻² → 10⁻³ → 10⁻⁴ → 10⁻⁵ |
| Other solver settings | lindep 10⁻¹², level shift 0, deg_eia_thresh 10⁻³ |
| Search space per `real_eig` call | space_inc = 30 (max(nroots, min(MAX_SPACE_INC = 20, A/2))); **max_space = 2862** vectors from td.max_memory = 2100 MB (four holders of A × max_space ≈ 1.04 GB). An uninterrupted solve adds at most 30 vectors per cycle and collapses to the Ritz vectors only when the space is full (after about 95 cycles). |
| Production chunking | 2 cycles per `td.kernel` call (`--td-chunk 2`); budget 400 nominal cycles per functional |
| Production resources | 3000 MB reservation, **2 threads reserved, but 4 used**: `qc.set_threads()` overrode the budget's `OMP_NUM_THREADS` with `os.cpu_count()` (fixed in this investigation; affects timing, not numerics); container cgroup limit 13.68 GB, 4 cores |
| Ground state | SCF conv_tol 10⁻⁹ (energy); the production log shows \|g\| about 3×10⁻⁵ at convergence for CAM-B3LYP (PySCF default gradient threshold √conv_tol). A first diagnostic value of 0.044 was an artefact: PySCF's default `get_grad` omits the PCM potential for solvated objects; corrected in section 3. |

## 2. What the checkpoint scheme preserves, and what it discards

**Code:** `phase1/tdcheckpoint.py` (`solve`, `run_chunk`, `start_vectors`) driving `pyscf.tdscf.rhf.TDHF.kernel`, which for full TD-DFT with a hybrid functional calls `pyscf.tdscf._lr_eig.real_eig` (installed PySCF 2.14.0, inspected in place).

**Saved after every chunk** (`td_ckpt/<F>_RPA_30/latest.npz`, previous copy in `latest.prev.npz`):
- the 30 Ritz vectors `[X, Y]` (from `td.xy`, i.e. renormalised to X·X − Y·Y = ½) and their Ritz values;
- per-root converged flags at the current stage tolerance;
- the MO coefficients, occupations and energies they are expressed in;
- metadata: operator fingerprint hash, nominal cycle count, completed-stage tolerance, stage history.

**Not saved, and therefore discarded at every chunk boundary.** These are the entire internal state of `real_eig`:
- the trial basis `V, W` (up to `max_space` vectors per half; see section 1) and the operator products `U1 = A V + B W`, `U2 = A W + B V`;
- the projected matrices `a, b, σ, π` of the subspace problem;
- the correction (preconditioned residual) vectors generated in the **last** cycle of a chunk. They are computed and never used: with `max_cycle = 2`, cycle 0 builds the space from the start set, cycle 1 adds one block of corrections, and the corrections from cycle 1 are dropped when the loop ends;
- the previous Ritz values/vectors that `real_eig` uses to match roots between iterations (`_sort_elast`).

**What a 2-cycle chunk therefore is.** Each `td.kernel(x0=...)` call is a **fresh start**.
1. **Start set:** `start_vectors()` builds it from the 30 saved Ritz vectors plus the solver's 30 standard initial-guess vectors (lowest orbital-energy-difference excitations), orthonormalised by SVD in the Euclidean [X, Y] sense.
2. **Cycle 0:** the response operator is applied to all of them (about 60 vectors).
3. **Cycle 1:** the operator is applied to one block of corrections (one per unconverged root), and the 30 lowest Ritz pairs of that space are returned.

So every chunk performs **one subspace expansion** at the cost of about 60 + n_unconverged operator applications. An uninterrupted solve instead performs one expansion per cycle at a cost of about n_unconverged applications, keeping all previous directions until the space reaches `max_space`. Only at that point does `real_eig` itself collapse to the 30 current Ritz vectors (`fresh_start`).

**Consequences:**
- The production scheme is a different iterative method: a restarted block-Davidson with a restart length of one expansion. It is not a checkpointed copy of one uninterrupted solve, and "400 cycles" in the 2-cycle scheme is about 200 expansions, each from a nearly empty space.
- PySCF treats `x0` as an initial subspace only (`real_eig` builds `V, W` from it on a fresh start). Nothing in the installed code restores a previous search space.
- Restarted methods with short restart lengths are known to converge slowly for clustered or high-lying roots, but that is a hypothesis for this system and is tested in section 4. The 30 re-added guess vectors are the same at every chunk; after the first chunk they mostly re-cost operator applications without adding new directions.
- The nominal "cycles" counter therefore understates cost differences. Section 4 compares operator applications and wall time.

## 3. Per-root evidence

(pending)

## 4. Restart vs uninterrupted segments

(pending)

## 5. Remedies

(pending)

## 6. Terminal-job handling

(pending)

## 7. Stopping rule and resources

(pending)
