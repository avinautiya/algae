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

**Method.**
- **Tool:** `convergence_diag.py residuals` (`phase1/results/convergence_diag/res_B3LYP/residuals.json`, `res_CAM/residuals.json`).
- **Inputs:** every saved checkpoint copy is mapped into freshly converged orbitals (rotation from the AO overlap; orbital-energy differences ≤ 1×10⁻⁵ Eh; occupied/virtual subspaces unchanged to 10⁻⁸).
- **Residuals:** one operator application gives the true residual of each saved Ritz pair in the solver's own convention.
- **State tracking:** by RPA-metric amplitude overlap between consecutive checkpoints, plus principal angles of near-degenerate clusters (gaps < 0.05 eV).
- **Checkpoints available:**
  - B3LYP: cycles 16, 18, 26, 28 (backups), 400 and 402;
  - CAM-B3LYP: cycles 58 and 60.
  - The production logs record only converged-root counts, not residuals.

**Ground state.**
- Orbital gradients with the solvated Fock matrix: B3LYP 2.6×10⁻⁷, CAM-B3LYP 4.0×10⁻⁶.
- The saved Ritz values agree with the Rayleigh quotients on the current operator to ≤ 9×10⁻⁵ eV.
- So the stored states belong to the operator they are tested against. Internal stability: see the end of this section.

**B3LYP (exhausted at 400 cycles), true residuals:**

| Root (index) | E (eV) at 402 | Residual at cycles 16 / 18 / 26 / 28 / 400 / 402 |
|---|---|---|
| 0–28 (29 roots) | 2.861 … 7.629 | all 4×10⁻³ – 1.0×10⁻² from cycle 400; most were already 6–9×10⁻³ at cycle 16 and moved by ≤ 10 % over 380 cycles |
| 21 | 7.169 | 1.6 / 1.7 / 1.7 / 1.5×10⁻² / 9.8×10⁻³ / 9.8×10⁻³ |
| 25 | 7.466 | 1.2 / 1.2 / 1.2 / 1.1×10⁻² / 7.9×10⁻³ / 7.9×10⁻³ |
| 27 | 7.580 | 1.6 / 1.7 / 1.7 / 1.5×10⁻² / 9.8×10⁻³ / 9.8×10⁻³ |
| **29 (root 30)** | **7.674** | **3.0 / 2.2 / 2.1 / 2.4 / 1.6 / 1.6 ×10⁻²** |

- **Classification:** root 30 is improving very slowly, not oscillating and not being exchanged.
  - Its Ritz value decreases monotonically (7.697 → 7.674 eV) and its residual falls by a factor of about 1.5 over about 186 restarted expansions (cycles 28 → 400).
  - It keeps its index (overlap with the previous checkpoint > 0.9 at every step) and stays in a stable top cluster (roots 28–30; subspace cosines 0.97–1.0).
  - In practice this is stagnation under the current scheme, not intrinsic non-convergence.
- **Roots below tolerance are frozen:** once a root's residual falls below the stage tolerance, `real_eig` stops expanding for it. Their residuals then stay just under 0.01 for hundreds of cycles; they are "converged" only in the stage sense.
- **Near-degenerate cluster:** roots 19–21 (7.135–7.169 eV, two within 2 meV) rotated within their 3-dimensional subspace between cycles 28 and 400 (smallest principal cosine 0.63), without index exchange. Any per-state property in this cluster is ill-defined until it converges tightly; the cluster subspace is the well-defined object.

**CAM-B3LYP (cycle 60), true residuals:**
- 12 roots are above 0.01: indices 15, 18, 20–29 (7.37–8.61 eV; maximum 3.9×10⁻² at root 27 / index 26).
- Between cycles 58 and 60, 10 of these 12 improved by 5–25 % (e.g. index 22: 2.3 → 1.8×10⁻²; index 25: 1.7 → 1.3×10⁻²).
- **So the flat "18/30 converged" count since cycle 30 hides steady per-root progress.** CAM-B3LYP is improving slowly, not stalled.
- There are no index exchanges; clusters [3,4], [17,18], [20–22] and [25–27] are stable (cosines ≥ 0.998).

**Converged-root count is not a diagnostic.** The same count concealed both a single slowly creeping root (B3LYP) and broad, steady improvement (CAM-B3LYP).

## 4. Restart vs uninterrupted segments

**Design.**
- **Code path:** both arms run the production code path (`tdcheckpoint.solve`) on COPIES of the same checkpoint, with the same physics (full TD-DFT, PCM, functional, basis, 30 roots), the same staged tolerances (10⁻² → 10⁻⁵) and the same start set (Ritz vectors + standard guess).
- **The only difference is the chunk length:**
  - "restart": 2 cycles per `td.kernel` call, as in production;
  - "uninterrupted": up to 60 (small system: 200) cycles per call.
- **Resources:** 2 threads each, run concurrently on the 4-core container so the wall-clock comparison is fair.
- **Measures:** operator applications (vectors passed through the response operator), operator time and wall time; residuals from the solver's per-root DEBUG1 log.
- **Records:** `phase1/results/convergence_diag/{small,cam,b3}_*/` (`solver.log`, `calls.jsonl`, `segment.json`), summaries in `records/convergence_investigation/`.

**Small representative system** (same molecule and settings, 3-21G basis; fresh start; B3LYP):

| | Cycles | Operator vectors | Wall/operator time | Result |
|---|---|---|---|---|
| Uninterrupted | 148 | 945 | 3800 s | all 30 roots < 10⁻⁵ (status converged). All < 10⁻² after 207 vectors (864 s), < 10⁻³ after 459 (1884 s). |
| 2-cycle restart (stopped at equal wall time) | 30 | about 1060 | about 3830 s | 5 roots still > 10⁻² (max 1.4×10⁻²); first stage not reached |

**Stalled system: L1 CAM-B3LYP, 6-31G*, from the production cycle-60 checkpoint** (12 roots > 10⁻², max 3.9×10⁻²). Both arms ran 2 h wall:

| | Cycles | Operator vectors | Operator time | Roots > 10⁻² | Roots > 10⁻³ | Max residual |
|---|---|---|---|---|---|---|
| 2-cycle restart (production scheme) | 16 | 583 | 6870 s | 12 | 30 | 3.85×10⁻² (unchanged) |
| Uninterrupted | 77 | 514 | 6491 s | **0** (stage 10⁻² done at cycle 39, 220 vectors) | **0** (stage 10⁻³ done) | 9.1×10⁻⁴ (in stage 10⁻⁴) |

**Conclusion (restart-induced stagnation, confirmed on the stalled system):**
- At equal or lower operator cost, the uninterrupted solve converges where the 2-cycle scheme makes no net progress.
- The cause is the restart length. Each 2-cycle call spends about 61 operator applications re-applying the operator to its start set for a single expansion, and discards every search direction.
- The physics, tolerance and root count are not the cause: the same solver on the same operator converges once it is allowed to keep its subspace.
- **Correction to an interim remark:** an early look at maxima over every other cycle suggested that each restart actively increases residuals. A per-root check of the small-system run found no root whose residual rose by more than 5 % across any restart boundary. The scheme mainly loses progress; it does not degrade the solution.

## 5. Remedies

(pending)

## 6. Terminal-job handling

(pending)

## 7. Stopping rule and resources

(pending)
