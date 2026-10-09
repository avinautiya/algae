# Compute reliability: budget, recovery guarantees, scratch

## Shared resource budget (`common/resources.py`, `common/run_budgeted.py`)

- **Memory limit:** the effective limit is the memory cgroup of this container (cgroup v1 `.../claude-code-bash`, 14.35 GB = 13 680 MiB), read along the cgroup path. Host `MemTotal` (16 GB) is NOT the limit.
- **Unreclaimable usage:** anonymous RSS + mapped file + shmem. Unmapped page cache is reclaimed before an OOM kill.
- **Admission:**
  - Condition: unregistered usage + Σ max(reservation, current tree RSS) + request ≤ limit − 1200 MB; reserved threads + request ≤ cores; free disk ≥ scratch + 5 GB.
  - Every Phase 1 job (`phase1/jobs.py`) and every analysis launched through `common/run_budgeted.py` holds a reservation in the shared, file-locked ledger `phase1/results/resources/budget.json`.
  - Reservations of dead processes are pruned automatically.
- **Thread pools:** `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, `NUMEXPR_NUM_THREADS`, `VECLIB_MAXIMUM_THREADS` and `BLIS_NUM_THREADS` are all set to the reserved thread count.
- **Hard backstop:** `run_budgeted` sets `RLIMIT_AS` (default 2.5 × reserved memory), so a runaway analysis fails with MemoryError instead of pushing the cgroup into an OOM kill of chemistry. `nice` is not used as a memory control.
- **Peak-memory estimates** (`resources.ESTIMATES_MB`, conservative, from observed runs):

  | work | estimate |
  |---|---|
  | Level 2 TD | 7.5 GB |
  | Level 2 optimisation | 6 GB |
  | Level 1 | 3 GB |
  | Fe complexes | 4.5 GB |
  | AR(1) calibration | 1.5 GB |
  | emulator build (2 workers) | 1.5 GB |
  | held-out scoring | 2 GB |
  | scene product | 2.5 GB |

- **Tests:** `common/tests/test_resources.py` (7 passed) covers cgroup v1 and v2 limits on fake `/proc` and cgroup trees, over-commit refusal, stale-reservation pruning, disk admission, thread caps, and the budgeted runner (refusal without running, thread and scratch environment, ledger release, `RLIMIT_AS` failure of the offending process only). `phase1/tests/test_jobs.py::test_runner_waits_when_budget_refuses` checks that the runner never launches a refused job. No test pushes the real container to its limit.

## Per-job scratch

- Each Phase 1 job has `TMPDIR`/`PYSCF_TMPDIR` = `phase1/results/v2/<job>/scratch`, emptied when the job is (re)launched. That is the only time a job can be inactive under the runner's control.
- Each budgeted analysis has `phase1/results/resources/scratch/<name>`, removed after a successful run and kept after a failure.
- Checkpoints (`td_ckpt/`, `scf_*.chk`, `opt_*_last.xyz`, `*_steps.jsonl`) live outside scratch and are never cleaned automatically.

## Recovery guarantees after a kill or container restart (state them separately; none is "at most one chunk" overall)

| Phase | What survives | Maximum loss |
|---|---|---|
| Setup (molecule, PCM, grids) | nothing (recomputed) | minutes |
| Density-fitting integrals (scratch, ~3.8 GB at Level 2) | nothing (scratch is recreated) | the full DF build (tens of minutes under contention) |
| SCF | `scf_<functional>.chk` (orbitals): the restart begins from the converged/last orbitals | the SCF iterations since the last chk write; a converged SCF restarts in 1–2 cycles |
| Geometry optimisation | `opt_<F>_last.xyz` and `_steps.jsonl` after every gradient step; resumes from the last geometry | the current step (SCF + gradient) |
| TD/TDA Davidson | `td_ckpt/.../latest.npz` every `TD_CHUNK` = 2 iterations (validated by operator fingerprint and orbital mapping) | the setup + DF + SCF time above PLUS up to 2 iterations. **Before the first checkpoint, everything since launch is lost.** Under contention a Level 2 iteration has taken > 20 min. |
| Held-out scoring | physics tables (`cache/tables_*.npz`, fingerprinted) and per-fold baseline stages (`stages/baselines_<fold>.pkl`) | at most one fold's baseline fits |
| Phase 4 scene runs | `COMPLETE.json` written last; incomplete runs are redone | one scene |

## Incidents behind this design

P1-ENV-2 and P1-ENV-3 in `docs/repair_ledger.md`: two OOM kills of chemistry by unbudgeted analysis, and disk exhaustion from orphaned DF tensors in `/tmp`.
