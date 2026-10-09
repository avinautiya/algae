# Phase 1 jobs: how to run, resume and check them

All Phase 1 production calculations are defined in `phase1/jobs.py` and write to `phase1/results/v2/<JOB>/`.

```bash
python phase1/jobs.py status                      # state + validation of every job
python phase1/jobs.py run --max-jobs 2            # start/resume everything runnable (safe to re-run)
python phase1/jobs.py run --retry-failed          # also retry jobs that failed MAX_ATTEMPTS times
python phase1/jobs.py run --only L2_CAM_TDA15     # one job
python phase1/compare_states.py                   # approximation checks A/B/R/G/P/T once outputs exist
```

- **Resuming:** excited-state solves checkpoint every `--td-chunk` Davidson iterations under `<outdir>/td_ckpt/`, and optimisations resume from `opt_*_last.xyz`. After a container restart, run `jobs.py run` again; at most one chunk is lost.
- **Completion:** a job is complete only when `completion.validate_run` accepts its `COMPLETE.json` (checksummed artifacts). A `DIAGNOSTIC.json` marks unconverged or invalid results, which are never consumed downstream.
- **Memory:** jobs declare a peak memory (`mem_mb`). The runner keeps `MEM_RESERVE_MB` (3 GB) free for other work, so two 7 GB Level 2 jobs do not run at the same time on a 16 GB host.
- **Not scheduled here:** full TD-DFT CAM-B3LYP at Level 2 needs a host with ≥ 4 cores and ≥ 12 GB that stays up for about 2 days (command in `docs/repair_ledger.md`, P1-ENV-1).
- **Restart caveat:** a runner that is already running does not see jobs added to `JOBS` after it started (for example `L2_COO_OPT_TDA15`). Restart it with `jobs.py run`; running jobs are left alone, and their state is detected from their PID and output markers.
