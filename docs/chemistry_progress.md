# Durable chemistry progress (audited 2026-10-09 21:35 UTC from checkpoint metadata and job logs)

"Runner restarted" is not progress. **Durable progress** = a new validated checkpoint (`td_ckpt/*/latest.npz`: total Davidson cycles, roots converged at the current stage, completed stage tolerance) or a validated completion marker.

## Per job

| | L2_CAM_TDA15 (Level 2, CAM-B3LYP, TDA, 15 roots, tol 1e-5) | L1_FULL (Level 1, B3LYP + CAM-B3LYP, full TD, 30 roots, tol 1e-5) |
|---|---|---|
| Last valid checkpoint | `latest.npz` written 17:07:06 UTC; fingerprint 8d13072a0853 | `latest.npz` 18:12:05 (prev 18:00:23); fingerprint 303d4c157142 (B3LYP_RPA_30) |
| Saved state | 2 Davidson cycles; 6/15 roots at residual 0.01; no stage completed | 18 cycles; 24/30 roots at 0.01; no stage completed; CAM-B3LYP not started |
| Launches (UTC) | 10:58, 12:19, 14:52, 19:52, 20:57 | 12:19, 14:52, 19:52, 20:57 |
| Setup + DF + SCF before TD (last productive launch) | 44 min (14:52 → checkpoint lock 15:36:35; ground-state B3LYP SCF and CAM SCF) | optimisation completed in the 14:52 launch (final geometry 16:27), TD from 16:27 |
| TD chunk wall time (2 cycles) | 5431 s (90.5 min), under CPU contention | 491–790 s (8–13 min) |
| Time from a relaunch to the next durable checkpoint | about 134 min | before the fix: SCF + gradient + geomeTRIC re-check of an already-converged geometry, then about 10 min of TD. **Neither relaunch reached TD.** After the fix (`--skip-opt` from the converged geometry): SCF (minutes) + one chunk (about 10 min). |
| Advancement across the 19:52 and 20:57 relaunches | **none** (no checkpoint written) | **none** (no checkpoint written) |
| After the fixes (21:53 relaunch) | not launched (cannot reach a checkpoint here) | SCF, then the 18-cycle checkpoint was reused via the rounding-tolerant match (no stage credit). **New durable checkpoint 22:01:14 UTC: 20 cycles, 25/30 roots at the current stage; fingerprint 32a898f0c632 (now the reference for further restarts).** Further checkpoints followed; latest 22:25 UTC: 28 cycles, 25/30 roots at the current stage (the count fluctuates as Davidson refines). No stage is completed yet. The container rebooted at about 23:02 (uptime about 92 min); at 23:03 the runner detected the stale record through the identity check and relaunched from the 22:25 checkpoint. |

## Environment

- **Observed uptimes between reboots:** 19:52 → about 20:56 (about 64 min); 20:57 → about 21:29 (about 32 min); 21:30 → about 23:02 (about 92 min).
- **L2_CAM_TDA15 cannot complete here.** Time to its next checkpoint (about 134 min) exceeds the runtime between reboots. At about 90 min per 2-cycle chunk, the remaining stages (residual 0.01 → 10⁻⁵, 15 roots) are estimated at 15–30 chunks, i.e. **22–45 h of uninterrupted compute** (rough; the cycle count is not known in advance).
- **L1_FULL can advance** once re-optimisation is skipped: about 10 min to the first checkpoint after SCF.
  - Remaining: B3LYP 30 roots to 10⁻⁵ (estimated 20–40 more chunks, i.e. 3–8 h), then CAM-B3LYP 30 roots from scratch (similar).
  - These are estimates, not measurements.

## Defects found and fixed (P8-REL-1)

1. **L1_FULL re-ran a converged geometry optimisation on every relaunch.** `jobs.command` now passes `--skip-opt --start-xyz opt_B3LYP_final.xyz` when the job's own optimisation record says converged. It writes the acceptance sidecar (`.geometry.json`, with thresholds, energy and xyz sha256) required by `completion.geometry_status`.
2. **The checkpoint fingerprint was not reproducible from the saved geometry.** Re-reading the `.xyz` flips the 6th decimal (Bohr) of one coordinate (difference 1.0×10⁻⁶ Bohr). The next TD start would have overwritten `fingerprint.json`, rejected the 18-cycle checkpoint as "a different operator" and started fresh.
   - **Now:** the stored fingerprint is read before it is overwritten. If every non-geometry field is identical and the geometry agrees within 10⁻⁵ Bohr, the restart vectors are reused without stage credit. If the orbitals also differ marginally, the vectors are used only as an unvalidated initial guess, and the solve reconverges on the current operator to the full tolerance.
   - The benefit (fewer cycles) is **not demonstrated**: on the formaldehyde test system the guess saved no cycles.
3. **Process adoption used PID existence only.** After a reboot, PIDs are reused, so a stale "running" record could block a relaunch forever.
   - **Now:** jobs and reservations record boot ID, process start time and command-line hash. Adoption requires all to match, and the command line must name the job's output directory.
   - **Legacy records** without identity are accepted only if the live process started before the recorded launch.

## Not changed

Root counts, TD vs TDA, residual tolerances, chunk size (2 cycles), basis, functional, solvent.

**One-iteration chunks** were NOT adopted:
- each chunk restarts Davidson from Ritz vectors, which can increase total work (`docs/compute_reliability.md`);
- it would cut the Level 2 lost-work window from about 90 to about 45 min, still above the 32-min window observed;
- it would need a matched small-system convergence/timing test before use.

## Recommendation

Run L2_CAM_TDA15 (and the queued Level 2 / Fe jobs) on a stable machine or scheduler: `docs/portable_chemistry_package.md`. No paid compute has been provisioned.
