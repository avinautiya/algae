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
| After the fixes (21:53 relaunch) | not launched (cannot reach a checkpoint here) | SCF, then the 18-cycle checkpoint was reused via the rounding-tolerant match (no stage credit). **New durable checkpoint 22:01:14 UTC: 20 cycles, 25/30 roots at the current stage; fingerprint 32a898f0c632 (now the reference for further restarts).** Further checkpoints followed; latest 22:25 UTC: 28 cycles, 25/30 roots at the current stage (the count fluctuates as Davidson refines). No stage is completed yet. The container rebooted at about 23:02 (uptime about 92 min); at 23:03 the runner detected the stale record through the identity check and relaunched. **That 23:05 relaunch discarded the 28-cycle checkpoint** (defect 4 below); the reboot at about 00:05 struck before a fresh run could overwrite the file. After the fix, the 00:09 relaunch reused the 28-cycle vectors as an unvalidated initial guess (no stage credit). **The guess saved work:** the first chunk after it stood at 25/30 roots (cycle 30), whereas a fresh start reached only 19/30 after 8 cycles (about 90 min). Checkpoints then came every about 7 min; by 01:01 the job was at 42 cycles, 26/30 roots; at 01:34, 54 cycles, 27/30. The container rebooted again by 02:10 and the runner was relaunched at 02:11; one more checkpoint at 02:15 (56 cycles, 27/30), then the next restart before 03:14. The 03:14 relaunch was reclaimed before its first chunk. At 04:16 the job was relaunched with a pending background wait keeping the session active: 58 → 64 cycles by 04:33 (about 4 min per chunk), 27/30 roots, stage not completed. Continuous running to 06:30: **126 cycles, 29/30 roots** at the residual-0.01 stage; one root has not converged since about cycle 86. No stage completed. At 08:27: 184 cycles, still 29/30. |

## Environment

- **Correction to "random reboots" (03:14 UTC).** The restarts coincide with the Claude Code session going idle; the cloud container is reclaimed after a period of inactivity.
  - The last checkpoint before each restart came a few minutes after the session's last activity: 01:34 after a turn ending about 01:31; 02:15 after one ending about 02:12.
  - Each restart was found when a scheduled wake-up reactivated the session (00:05, 02:10, 03:13).
  - **Keeping the session active:** a pending background task keeps the container alive (observed 04:21–04:33), so long waits on the job's own progress let it run continuously.
  - **Consequence:** chemistry runs only while the session is active. An hourly watchdog yields about one 2-cycle chunk per hour, after SCF setup. Unattended long jobs cannot complete here, so the portable package is the route for completion.

- **Observed uptimes between reboots:** 19:52 → about 20:56 (about 64 min); 20:57 → about 21:29 (about 32 min); 21:30 → about 23:02 (about 92 min); 23:03 → about 00:05 (about 62 min); 00:06 → before 02:10 (last checkpoint 01:34).
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

4. **A checkpoint with an exactly matching operator was discarded when the restarted SCF converged to marginally different orbitals** (orbital energies beyond the 10⁻⁶ Eh validation tolerance, with the total energy identical to 10⁻¹¹ Eh). At 23:05 this discarded the 28-cycle L1_FULL checkpoint, and the next fresh save would have overwritten it.
   - **Now:** whenever the operator matches (exactly, or up to coordinate rounding) and the stored vectors have the right shape, they are used as an unvalidated initial guess, with no stage credit; the solve still converges to the full tolerance on the current operator.
   - A checkpoint that is truly unusable is copied to `invalidated_<time>_latest.npz` before any fresh start.
   - The validation tolerance itself was not relaxed. Tests cover both paths.
   - The 28-cycle checkpoint is also backed up in `phase1/results/ckpt_backup_2310_L1_28cyc/`.

## Convergence observation (01:08 UTC, not acted on)

- **The residual-0.01 stage of L1_FULL (B3LYP, full TD, 30 roots) has stagnated.** After cycle 10 the converged-root count oscillates between 21 and 27 of 30 over 44 cycles (24/30 at cycle 18, 26/30 at cycle 42, 27/30 at cycle 54). No stage has completed.
- **A plausible cause is the 2-cycle chunking itself:** each chunk restarts Davidson from the 30 Ritz vectors only, discarding the rest of the subspace (`docs/compute_reliability.md`). The higher roots of a dense manifold may then not converge.
- **Not changed:** chunk size, root count, tolerance and TD/TDA are held fixed, as required. A larger chunk lowers restart overhead but raises the work lost per reboot.
- **Decision needed:** whether to test a larger chunk (e.g. 6–10 cycles) in a matched comparison, or to accept the current rate. With 30 roots, two functionals and three stages remaining, the time to completion cannot be estimated while the first stage is not converging.

## Projected outcome of L1_FULL as configured (08:28 UTC)

- **The unconverged root is the highest requested one** (root 30, about 7.67 eV; roots 25–30 span 7.37–7.67 eV). It has not reached residual 0.01 in about 100 cycles. This is typical when the top requested root couples to states just outside the window, made worse by restarting the subspace every 2 cycles.
- **Budget:** `--td-max-cycles` is 400 per solve. With 216 cycles left and the 10⁻³ … 10⁻⁵ stages still ahead, the B3LYP solve will very likely exhaust the budget and end "not converged". That status is diagnostic only, never a success (`test_max_cycle_without_convergence_is_not_success`). The CAM-B3LYP solve would follow.
- **Options (none applied; each changes the computation and needs approval):**
  1. Request extra roots (e.g. 34) and judge convergence on the lowest 30.
  2. Use a larger chunk.
  3. Raise the cycle budget.

  All affect cost or definitions, not the visible-band physics. The visible-relevant roots are converged at the 0.01 stage.

## Not changed

Root counts, TD vs TDA, residual tolerances, chunk size (2 cycles), basis, functional, solvent.

**One-iteration chunks** were NOT adopted:
- each chunk restarts Davidson from Ritz vectors, which can increase total work (`docs/compute_reliability.md`);
- it would cut the Level 2 lost-work window from about 90 to about 45 min, still above the 32-min window observed;
- it would need a matched small-system convergence/timing test before use.

## Recommendation

Run L2_CAM_TDA15 (and the queued Level 2 / Fe jobs) on a stable machine or scheduler: `docs/portable_chemistry_package.md`. No paid compute has been provisioned.
