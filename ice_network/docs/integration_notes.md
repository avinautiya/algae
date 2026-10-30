# Integration notes: the Codex verification contribution

Source: repository avinautiya/algae, branch `codex/ice-network-verification`, commit `a11895a4d081e5cd15a367be7132861afe59258d`, directory `ice_network_verification/` (19 files, additive). It was developed on the older base `0fc0096...`; the merge base with this branch is `e87a6d14...`.

## How it was integrated
`git cherry-pick -x a11895a4` onto `claude/ice-network-m1` (commit c8892e8, authorship preserved). Nothing was reset or merged wholesale; no production file was touched by the cherry-pick; the independent reference calculations stay in `ice_network_verification/` and are not imported by the production package. Its 26 tests run unchanged and pass.

## What the contribution verified, stated exactly
* The seven analytical limiting-case **specifications** and a strict CSV **checker**. Its own identity-fixture tests test the harness, not a solver.
* Independently recomputed mass and thermal-energy **ledger** audits (the global ledger cannot validate geometry; geometry is checked separately here, see below).
* A reference-error order utility.
* A **numerical audit of the thermal kernel** as of production commit `7b1d59ad`, executed with the Python function bodies (no Numba JIT), so it says nothing about JIT behaviour or performance. Ordinary scenarios agreed with its 60-digit ODE reference to 1.4e-9 W; a stress setting (`nu_factor` = 1e-12, thermal exposure near 1e-13) showed cancellation in the equilibrium-temperature algebra.

## What was done with it here (all in this branch)
1. **Stable algebra ported** (`thermal.channel_heat.segment_heat`): the Codex formulae with `expm1` and a series for 1 - F/a, using the production constants and temperature datum (T_m = 0). The entrance dissipation is still added to segment 0 only; the earlier equilibrium-temperature form is kept only as a test oracle. No clipping. Zero conductance and zero flow are handled outside the formula (zero flow never reaches it). Numba compilation verified (it is `@njit`).
2. **Regression test** `tests/test_thermal_stability.py` (180 cases): agreement with a 60-digit Decimal evaluation for exposures from 0 to 50, with and without dissipation and with zero or positive excess temperature, nonuniform geometry, different segment lengths and entrance dissipation; reproduces the legacy defect (negative heat to ice at exposure 1e-13) and shows the new form fixed. The stress case is a numerical stress, not a plausible closure; default runs are not claimed to have been wrong.
3. **Ordinary results unchanged**: re-running the budget tables and the refinement study with the new algebra changed no non-residual quantity by more than 5e-10 (relative).
4. **Compatible limiting cases run** through test-only adapters (`tests/adapters.py`, `examples/export_codex_cases.py`): linear_reservoir, finite_drainage, symmetric_outlets, unequal_outlets (prescribed linear/Torricelli edge laws on the production storage/ODE/ledger code: **these do not test Manning/GVF hydraulics**), pre_overtopping (real production Manning/GVF edge), heat_to_melt (production segment kernel with the manifest's c_p = 4180 passed as an argument, because the production constant is 4217 J/kg/K; no production constant was changed). Result with the Codex checker at rtol 1e-6, atol 1e-9: **6 of 7 pass**, largest absolute error 1.5e-10 (`records/codex_cases_report_tight.json`). **thermal_front is incompatible** with the reduced model (no transport PDE); it is exercised through the independent resolved reference solver instead (`ice_network_reference/tests`).
5. Claims kept apart: (1) the checker works; (2) a numerical kernel reproduces a limiting solution (6 cases); (3) the coupled solver converges (`docs/verification_report.md`); (4) physical formulation matches experiments: **not shown**; (5) field prediction: **not shown**.

## Not done / caveats
* The seven cases are idealised and use prescribed constitutive laws; passing them is unit-level evidence.
* The Codex audit pinned an older production commit; the kernel has since changed (stable algebra), so its numbers describe the old kernel; the regression test is the current evidence.
