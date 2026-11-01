# Swiss original-data recovery and field comparison

This additive package resumes the ETH/Swiss field work. It contains executable analyses and original numerical records, not a literature-only review. It does **not** establish consistent field superiority of the complete glacier/cliff model.

## Original observations recovered

The accepted-stage Ogier et al. preprint discloses an author mirror at https://people.ee.ethz.ch/~werderm/pm-data/ . Selected archive members were recovered using bounded HTTP ranges and checked by original ZIP CRC32, uncompressed size and saved SHA256. `source/` contains the complete selected numerical working set, including the high-frequency processed temperature series. The entire 337,848,233-byte archive was not downloaded. The September mirror's identity with the later DOI deposit is not established.

Original study: Ogier et al. (2021), https://doi.org/10.5194/tc-15-5133-2021 . Data record and license attribution: https://doi.org/10.3929/ethz-b-000504956 . Original authors retain their rights; no new ownership or relicensing is asserted. See `original_mirror/ACCESS_AUDIT_RESUMED.md` and source manifests.

## What ran

- **Temperature processing:** verified calibration offsets; recovered numerical observations replace digitization as primary inputs. Native-time integration preserves gaps and changing sample cadence.
- **Controlled gap tests:** predeclared 1/6/24-hour temperature masking tests compare three reconstructors without incision targets. Linear bridging wins all six station-duration groups, but actual missingness need not behave like these selected observed blocks.
- **Thermal field comparison:** all frozen model families, scenarios, fitted coefficients, predictions, exclusions and numerical audits are retained. Direct survey elevation changes are the target; daily waterline marks are not relabelled as bed melt.
- **Local geometry development:** uses local stage-based hydraulic diameter and consistently recomputed Reynolds number. This is an explicitly retrospective development analysis with conditional hydraulic inputs, not an untouched independent validation.
- **Hydraulic correction:** two independent code inspections found an incorrect downstream-depth sign in the archived head expression. The source-level correction is established; transferable predictive superiority is not. Its later-velocity test is retained even though it failed.
- **Evidence gate:** conditional paired-error bounds prevent a pooled improvement or supplied uncertainty ranges being misrepresented as consistent superiority.

P3 and P5 belong to one channel/drainage event. Their intervals share survey epochs and forcings. These are not independent glaciers. All water-temperature-driven predictions here are forced retrospective hindcasts, not forecasts issued before future water conditions occurred. Temperature near the bed is not independently measured bulk-water temperature; source clocks and reconstructed depth/discharge remain qualified inputs.

## Read these first

1. `CLAUDE_HANDOFF.md` — integration instructions and current decision.
2. `conditional_benchmark/RESULTS.md` — strict native gate and author-daily reconstruction.
3. `conditional_benchmark/RESULTS_V4_NATIVE_GAPS.md` and `RESULTS_V5_DEVELOPMENT.md` when present — versioned native-time and local-geometry analyses.
4. `original_mirror/HEAD_DROP_AUDIT.md` and `velocity_transfer/REPORT.md` — correction versus prediction.
5. `OBSERVATION_OPERATOR.md` — why sensor pressure, full water depth and water contact are distinct.

Use the supplied frozen protocols and exact command arguments in each analysis folder. Save copies of existing outputs before rerunning: several scripts intentionally write results next to their code or into a specified output directory. Standard-library modules cover the integrations, correlations and arithmetic audits; plotting uses Matplotlib/NumPy. Do not execute the source notebooks to run these analyses.
