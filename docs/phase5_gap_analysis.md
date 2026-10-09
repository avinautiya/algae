# Gap analysis: current code and results vs H1–H4 (2026-10-09)

| Hypothesis | What exists | Gap | Status |
|---|---|---|---|
| H1 independent albedo | Forward emulator (BBA, bands) for all optics; `heldout_v2` forward test on nadir HCRF only. Physics models over-predict broadband HCRF by +0.1 to +0.3, and training-site band climatology beats every physics model by ~6 nats. | (a) No test against hemispherical albedo, although 87 S6 albedo spectra exist (47 with counts). (b) The source of the HCRF bias (geometry/k, ice state, optics, illumination) is not separated. (c) Nuisance state comes only from priors fitted on the same surfaces' HCRF. | not tested |
| H2 absorbed SW | SEB uses measured SW↓; the multi-scene product gives clear-sky forcing. | No prediction of station absorbed SW from satellite-derived albedo; KAN_L hourly data not ingested; no station-pixel extraction. | not tested |
| H3 ablation | SEB repaired and verified numerically; validated against ablation with MEASURED albedo (over-predicts 1.5–1.7× in 2016–2018); 2019 sensors disagree 2.5×. | The prerequisite (an SEB that reproduces ablation with measured albedo) is not met at KAN_M; χ fitted, not independent. | expected UNTESTED at KAN_M |
| H4 beyond empirical satellite albedo | Retrieval and forward model exist. | No empirical narrowband-to-broadband baseline; no station albedo target pipeline. | not tested |
| Uncertainty | Joint draws exist in Phase 3 only. | Held-out and emulator use posterior-mean optics (deviation D1); no draw-consistent chain. | gap |
| Satellite product | `run_phase4` maps (abundance, BBA, RF) plus summaries. | No interpretation classes (snow, water, dust vs algae ambiguity), no versioned schema, no NetCDF/COG with validated round trip, no per-scene provenance sidecar. | gap |
| Host coupling | Offline SEB only. Branch `codex/glacier-model-integration` exists (not yet reviewed). | Modes A/B not defined in code; double-counting rules not enforced. | gap |
| Surrogate | Code with qualification gate; benchmark running (tolerances set before outcomes). | High SZA, boundary states and stale-qualification binding to source/config to be verified; rejection outside domain. | in progress |
| Compute | Per-job scratch, 2-iteration TD checkpoints. | Memory budget uses host MemTotal (16 GB), not the cgroup limit (14.35 GB, `claude-code-bash`); analysis jobs not budgeted; thread pools only partly limited; scoring not resumable. | gap (first) |

## Resource-safe implementation sequence

1. **Shared resource budget** (`common/resources.py`): cgroup-aware memory limit and usage, thread caps for OMP, OpenBLAS and MKL, disk admission. Used by the Phase 1 runner and by every analysis launcher (`common/run_budgeted.py`). Tested with fixtures.
2. **Forward-optics diagnostics on H1 data** (cheap: forward emulator evaluations only, under 1 GB).
   - Optics-only with known abundance, one factor varied at a time.
   - Run only after step 1, and only while the chemistry headroom allows.
3. **H1 primary analysis** per the protocol (leave-one-day-out). Small.
4. **Ingest KAN_L hourly; extract station pixels** from Sentinel-2 scenes (network, small) for H2/H4.
5. **Scene-interpretation product and output schema**, validated on ONE scene first (moderate: one emulator set, about 1.5 GB peak with chunking).
6. **Host-coupling modes:** review the codex branch; offline reference-SEB coupling example.
7. **Joint-draw uncertainty follow-up** (expensive: emulator rebuild per draw). Admitted only under the budget, and after the Level 2 chemistry finishes or is checkpointed.
