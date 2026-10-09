# Project status (2026-10-09, updated with each result)

## Completed evidence

| Item | Result | File |
|---|---|---|
| `heldout_v2` abundance retrieval (frozen) | formal rule met; weak evidence (log score only; D1–D3) | `records/heldout_v2/`, P5-CORR-1 |
| **H1** plot albedo, pigment/cell vs established optics | **NOT SUPPORTED** (−0.0019 MAE; algae vs none +0.031) | `records/h1_albedo/`, P6-H1-1 |
| **H2** station absorbed SW | **NOT SUPPORTED** (Tier A best: MAE 33.3 W m⁻²; M3 biased bright) | `records/h2h4/`, P6-H2H4-1 |
| **H4** vs empirical satellite albedo | **NOT SUPPORTED** (+0.004 albedo, interval spans 0) | same |
| **H5** station albedo from biology under the radiometer (PROMBIO) | **insufficient evidence** (2 test station-days; optics indistinguishable) | `records/h5/`, P6-H5-1 |
| **H3** ablation | **UNTESTED** (prerequisites met at ≤ 1 usable station-year; rule-dependent) | P5-H3-1, P6-SW-1 |
| Matched four-band forward, physics vs ML | every physics variant worse than the training mean; visible too bright | `records/ml_comparison_physics/`, P6-FWD-1 |
| Data-use register, H1 lineage audit, D3 correction | — | `docs/data_use_register.md` |
| PROMBIO roles frozen and matching manifest | 12 development, 5 test station-days; usable for radiation: 3 / 2 | `docs/new_data_inventory.md`, `records/prombio_matching/` |

## Running (all under the shared budget)

| Job | Class | Restart boundary |
|---|---|---|
| Chemistry L2_CAM_TDA15, L1_FULL (`phase1/jobs.py` runner) | foreground, 7000 + 3000 MB | TD checkpoints every 2 iterations |
| Forward-optics diagnostics (`phase4/forward_diagnostics.py`) | background, 1500 MB | per optics model; calibrations cached on disk |

## Queued

| Job | Waiting for |
|---|---|
| Surrogate re-benchmark under the new qualification binding | memory (1500 MB) |
| Scene product on S2A_22WEV_20190723 (3 km) | after diagnostics |

## Missing inputs (require external action or new data)

| Input | Status |
|---|---|
| Naegeli et al. (2017) published coefficients | blocked (403 / bot protection); M2 remains the stated fallback |
| UPE_U 2018 random-grid counts and UAS | not acquired (BAS PDC) |
| Hyperspectral scenes with co-located labels | none identified |
| Joint calibration draws through emulators | expensive (emulator rebuild per draw); waits for chemistry headroom |

## Can finish without expensive compute

- The H1 spectral secondaries come from the running diagnostics.
- k-free `heldout_v2` sensitivity (reuses the cached physics tables).
- Host-coupling example (reference SEB, Mode A/B).
- Effect-ranking example on the scene product (needs the posterior-mean draw; flagged single-draw).
- Claims/limitations audit and reproduction commands.
