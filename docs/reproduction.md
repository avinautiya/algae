# Reproduction commands (all compute through the shared budget)

**Environment:** `requirements.txt` (pinned versions), plus a BioSNICAR checkout at revision fe74eee (`BIOSNICAR_PATH=...`). Python 3.13. One numerical thread per job.

**Budget wrapper.** Every command below is prefixed by `python common/run_budgeted.py --name <n> --mem-mb <m> [--background|--short] --`. It admits the job against the container's memory cgroup, caps thread pools, applies an RLIMIT_AS backstop and gives each job its own scratch directory (`docs/compute_reliability.md`).

| Evidence | Command | Inputs (checksums) | Output |
|---|---|---|---|
| Held-out v2 (frozen) | `python3 phase4/heldout.py --outdir <new>` | `records/heldout_v2/FROZEN_PROVENANCE.json` | `records/heldout_v2/` (read-only) |
| k-prior sensitivity | `python3 phase4/heldout_kfree.py --variant independent` (and others) | as above | `phase4/results/heldout_v2_kfree/<variant>/` |
| H1 | `python3 phase4/h1_albedo.py --outdir <dir>` | BioSNICAR `Albedo_master.csv` (b4f645cb…), counts | `records/h1_albedo/` |
| H2/H4 population | `python3 phase4/station_pixels.py --outdir <dir>` then `python3 phase4/audit_h2h4.py` | `data/promice_raw/SHA256SUMS`; S2 L2A via earth-search | `records/station_pixels/`, `records/h2h4/audit.json` |
| H2/H4 scores | `python3 phase4/h2h4_station.py --pixels records/station_pixels/station_pixels.csv --outdir <dir>` | as above | `records/h2h4/` |
| H3 prerequisites | `python3 phase4/h3_prerequisites.py --outdir <dir>` | `data/promice_raw` | ledger P5-H3-1 |
| Low-sun SW audit | `python3 phase4/audit_low_sun.py` | `data/promice_raw` | `records/low_sun_audit/` |
| PROMBIO matching | `python3 -m agency_design.prombio_manifest --outdir records/prombio_matching` | provider MD5 verified in code | `records/prombio_matching/` |
| H5 | `python3 phase4/h5_station_biology.py --split development` then `--split test` | manifest above | `records/h5/` |
| Physics vs ML forward | `python3 phase4/physics_forward_bands.py` | `records/ml_comparison/` | `records/ml_comparison_physics/` |
| ML study (imported) | `python3 -m comparison_study.study --outdir <new>` | field CSVs | `records/ml_comparison/` |
| Forward diagnostics | `python3 phase4/forward_diagnostics.py --outdir <dir>` | `Albedo_master`, HCRF, counts | `phase4/results/forward_diagnostics/` |
| Host coupling | `python3 -m model_integration.example_reference_seb` | station pixels, `data/promice_raw` | `records/host_coupling/` |
| Scene product | `python3 phase4/scene_product.py --scene-id S2A_22WEV_20190723_0_L2A --size-km 3 --outdir <new>` | STAC item (hash in sidecar) | COG/NetCDF/sidecar/report |

**Tests (short class):**
```
python -m pytest common/tests phase1/tests phase4/tests
python -m unittest discover -s model_integration/tests
python -m unittest discover -s agency_design/tests
python -m pytest comparison_study/tests
```
