# Research handoff: new datasets and effect interpretation

This additive branch extends `codex/ml-comparison` with:

- `DATASETS.md`: prioritized verified sources and explicit acquisition/overlap limitations.
- `MODEL_DESIGN.md`: proposed physics-constrained hierarchical retrieval, observation operators, validation and novelty hypotheses. **Not a trained operational model.**
- `data/`: original CC0 PROMBIO 2024 CSV and 2021–2023 workbook, plus provider metadata. Both original data files match provider MD5 checksums.
- `prombio.py`: read-only 2024 ingestion preserving source values, missingness, censoring and unknown positions.
- `effects.py`: paired joint-ensemble radiative-effect summaries/rankings with mandatory evidence gates and unresolved comparisons.

Run:

```sh
python3 -m unittest discover -s agency_design/tests -v
python3 -m agency_design.prombio
```

Twelve tests passed locally. They verify software contracts and data ingestion, not scientific prediction. Effect-ranking tests use analytic fixtures; no fabricated scientific map or trained model is supplied. Real effects require validated paired radiative-transfer ensembles.

Claude should inspect these files and reconcile concurrent changes before merging. Do not overwrite frozen physics evidence or train on all new campaigns: first define acquisition/support matching and reserve independent campaign/site-year validation. Keep qualitative mineral ratios qualitative, station coordinates distinct from sample coordinates, and public imagery distinct from biological ground truth.
