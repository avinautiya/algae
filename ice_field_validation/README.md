# Public field hindcast benchmark

This additive package retains mixed field results; it does not establish field accuracy of the coupled stream/debris/cliff system. See FIELD_STATUS.md for scores and measurement limitations.

Inputs: Petersen et al. Kennicott data, DOI 10.1002/gdj3.70032, Zenodo 14625691. Original archive: 9,948,630 bytes, MD5 f740fc1d0e6ff19c1b8049b427db0c3a, SHA256 a243a485eb83c51f54bf135cfecca9ca73e59ecc052fbe0ef94fbb70cca8044f. Included critical text files are unchanged source subsets; SHA_MANIFEST.json records their hashes. Derived observed_intervals.csv is an observation transformation, not raw data. Source rights and attribution remain with original authors.

Reproduce in an isolated Python environment with requirements.txt installed:

```sh
export ALGAE_FIELD_DATA_ROOT="$PWD/kennicott"
python -m unittest discover -s event_benchmark -p 'test_*.py'
python audit_saved_metrics.py
python make_field_figures.py
```

Model reruns overwrite saved outputs: first make a separate copy. Exact executed 2022 source is event_benchmark/run_2022_frozen.py; exact executed 2023 source is secondary2023_corrected.py. run.py and secondary2023.py are the reviewed equivalent helper versions; the latter normalizes year-dependent merged weather-column names. Run the historical source snapshots in a copied directory to reproduce their code hashes. 2023 imports run.py; use the supplied version. Melt fitting is 2021 only; weather transfer is 2021–2022. All models are retained; 2023 was a secondary test after inspecting 2022. No retuning occurred here.

Missing stream temperature/discharge/contact at the stake cliffs is unknown, not observed zero. Site histories and adjacent event times are dependent. The shortwave sensitivity is an engineering scenario, not a confidence interval. Calibration albedo reaches its bound. No accuracy or award guarantee is implied.

Additional Ogier public-source matching is ongoing in the parent workspace; it is not included as a completed field score in this core release. No Claude worktree or existing simulation results are modified by this package.
