# Fixed physical footprint and float64 run bundles

`capture` fixes a rectangular XY footprint on cell-aligned coordinates. All selected cells must be ice and ghosts are excluded. `observe` computes geometric bed-volume loss independently of the model's heat accumulator. Optional heat-implied volume must use the identical mask and cumulative interval; mismatches refuse. This is a floor-melt observation, not vertical face melt or cliff lifetime.

`compare` requires equal footprint, coordinate frame, area and physical-time endpoint. It always reports signed absolute error and reports relative error only above a predeclared positive reference-volume threshold. It never converts a small absolute difference into a claim of negligible error without physical context.

`save_bundle` preflights configuration, source identifiers, arrays, geometry, evidence class and strict JSON before writing. It preserves float64 arrays and source/configuration hashes. The JSON manifest is the completion marker, written after the data. A crash can leave an orphan; such a partial bundle must be retained/refused or explicitly investigated, never silently treated as complete. No-overwrite publication uses exclusive hard links.

`load_bundle` refuses NaN, duplicate JSON keys, wrong masks/areas/bounds, unsupported evidence labels, changed configuration, corrupted bytes, mismatched geometric/heat volumes and incompatible schemas. Its hashes establish byte integrity only. They do not independently authenticate the runner, source truth or model applicability. No pickled arrays are loaded. Legacy float32 snapshots and reports cannot be upgraded into authenticated v4 data by this exporter.

The final short-run results are under `results/final_physical_smoke/`. Their source hashes match the delivered solver, vendored kernel, runner and observation module. Earlier pre-repair smoke files are kept separately outside this release; no result was overwritten. Thirteen observation tests pass. The three short numerical runs establish neither long-term convergence nor field accuracy.
