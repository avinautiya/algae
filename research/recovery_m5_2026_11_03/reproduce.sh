#!/bin/bash
# Reproduce every number in EXPERIMENT_STATUS.md (total CPU time < 2 min). Run from this directory inside the isolated checkout.
set -e
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
REPO=$(git rev-parse --show-toplevel); S=${SCRATCH:-/tmp/recovery_m5_scratch}; mkdir -p "$S"
# 0. packet manifest + helper tests (copy helpers first: screen.py writes RESULTS.json)
python3 - <<PY
import json,hashlib,os
m=json.load(open("$REPO/research/glacier_campaign_20261103/MANIFEST.json"))
for p,h in m["sha256"].items():
    f="$REPO/research/glacier_campaign_20261103/"+p
    assert hashlib.sha256(open(f,"rb").read()).hexdigest()==h,p
print("manifest ok")
PY
rm -rf "$S/camp"; mkdir -p "$S/camp"; cp -r "$REPO/research/glacier_campaign_20261103/"{baseline,execution,scoring} "$S/camp/"
(cd "$S/camp/baseline" && python3 screen.py && python3 -m unittest -q test_screen)
(cd "$S/camp/scoring" && python3 -m unittest -q test_score_trace)
(cd "$S/camp/execution" && python3 -m unittest -q test_preflight)
# 1. protected preflight (66 pinned tests, installed-copy hashes, covered-source refusal)
rm -rf "$S/pf"; python3 "$REPO/research/glacier_campaign_20261103/execution/preflight.py" --repo "$REPO" --out "$S/pf" --mode physical_bare | tail -3
# 2. independent viability check, refinement, analytic lane, figures
python3 viability/independent_check.py "$S/viability_check.json"
python3 verification/refinement.py "$S/refinement.json"
python3 analytic_lane/run_analytic_lane.py "$S/analytic_lane_results.json"
python3 -m pytest -q "$REPO/ice_integration" "$REPO/research/readiness_2026_11_01"
