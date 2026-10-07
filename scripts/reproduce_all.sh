#!/usr/bin/env bash
# Regenerate every table and figure (600-dpi PNG + vector PDF) from the production Phase 1 output.
#
#   BIOSNICAR_PATH=/path/to/biosnicar-py bash scripts/reproduce_all.sh [phase1_results_dir]
#
# Phase 1 itself (DFT/TD-DFT, hours on a CPU) is not rerun here; produce it with
#   python phase1/run_phase1.py --level level2 --max-memory 4000   (and --level level1)
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
P1="${1:-$ROOT/phase1/results}"
W="${WORKERS:-$(nproc)}"
test -f "$P1/level2/summary.json" || { echo "missing $P1/level2 (run Phase 1 first)"; exit 1; }
L1=()
[ -f "$P1/level1/summary.json" ] && L1=(--phase1-l1 "$P1/level1")

echo "== tests";             python -m pytest -q "$ROOT"/phase{1,2,3,4}/tests
echo "== Phase 2 (Fig. 2A-C, S0-S2)"
python "$ROOT/phase2/run_phase2.py" --phase1-l2 "$P1/level2" "${L1[@]}" --outdir "$ROOT/phase2/results"
echo "== Phase 3 (Fig. 3A-C, S3)"
python "$ROOT/phase3/run_phase3.py" --phase1-l2 "$P1/level2" --workers "$W" --outdir "$ROOT/phase3/results"
echo "== Phase 4 synthetic + real scene (Fig. 4A-C, S4, S5)"
python "$ROOT/phase4/run_phase4.py" --source synthetic --phase1-l2 "$P1/level2" --workers "$W" \
       --outdir "$ROOT/phase4/results/synthetic"
python "$ROOT/phase4/run_phase4.py" --source s2 --phase1-l2 "$P1/level2" --workers "$W" --resolution 40 \
       --outdir "$ROOT/phase4/results/real"
echo "== Phase 4 field-site bias study"
python "$ROOT/phase4/bias_study.py" --phase1-l2 "$P1/level2" --phenol tddft --workers "$W" \
       --outdir "$ROOT/phase4/results/bias_study_tddft"
python "$ROOT/phase4/bias_study.py" --phase1-l2 "$P1/level2" --phenol williamson2020 --workers "$W" \
       --outdir "$ROOT/phase4/results/bias_study_williamson2020"
