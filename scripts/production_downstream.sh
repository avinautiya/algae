#!/usr/bin/env bash
# Downstream production runs on the Phase 1 Level 2 B3LYP spectrum (tier D default), sequential and
# restart-safe (finished steps are skipped). Re-run after an interruption.
#   setsid nohup scripts/production_downstream.sh > results_downstream.log 2>&1 &
cd "$(dirname "$0")/.." || exit 1
W=${WORKERS:-2}
mkdir -p phase4/results phase3/results
if [ ! -s phase4/results/satellite_validation_tddft/satellite_validation_summary.json ]; then
  python3 phase4/satellite_validation.py --phenol tddft --workers $W --outdir phase4/results/satellite_validation_tddft
fi
python3 phase4/multi_scene.py --phenol tddft --workers $W --outdir phase4/results/multi_scene_tddft
if [ ! -s phase3/results/full/tables/mc_summary.csv ] || [ ! -s phase3/results/full/figures/FigS3_sobol_convergence.png ]; then
  python3 phase3/run_phase3.py --n-mc 1000 --n-sobol 1024 --workers $W --reuse --outdir phase3/results/full
fi
python3 phase3/literature_comparison.py --n 200 --workers $W --outdir phase3/results/literature_comparison
echo "DOWNSTREAM DONE $(date -u)"
