#!/usr/bin/env bash
# Production Phase 1 queue (all TD-DFT jobs of the paper), restart-safe: a step whose summary.json exists
# is skipped, and an interrupted optimisation resumes from its last trajectory frame. Re-run this script
# after any interruption (container restart, OOM); it carries on where it stopped.
#   setsid nohup scripts/phase1_production.sh > phase1/results/queue.log 2>&1 &
# Memory: --max-memory 6000 on a 16 GB machine (PySCF overshoots the limit by ~35 %; 7000 left too little
# headroom for side jobs: the VM restarted twice).
cd "$(dirname "$0")/../phase1" || exit 1
R=results
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
run() {   # run OUTDIR OPT_FUNCTIONAL ARGS...: up to 3 attempts, resuming from the last frame
  local out=$1 f=$2; shift 2; mkdir -p "$out"
  for i in 1 2 3; do
    [ -s "$out/summary.json" ] && return 0
    local start=()
    [ -s "$out/opt_${f}_last.xyz" ] && start=(--start-xyz "$out/opt_${f}_last.xyz" --no-xtb)
    echo "=== attempt $i $(date -u) $* ${start[*]}" >> "$out/run.log"
    python3 run_phase1.py "$@" "${start[@]}" --max-memory 6000 --verbose 4 --outdir "$out" >> "$out/run.log" 2>&1
    echo "exit $? $(date -u)" >> "$out/run.log"
  done
  [ -s "$out/summary.json" ]
}
spin() {  # spin-state check of an Fe complex at its final geometry
  local out=$1 mode=$2
  [ -s "$out/spin_check.json" ] || [ ! -s "$out/opt_B3LYP_final.xyz" -a ! -s "$out/xtb_geometry.xyz" ] && return 0
  local g=$out/opt_B3LYP_final.xyz; [ -s "$g" ] || g=$out/xtb_geometry.xyz
  python3 -c "
import json, qc, fe_complex as F
out = F.spin_check(qc.xyz_file_to_atom_block('$g'), '$mode')
json.dump({str(k): dict(E_Eh=v[0], S2=v[1], converged=v[2]) for k, v in out.items()},
          open('$out/spin_check.json', 'w'), indent=1)" >> "$out/spin_check.log" 2>&1
}
# Level 2 B3LYP optimisation was stopped at step 16 (results/level2/stop_criterion.json): single point + TD-DFT
run $R/level2 none --level level2 --skip-opt --start-xyz $R/level2/opt_B3LYP_final.xyz --tddft-functionals B3LYP
echo "LEVEL2 DONE $(date -u)" >> $R/level2/run.log
run $R/level2_cam none --level level2 --skip-opt --start-xyz $R/level2/opt_B3LYP_final.xyz --tddft-functionals CAM-B3LYP
run $R/level1 B3LYP --level level1 --tddft-functionals B3LYP CAM-B3LYP
run $R/level3_catecholate B3LYP --level level3_catecholate --tddft-functionals B3LYP --nstates 50
spin $R/level3_catecholate catecholate
run $R/level3_tropolonate none --level level3_tropolonate --tddft-functionals B3LYP --nstates 50 --xtb-geometry
spin $R/level3_tropolonate tropolonate
echo "QUEUE DONE $(date -u)"
