"""Matched-accuracy cost comparison from the serial timing set (rule in cases/FROZEN_ACCEPTANCE.md section 5). Accuracy from the frozen-set records (deterministic runs);
cost = median of the serial repeats of solver_s. Failed/slow repeats are kept and flagged (max/min repeat ratio > 1.25)."""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import analyze_campaign as A
S = json.load(open(os.path.join(HERE, "records", "campaign_summary_frozen.json")))
RU = json.load(open(os.path.join(HERE, "records", "reference_uncertainty.json")))
runs = A.load(lambda r: A.case_set(r["spec"]) == "timing")
tm = {}; bad = []
for r in runs:
    s = r["spec"]; key = (s["scn"], ("ref%d" % s["ncell"]) if s["kind"] == "reference" else "r%d" % s["stations"])
    if r["status"] != "completed": bad.append((r["id"], r["status"])); continue
    tm.setdefault(key, []).append(r["timing"])
rows = []; out = {}
for scn in ("S1", "S2", "S5"):
    v = S[f"{scn}|{{}}"]; e_ref = {int(n): x["hydrograph_rmse_over_I"] for n, x in v["reference_cost_vs_error"].items()}      # error of 100/200-cell refs vs the 400-cell target
    cost = {k: np.median([t["solver_s"] for t in tm[(scn, k)]]) for (s, k) in tm if s == scn}
    spread = {k: (max(t["solver_s"] for t in tm[(scn, k)]) / min(t["solver_s"] for t in tm[(scn, k)]), len(tm[(scn, k)])) for (s, k) in tm if s == scn}
    if "ref100" not in cost or "ref200" not in cost: continue
    lx = np.log([e_ref[200], e_ref[100]]); ly = np.log([cost["ref200"], cost["ref100"]])
    e400 = RU.get(f"convergence|{scn}", {}).get("hydrograph_rmse_over_I", {}).get("err_coarser_of_last_pair")
    for k, key in (("r9", "r9_dt900"), ("r17", "r17_dt450"), ("r33", "r33_dt225")):
        if k not in cost: continue
        err = v["reduced"][key]["hydrograph"]["hydrograph_rmse_over_I"]
        matched = e_ref[200] <= err <= e_ref[100]
        t_ref = float(np.exp(np.interp(np.log(err), lx, ly))) if matched else None
        rec = dict(err=err, solver_median_s=float(cost[k]), repeats=spread[k][1], max_over_min=float(spread[k][0]), matched=matched, ref_time_at_equal_error_s=t_ref, speedup=None if t_ref is None else t_ref / cost[k],
                   warmup_s=float(np.median([t.get("warmup_s", np.nan) for t in tm[(scn, k)]])), sampling_run_s=float(np.median([t.get("sampling_run_s", np.nan) for t in tm[(scn, k)]])), e400_hydro=e400)
        out[f"{scn}|{k}"] = rec
        rows.append([scn, k, f"{err:.4f}", f"{cost[k]:.1f} (n={spread[k][1]}, max/min {spread[k][0]:.2f})", "yes" if matched else "no (outside %.4f-%.4f)" % (e_ref[200], e_ref[100]), "-" if t_ref is None else f"{t_ref:.1f}", "-" if t_ref is None else f"{t_ref / cost[k]:.2f}", f"{rec['warmup_s']:.1f} / {rec['sampling_run_s']:.1f}"])
    for k in ("ref100", "ref200"):
        if k in cost: rows.append([scn, k, f"{e_ref[int(k[3:])]:.4f}", f"{cost[k]:.1f} (n={spread[k][1]}, max/min {spread[k][0]:.2f})", "-", "-", "-", "-"])
head = ["scn", "config", "hydrograph RMSE/I vs 400-cell ref", "serial solver s (median)", "matched to reference 100-200 cell interval", "reference s at equal error", "speed-up", "warm-up / sampling-run s"]
json.dump(dict(head=head, rows=rows, detail=out, failed_or_incomplete=bad), open(os.path.join(HERE, "records", "timing_serial_summary.json"), "w"), indent=1, default=float)
for r in rows: print(r)
print("failed/incomplete:", bad)
