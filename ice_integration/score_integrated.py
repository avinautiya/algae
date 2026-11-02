"""Integrated evaluation scoring (spec: INTEGRATED_EVALUATION_SPEC.md, committed first). `python -m ice_integration.score_integrated OUT_DIR`
1) equivalence: integrated point model reproduces every saved V5 prediction (3 imputers x 27 scenarios x all accepted intervals), 2) independent recomputation of MAE/RMSE/bias per station and pooled, 3) consistency verdict per the fixed acceptance rule, 4) per-baseline win/loss counts. Retrospective evaluation on inspected development data."""
import json, math, os, subprocess, sys
from . import point_model as PM
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
V5 = os.path.join(ROOT, "ice_swiss_field_validation", "conditional_benchmark", "v5_local_geometry_development")
CANDIDATE = "db04_scale_fitted"; BASELINES = ["constant_rate_fitted", "constant_h_fitted", "db04_unfitted", "db03_unfitted", "gnielinski_smooth_unfitted"]
def main(out):
    os.makedirs(out, exist_ok=True); eq_max = 0.0; eq_n = 0; metric_gap = 0.0; table = {}; nominal = {}
    for imp in ("linear", "left", "right"):
        d = json.load(open(os.path.join(V5, imp, "CONDITIONAL_RESULTS.json")))
        for name, sc in d["scenarios"].items():
            if sc.get("status") != "completed": table[(imp, name)] = dict(status=sc.get("status")); continue
            cal = sc["calibration"]; test = [i for i in sc["accepted_intervals"] if i["role"] == "test"]; errs = {}
            for it in sc["accepted_intervals"]:
                p = PM.predict(it, cal)
                for fam in PM.FAMILIES:
                    ref = it["predictions_m"][fam]; rel = abs(p[fam] - ref) / max(abs(ref), 1e-300); eq_max = max(eq_max, rel); eq_n += 1
            for fam in PM.FAMILIES:
                for st in ("P3", "P5", "pooled"):
                    sel = [i for i in test if st == "pooled" or i["station"] == st]; errs[(fam, st)] = PM.metrics([PM.predict(i, cal)[fam] - i["observed_lowering_m"] for i in sel]) if sel else None
            for st in ("P3", "P5"):
                for fam in PM.FAMILIES:
                    saved = sc["test_metrics"].get(st, {}).get(fam)
                    if saved and errs[(fam, st)]: metric_gap = max(metric_gap, abs(saved["MAE_m"] - errs[(fam, st)]["MAE_m"]), abs(saved["RMSE_m"] - errs[(fam, st)]["RMSE_m"]), abs(saved["bias_m"] - errs[(fam, st)]["bias_m"]))
            table[(imp, name)] = dict(status="completed", errs={f"{k[0]}|{k[1]}": v for k, v in errs.items()}, n_test=len(test), excluded=len(sc.get("excluded_intervals", [])))
            if imp == "linear" and name == "clock+0h_T+0.00C_W2.00m": nominal = table[(imp, name)]
    # fixed acceptance rule
    cmp = {}; win = {b: {st: 0 for st in ("P3", "P5", "pooled")} for b in BASELINES}; total = 0
    for key, row in table.items():
        if row["status"] != "completed": continue
        total += 1
        for b in BASELINES:
            for st in ("P3", "P5", "pooled"):
                c = row["errs"].get(f"{CANDIDATE}|{st}"); o = row["errs"].get(f"{b}|{st}")
                if c and o and c["MAE_m"] < o["MAE_m"]: win[b][st] += 1
    nominal_ok = all(nominal["errs"][f"{CANDIDATE}|{st}"]["MAE_m"] < nominal["errs"][f"{b}|{st}"]["MAE_m"] for b in BASELINES for st in ("P3", "P5", "pooled"))
    consistent = all(win[b][st] == total for b in BASELINES for st in ("P3", "P5", "pooled")) and nominal_ok
    res = dict(spec="INTEGRATED_EVALUATION_SPEC.md", status="retrospective evaluation on inspected development data", equivalence_max_relative_difference=eq_max, equivalence_n=eq_n, equivalence_pass=eq_max < 1e-9, saved_metric_vs_recomputed_max_abs_difference_m=metric_gap,
               scenarios_scored=total, nominal=nominal, candidate=CANDIDATE, candidate_wins_over_baseline_by_station=win, nominal_candidate_beats_every_baseline_everywhere=nominal_ok, consistent_improvement=consistent,
               table={f"{k[0]}|{k[1]}": v for k, v in table.items()}, git=subprocess.run(["git", "-C", ROOT, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip())
    json.dump(res, open(os.path.join(out, "EVALUATION_RECORD.json"), "w"), indent=1, default=float)
    print("equivalence max rel diff %.2e over %d predictions; saved-vs-recomputed metric gap %.2e m" % (eq_max, eq_n, metric_gap)); print("scenarios scored:", total)
    print("nominal MAE (m) [P3 / P5 / pooled]:")
    for f in PM.FAMILIES: print("  %-28s %s" % (f, " / ".join("%.4f" % nominal["errs"][f"{f}|{st}"]["MAE_m"] for st in ("P3", "P5", "pooled"))))
    print("candidate wins vs baseline (of %d scenarios) [P3 / P5 / pooled]:" % total)
    for b in BASELINES: print("  %-28s %s" % (b, " / ".join(str(win[b][st]) for st in ("P3", "P5", "pooled"))))
    print("CONSISTENT IMPROVEMENT:", consistent)
if __name__ == "__main__": main(sys.argv[1])
