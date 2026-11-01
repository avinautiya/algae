"""Evaluation analysis (M3): results table incl. failures, contact-threshold sensitivity, comparisons with the committed reduced forecasts, numerical-error checks (coarse vs medium), and a machine-readable summary.
`python -m ice_cliff_survival.reactivation.analyze_eval EVAL_DIR SPEC_DIR OUT_DIR`"""
import glob, json, os, sys
import numpy as np
from . import metrics as M
def load(d, suffix):
    out = {}
    for f in sorted(glob.glob(os.path.join(d, "E*_%s.json" % suffix))):
        r = json.load(open(f)); n = os.path.basename(f)[:-5].rsplit("_", 1)[0]; out[n] = r
    return out
def err_stats(f_cfd_days, f_cfd, f_pred_days, f_pred):
    d = np.arange(0.0, min(f_cfd_days[-1], f_pred_days[-1]) + 1e-9, 0.25); a = np.interp(d, f_cfd_days, f_cfd); b = np.interp(d, f_pred_days, f_pred); return dict(mean_abs=float(np.mean(np.abs(a - b))), rms=float(np.sqrt(np.mean((a - b) ** 2))), max_abs=float(np.max(np.abs(a - b))))
def main(ev, specdir, out):
    os.makedirs(out, exist_ok=True); spec = json.load(open(os.path.join(specdir, "production_spec.json"))); U = spec["numerical_uncertainty"]; Uloss = U["contact_loss_day"]["U_used"]; Uf = U["flow_partition_vs_independent_column_weir"]["U_f_used"]
    med, coa = load(ev, "medium"), load(ev, "coarse"); table = []; comp = []
    for rr in spec["evaluation_runs"]:
        n = rr["name"]; r = med.get(n); c = coa.get(n); row = dict(name=n, layout=rr["layout"], removal_day=rr["removal_day"], status=(r or {}).get("status", "MISSING"))
        if not r or r.get("status") != "complete": row["error"] = (r or {}).get("error"); table.append(row); continue
        r["name"] = n; s = M.summarize(r); row.update({k: s[k] for k in ("loss", "rewet_episodes", "sustained_reconnection", "f0", "f_end", "post_removal_mean_f", "delivered_fraction", "toe_melt_reach_m3", "byp_melt_reach_m3", "crest_toe_end", "crest_byp_end", "max_water_resid", "max_energy_resid", "wall_s", "maxrss_MB", "cells")})
        row["theta_sensitivity"] = {str(t): {k: M.summarize(r, t)[k] for k in ("loss", "rewet_episodes", "sustained_reconnection")} for t in M.THETAS}
        if c and c.get("status") == "complete":
            c["name"] = n; sc = M.summarize(c); a, b = M.arrays(r), M.arrays(c); row["coarse_vs_medium"] = dict(loss_diff=None if (s["loss"] is None or sc["loss"] is None) else sc["loss"] - s["loss"], f_end_diff=sc["f_end"] - s["f_end"], mean_abs_f_diff=float(np.mean(np.abs(np.interp(a["day"], b["day"], b["f"]) - a["f"]))), toe_melt_rel_diff=(sc["toe_melt_reach_m3"] - s["toe_melt_reach_m3"]) / max(s["toe_melt_reach_m3"], 1e-12), reconnection_same=(sc["sustained_reconnection"] is None) == (s["sustained_reconnection"] is None), coarse_wall_s=sc["wall_s"])
        table.append(row); a = M.arrays(r)
        for m in spec["reduced_predictions"].get(n, []):
            es = err_stats(a["day"], a["f"], np.array(m["series_days"]), np.array(m["series_f"])); lo = m["loss"]; cl = s["loss"]
            comp.append(dict(run=n, method=m["method"], pred_loss=lo, cfd_loss=cl, loss_err=(None if (lo is None or cl is None) else lo - cl), loss_outcome_match=((lo is None) == (cl is None)), pred_reconnection=m["sustained_reconnection"], cfd_reconnection=s["sustained_reconnection"],
                             reconnection_outcome_match=((m["sustained_reconnection"] is None) == (s["sustained_reconnection"] is None)) if (rr["removal_day"] is not None) else None, f_err=es, pred_f_end=m["f_end"], cfd_f_end=s["f_end"], toe_melt_pred=m["toe_melt_m3"], toe_melt_cfd=s["toe_melt_reach_m3"], toe_melt_rel_err=(m["toe_melt_m3"] - s["toe_melt_reach_m3"]) / max(s["toe_melt_reach_m3"], 1e-12), within_numerical_uncertainty_loss=(None if (lo is None or cl is None) else abs(lo - cl) <= Uloss), within_numerical_uncertainty_f=es["mean_abs"] <= Uf))
    json.dump(dict(table=table, comparisons=comp, U_loss_day=Uloss, U_f=Uf), open(os.path.join(out, "eval_analysis.json"), "w"), indent=1, default=float)
    L = ["| run | status | loss day (θ=1e-3) | rewet episodes | sustained reconnection | f0 → f_end | post-removal mean f | delivered toe-water fraction | toe-reach melt m³ | bypass melt m³ | max water resid | max energy resid | wall s | RSS MB |", "|" + "---|" * 14]
    fm = lambda x, p=3: "—" if x is None else (f"{x:.{p}f}" if isinstance(x, float) else str(x))
    for r in table:
        if r["status"] != "complete": L.append(f"| {r['name']} | **{r['status']}** | {r.get('error')} |" + " |" * 11); continue
        L.append(f"| {r['name']} | complete | {fm(r['loss'])} | {r['rewet_episodes']} | {fm(r['sustained_reconnection'])} | {r['f0']:.3f} → {r['f_end']:.3f} | {fm(r['post_removal_mean_f'])} | {r['delivered_fraction']:.4f} | {r['toe_melt_reach_m3']:.4f} | {r['byp_melt_reach_m3']:.4f} | {r['max_water_resid']:.1e} | {r['max_energy_resid']:.1e} | {r['wall_s']:.0f} | {r['maxrss_MB']:.0f} |")
    L += ["", "### Comparison with committed reduced forecasts", f"Numerical uncertainty: contact-loss day ±{Uloss:.2f} d; toe-flow fraction ±{Uf:.2f}", "| run | method | loss pred / CFD (d) | loss outcome | reconnection outcome | mean abs f error | toe melt rel err |", "|---|---|---|---|---|---|---|"]
    for c in comp: L.append(f"| {c['run']} | {c['method'][:26]} | {fm(c['pred_loss'],2)} / {fm(c['cfd_loss'],2)} | {'match' if c['loss_outcome_match'] else 'MISMATCH'}{'' if c['loss_err'] is None else ' (err %+.2f d)' % c['loss_err']} | {fm(c['reconnection_outcome_match'])} | {c['f_err']['mean_abs']:.3f} | {c['toe_melt_rel_err']:+.2%} |")
    open(os.path.join(out, "results_table.md"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
if __name__ == "__main__": main(*sys.argv[1:4])
