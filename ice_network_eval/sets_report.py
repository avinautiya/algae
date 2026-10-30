"""Generates records/SETS_REPORT.md: interpolation (param_eval) and extrapolation (challenge) results, separately, from set_summary_*.json (frozen rules in cases/FROZEN_ACCEPTANCE.md)."""
import json, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); R = lambda n: os.path.join(HERE, "records", n)
L = []; P = L.append
def q(a, k=3): return "-" if not len(a) else "%s / %s / %s" % tuple(f"{float(np.percentile(a, p)):.{k}g}" for p in (50, 90, 100))
def viol(r):
    fl = r.get("flags") or {}; sh = r.get("shared_assumption_flags") or {}
    v = [n for n, ok in fl.items() if not ok] + (["bank_overtopped"] if sh and not sh.get("bank_not_overtopped") else []) + (["Re<2000 (Manning)"] if sh and not sh.get("manning_regime_plausible") else [])
    return ", ".join(v) or "none"
for name, title in (("param_eval", "Interpolation: parameter-evaluation set (64 Sobol points inside the development ranges)"), ("challenge", "Extrapolation: challenge set (32 Sobol points outside the development ranges on every axis)")):
    S = json.load(open(R(f"set_summary_{name}.json"))); pts = S["points"]; A = S["aggregate"]
    P(f"## {title}\n")
    for st in ("17", "9"):
        rows = [v["reds"][st] for v in pts.values() if "reds" in v and st in v["reds"]]; ok = [r for r in rows if r["status"] == "completed"]
        n = len(rows); cls = {}
        for r in rows: cls[r["class_"]] = cls.get(r["class_"], 0) + 1
        P(f"### r{st} ({'primary' if st == '17' else 'secondary'} configuration)\n")
        P(f"* **Execution success**: {len(ok)} of {n} runs completed ({n - len(ok)} failed/stalled/timed out). **Acceptable predictive performance** (class ok): {cls.get('ok', 0)} of {n} ({cls.get('ok', 0) / n:.0%}); {cls.get('failure_unexplained', 0)} unexplained, {cls.get('failure_outside_regime', 0)} outside regime, {cls.get('failure_numerical', 0)} numerical (execution failure).")
        P(f"* Within both models' assumptions (reduced envelope + shared checks): {sum(r['within_both_models'] for r in ok)} completed cases; classes there {A['r' + st]['classes_within_both']}; outside: {A['r' + st]['classes_outside_both']}.")
        deg = [r for k, v in pts.items() if "reds" in v and st in v["reds"] and v["route_ref"] == "none" and (r := v["reds"][st])["status"] == "completed"]
        cl = [r for r in ok if not r["near_tie"]]; nt = [r for k, v in pts.items() if "reds" in v and st in v["reds"] and v["route_ref"] != "none" and (r := v["reds"][st])["status"] == "completed" and r["near_tie"]]
        P(f"* **Degenerate cases** (the reference delivers no water through either route in 8 h; agreement is trivial and carries no information about routing): {len(deg)} of {len(ok)} completed cases; the reduced model also delivers nothing in {sum(r['route'] == 'none' for r in deg)}.")
        P(f"* Route (outlet) choice, clear cases (reference margin > 10 %): {sum(r['route_match'] for r in cl)} of {len(cl)} correct; genuine near-ties (reference margin <= 10 %, not scored): {len(nt)}; of those the reduced model picks the reference route in {sum(r['route_match'] for r in nt)}.")
        sw = [r for r in ok if r["switch_ref"] or r["switch_red"]]
        P(f"* Switching (R2 overtakes R1): reference {sum(r['switch_ref'] for r in ok)} cases, reduced {sum(r['switch_red'] for r in ok)}; agreement on existence {sum(r['switch_ref'] == r['switch_red'] for r in sw)} of {len(sw)} cases where either has it; time error where both: " + ("|err| median / max " + "%.0f / %.0f s" % (np.median([abs(r['t_dom2_err_s']) for r in ok if r['t_dom2_err_s'] is not None]), max(abs(r['t_dom2_err_s']) for r in ok if r['t_dom2_err_s'] is not None)) + f" (n = {sum(r['t_dom2_err_s'] is not None for r in ok)})" if any(r['t_dom2_err_s'] is not None for r in ok) else "none"))
        P("* Error distributions over completed cases (median / 90th percentile / maximum):")
        P("  * hydrograph RMSE / inflow: " + q([r["hydro_rmse_over_I"] for r in ok]) + "; cumulative delivered-volume error (relative to total): " + q([r["cumvol_err_rel_total"] for r in ok]))
        P("  * |t_act2 error| (s): " + q([abs(r["t_act2_err_s"]) for r in ok if r["t_act2_err_s"] is not None], 4) + f" (n = {sum(r['t_act2_err_s'] is not None for r in ok)})")
        P("  * |sill-lowering error| R1 / R2 (mm): " + q([abs(r["sill_err_m"][0]) * 1e3 for r in ok], 3) + "  /  " + q([abs(r["sill_err_m"][1]) * 1e3 for r in ok], 3))
        P("  * reduced solver time (s; two concurrent workers): " + q([r["solver_s"] for r in ok], 3))
        P("")
    fails = [(k, st, r) for k, v in pts.items() if "reds" in v for st, r in v["reds"].items() if r["class_"] != "ok"]
    P("Non-ok reduced runs (all failures kept):\n"); P("| point | config | class | hydrograph RMSE/I | route ok | near-tie | violated assumptions | note |"); P("|---|---|---|---|---|---|---|---|")
    for k, st, r in fails:
        note = (r.get("error") or "")[:90] if r["status"] != "completed" else ""
        if r["status"] != "completed" and not note: lp = r.get("last_progress") or {}; note = f"{r['status']} at simulated t = {lp.get('sim_t_s')} s"
        P(f"| {k} | r{st} | {r['class_']} | {('%.3f' % r['hydro_rmse_over_I']) if r['status'] == 'completed' else '-'} | {r.get('route_match', '-')} | {r.get('near_tie', '-')} | {viol(r) if r['status'] == 'completed' else 'not diagnosed (no completed run)'} | {note} |")
    P("")
    if name == "challenge":
        P("### Every challenge case (r17): is it within the physical assumptions of both models?\n")
        P("Reduced-model envelope: quasi-steady (p90 r_qs <= 1), subcritical (<= 5 % active snapshots with Froude >= 1), channel storage neglected (p90 <= 5 % of node storage). Shared: banks not overtopped (< 1 m), Manning plausible (median max Re >= 2000). The reference additionally represents supercritical flow and jumps.\n")
        P("| point | inflow m3/s | pulse amp / period h | deficit m | R2-sill offset m | slopes R1/R2 | ref outlet (margin) | r17 class | hydrograph RMSE/I | switching ref/red | violated assumptions | within both models |"); P("|---|---|---|---|---|---|---|---|---|---|---|---|")
        for k, v in sorted(pts.items()):
            p = v["params"]; r = v["reds"]["17"]
            P(f"| {k} | {p['inflow_m3_s']:.3f} | {p['pulse_amp']:.2f} / {p['pulse_period_h']:.1f} | {p['deficit_m']:+.2f} | {p['sill_offset_m']:+.3f} | {p['slope1']:.4f} / {p['slope2']:.3f} | {v['route_ref']} ({v['route_margin']:.2f}) | {r['class_']} | {('%.3f' % r['hydro_rmse_over_I']) if r['status'] == 'completed' else '-'} | " + (f"{r['switch_ref']}/{r['switch_red']}" if r['status'] == 'completed' else '-') + f" | {viol(r) if r['status'] == 'completed' else 'n/a (run did not complete)'} | {r.get('within_both_models', 'n/a')} |")
        P("")
open(R("SETS_REPORT.md"), "w").write("\n".join(L)); print("ok", len(L))
