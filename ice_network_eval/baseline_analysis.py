"""Existing-method-surrogate comparison on the frozen S1-S5 problems (rules: cases/FROZEN_ACCEPTANCE.md section 6). Targets: thermal reference 400 cells; no-melting reference 200 cells for T0."""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import analyze_campaign as A
CFG = json.load(open(os.path.join(HERE, "competing_outlets_v1.json"))); I_OF = {s["id"]: s["inflow_m3_s"] for s in CFG["scenarios"]}
base = A.load(lambda r: A.case_set(r["spec"]) == "baselines"); fro = A.load(lambda r: A.case_set(r["spec"]) == "frozen" and r["status"] == "completed")
def find(rs, scn, **kw):
    for r in rs:
        s = r["spec"]
        if s["scn"] == scn and all(s.get(k) == v for k, v in kw.items()) and r["status"] == "completed": return r
rows = []; det = {}; failed = [(r["id"], r["status"]) for r in base if r["status"] != "completed"]
for scn in ("S1", "S2", "S3", "S4", "S5"):
    tgt = find(fro, scn, kind="reference", ncell=400); ts = A.series(tgt); tmx = tgt["metrics"]; I = I_OF[scn]
    for st in (17, 33):
        th = find(fro, scn, kind="reduced", stations=st); thm = th["metrics"]; the = A.ser_err(A.series(th), ts, I)["hydrograph_rmse_over_I"]
        entries = [("thermal (this model)", th)]
        for mn in ("T0nomelt", "T2dissall", "T3dissall_fixedaspect"):
            r = next((x for x in base if x["spec"]["scn"] == scn and x["spec"]["kind"] == "reduced" and x["spec"]["stations"] == st and x["spec"]["label"].endswith(mn) and x["status"] == "completed"), None)
            if r: entries.append((mn, r))
        for name, r in entries:
            m = r["metrics"]; e = A.ser_err(A.series(r), ts, I)["hydrograph_rmse_over_I"]
            relinc = [(m["sill_lowering%d_m" % k] - tmx["sill_lowering%d_m" % k]) / max(tmx["sill_lowering%d_m" % k], 1e-9) if tmx["sill_lowering%d_m" % k] > 1e-4 else None for k in (1, 2)]
            sw_ref, sw_b = tmx["t_dom2_s"] is not None, m["t_dom2_s"] is not None
            flags = []
            if A.route(m) != A.route(tmx): flags.append("outlet differs")
            if sw_ref != sw_b: flags.append("switching %s" % ("lost" if sw_ref else "spurious"))
            elif sw_ref and abs(m["t_dom2_s"] - tmx["t_dom2_s"]) > max(0.15 * tmx["t_dom2_s"], abs((thm["t_dom2_s"] or tmx["t_dom2_s"]) - tmx["t_dom2_s"])): flags.append("t_dom2 off %.0f s" % (m["t_dom2_s"] - tmx["t_dom2_s"]))
            if any(x is not None and abs(x) > 0.25 for x in relinc): flags.append("incision %s" % "/".join("n/a" if x is None else "%+.0f%%" % (100 * x) for x in relinc))
            rows.append([scn, f"r{st}", name, A.route(m), "-" if m["t_dom2_s"] is None else f"{m['t_dom2_s']:.0f}", f"{e:.4f}", "%.1f / %.1f" % (m["sill_lowering1_m"] * 1e3, m["sill_lowering2_m"] * 1e3), "/".join("n/a" if x is None else "%+.0f%%" % (100 * x) for x in relinc), f"{m['delivered1_m3']:.0f} / {m['delivered2_m3']:.0f}", "; ".join(flags) or "-" if name != "thermal (this model)" else "(reference for flags)"])
    r0 = find(base, scn, kind="reference")
    if r0: m = r0["metrics"]; rows.append([scn, "ref200", "T0 no-melt reference", A.route(m), "-" if m["t_dom2_s"] is None else f"{m['t_dom2_s']:.0f}", f"{A.ser_err(A.series(r0), ts, I)['hydrograph_rmse_over_I']:.4f}", "%.1f / %.1f" % (m["sill_lowering1_m"] * 1e3, m["sill_lowering2_m"] * 1e3), "-", f"{m['delivered1_m3']:.0f} / {m['delivered2_m3']:.0f}", "fixed-terrain resolved routing"])
    rows.append([scn, "ref400", "thermal reference (target)", A.route(tmx), "-" if tmx["t_dom2_s"] is None else f"{tmx['t_dom2_s']:.0f}", "0", "%.1f / %.1f" % (tmx["sill_lowering1_m"] * 1e3, tmx["sill_lowering2_m"] * 1e3), "-", f"{tmx['delivered1_m3']:.0f} / {tmx['delivered2_m3']:.0f}", "target"])
head = ["scn", "cfg", "incision law", "outlet (larger volume)", "t_dom2 (s)", "hydrograph RMSE/I vs thermal ref", "sill lowering R1/R2 (mm)", "rel. diff vs ref (R1/R2)", "delivered R1/R2 (m3)", "flags (rule section 6)"]
json.dump(dict(head=head, rows=rows, failed=failed), open(os.path.join(HERE, "records", "baseline_summary.json"), "w"), indent=1)
for r in rows: print(r)
print("failed:", failed)
