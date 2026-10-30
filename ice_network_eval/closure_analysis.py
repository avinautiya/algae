"""Closure sensitivity from per-run records: (a) numerical error of the reduced model under MATCHED closures (reduced and reference use the same variant), (b) the change of the
reference prediction between closure variants (physical-assumption sensitivity), (c) whether route choice / event times change.  Reference = 200 cells (not extrapolated; its own
mesh uncertainty is shown by the frozen-set 200->400 change)."""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import analyze_campaign as A
VAR = {"base": {}, "nu0.5": dict(nu_factor=0.5), "nu2": dict(nu_factor=2.0), "n0.03": dict(n_man=0.03), "n0.05": dict(n_man=0.05)}


def run(out=os.path.join(HERE, "records", "closure_summary.json")):
    runs = A.load(lambda r: r["status"] == "completed" and (A.case_set(r["spec"]) == "closure" or (A.case_set(r["spec"]) == "frozen" and (r["spec"].get("ncell") == 200 or r["spec"].get("stations") in (9, 17) and r["spec"]["dt_geom"] in (900, 450)))))
    idx = {}
    for r in runs:
        s = r["spec"]; ov = s.get("overrides", {}); v = next((k for k, o in VAR.items() if o == ov), None)
        if v is None: continue
        kind = "ref200" if s["kind"] == "reference" and s["ncell"] == 200 else ("r%d" % s["stations"] if s["kind"] == "reduced" else None)
        if kind: idx[(s["scn"], v, kind)] = r
    res = {}
    for scn in ("S1", "S2", "S3", "S4", "S5"):
        I = None; rows = {}
        for v in VAR:
            ref = idx.get((scn, v, "ref200"))
            if ref is None: continue
            m = ref["metrics"]; rs = A.series(ref)
            I = I or max(rs["Q1"].max(), rs["Q2"].max()) * 0 + json.load(open(os.path.join(HERE, "competing_outlets_v1.json")))["scenarios"][int(scn[1]) - 1]["inflow_m3_s"]
            row = dict(route=A.route(m), t_act2_s=m["t_act2_s"], t_dom2_s=m["t_dom2_s"], sill_lowering_m=[m["sill_lowering1_m"], m["sill_lowering2_m"]], delivered_m3=[m["delivered1_m3"], m["delivered2_m3"]], heat_to_ice_J=m["heat_to_ice_J"])
            for k in ("r9", "r17"):
                rr = idx.get((scn, v, k))
                if rr:
                    e = A.ser_err(A.series(rr), rs, I)
                    row[k + "_matched_closure_error"] = dict(hydrograph_rmse_over_I=e["hydrograph_rmse_over_I"], t_act2_err_s=None if None in (rr["metrics"]["t_act2_s"], m["t_act2_s"]) else rr["metrics"]["t_act2_s"] - m["t_act2_s"],
                                                           route_match=A.route(rr["metrics"]) == A.route(m), sill_lowering_err_m=[rr["metrics"]["sill_lowering1_m"] - m["sill_lowering1_m"], rr["metrics"]["sill_lowering2_m"] - m["sill_lowering2_m"]])
            rows[v] = row
        if "base" in rows:
            b = rows["base"]
            for v, row in rows.items():
                if v == "base": continue
                row["vs_base"] = dict(route_changed=row["route"] != b["route"], t_act2_shift_s=None if None in (row["t_act2_s"], b["t_act2_s"]) else row["t_act2_s"] - b["t_act2_s"],
                                      t_dom2_shift_s=("none" if row["t_dom2_s"] is None and b["t_dom2_s"] is None else (None if None in (row["t_dom2_s"], b["t_dom2_s"]) else row["t_dom2_s"] - b["t_dom2_s"])),
                                      t_dom2_appears_or_disappears=(row["t_dom2_s"] is None) != (b["t_dom2_s"] is None),
                                      sill_lowering_ratio=[x / max(y, 1e-12) for x, y in zip(row["sill_lowering_m"], b["sill_lowering_m"])], heat_ratio=row["heat_to_ice_J"] / max(b["heat_to_ice_J"], 1e-12))
        res[scn] = rows
    json.dump(res, open(out, "w"), indent=1, default=float); return res


if __name__ == "__main__":
    r = run()
    for scn, rows in r.items():
        for v, row in rows.items(): print(scn, v, row["route"], "t_act2 %.0f" % row["t_act2_s"], "t_dom2", None if row["t_dom2_s"] is None else round(row["t_dom2_s"]), "sill mm", [round(x * 1e3, 1) for x in row["sill_lowering_m"]], row.get("vs_base", {}).get("sill_lowering_ratio"), {k: round(v2["hydrograph_rmse_over_I"], 4) for k, v2 in row.items() if k.endswith("matched_closure_error")})
