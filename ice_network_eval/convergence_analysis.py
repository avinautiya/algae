"""One-factor-at-a-time convergence about the baseline (r17, dt_geom 450 s, msub 16, rtol 1e-7) and reference mesh/CFL study, from per-run records.
Target = finest reference (800 cells); its own uncertainty is reported as the 400->800 change (first-order scheme: error of the 800-cell run is of that size)."""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import analyze_campaign as A
BASE = dict(stations=17, dt_geom=450.0, msub=16, rtol=1e-7)


def run(out=os.path.join(HERE, "records", "convergence_summary.json")):
    runs = A.load(lambda r: r["status"] == "completed" and A.case_set(r["spec"]) == "convergence")
    res = {}
    cfg = json.load(open(os.path.join(HERE, "competing_outlets_v1.json")))
    for scn in ("S1", "S2", "S3", "S4"):
        I = next(s for s in cfg["scenarios"] if s["id"] == scn)["inflow_m3_s"]
        refs = {}; reds = []
        for r in runs:
            s = r["spec"]
            if s["scn"] != scn: continue
            if s["kind"] == "reference": refs[(s["ncell"], s.get("overrides", {}).get("cfl", 0.4))] = r
            else: reds.append(r)
        if (800, 0.4) not in refs: continue
        tgt = refs[(800, 0.4)]; ts = A.series(tgt); tm = tgt["metrics"]
        ref_tab = {}
        for (n, cfl), r in sorted(refs.items()):
            e = A.ser_err(A.series(r), ts, I) if (n, cfl) != (800, 0.4) else None
            ref_tab[f"{n}_cfl{cfl}"] = dict(hydrograph_rmse_over_I_vs800=None if e is None else e["hydrograph_rmse_over_I"], solver_s=r["timing"]["solver_s"], nsteps=r.get("nsteps"),
                                             **{k: r["metrics"][k] for k in ("t_act2_s", "t_dom2_s", "sill_lowering1_m", "sill_lowering2_m", "delivered1_m3", "delivered2_m3")})
        sweeps = {}
        for r in reds:
            s = r["spec"]; ov = s.get("overrides", {}); par = dict(stations=s["stations"], dt_geom=float(s["dt_geom"]), msub=ov.get("msub", 16), rtol=ov.get("rtol", 1e-7))
            diff = [k for k in BASE if par[k] != BASE[k]]
            factor = "baseline" if not diff else diff[0]
            if len(diff) > 1: continue
            e = A.ser_err(A.series(r), ts, I); m = r["metrics"]
            ev = lambda k: None if (m[k] is None or tm[k] is None) else m[k] - tm[k]
            sweeps.setdefault(factor, []).append(dict(value=par[factor] if factor != "baseline" else None, hydrograph_rmse_over_I=e["hydrograph_rmse_over_I"], cumvol_err_rel_total=e["cumvol_err_rel_total"], t_act2_err_s=ev("t_act2_s"), t_dom2_err_s=ev("t_dom2_s"),
                                                      sill_err_m=[m["sill_lowering1_m"] - tm["sill_lowering1_m"], m["sill_lowering2_m"] - tm["sill_lowering2_m"]], delivered_err_m3=[m["delivered1_m3"] - tm["delivered1_m3"], m["delivered2_m3"] - tm["delivered2_m3"]],
                                                      solver_s=r["timing"]["solver_s"], law_evals=r.get("discharge_law_evaluations"), conservation=r.get("conservation")))
        base = sweeps.pop("baseline", [None])[0]
        for f, rows in sweeps.items():
            rows.sort(key=lambda x: x["value"])
            if base: rows.append(dict(base, value=BASE[f])); rows.sort(key=lambda x: x["value"])
            for a, b in zip(rows[:-1], rows[1:]):     # observed order of the hydrograph error between successive refinements (ratio of the parameter)
                ratio = max(b["value"], a["value"]) / min(a["value"], b["value"]); ea, eb = a["hydrograph_rmse_over_I"], b["hydrograph_rmse_over_I"]
                fine, coarse = (b, a) if (f in ("stations", "msub")) else ((a, b) if f in ("dt_geom", "rtol") else (b, a))
                b["observed_order_vs_previous"] = None if min(ea, eb) <= 0 else float(np.log(max(ea, eb) / min(ea, eb)) / np.log(ratio))
        res[scn] = dict(reference=ref_tab, reference_uncertainty_400_to_800=dict(
            hydrograph_rmse_over_I=ref_tab["400_cfl0.4"]["hydrograph_rmse_over_I_vs800"], **{k: abs(ref_tab["800_cfl0.4"][k] - ref_tab["400_cfl0.4"][k]) if ref_tab["800_cfl0.4"][k] is not None and ref_tab["400_cfl0.4"][k] is not None else None for k in ("t_act2_s", "t_dom2_s", "sill_lowering1_m", "sill_lowering2_m")}),
            sweeps=sweeps, baseline=base)
    json.dump(res, open(out, "w"), indent=1, default=float); return res


if __name__ == "__main__":
    r = run()
    for scn, v in r.items():
        print(scn, "ref unc", {k: (round(x, 5) if x else x) for k, x in v["reference_uncertainty_400_to_800"].items()})
        for f, rows in v["sweeps"].items(): print("  ", f, [(x["value"], round(x["hydrograph_rmse_over_I"], 4), None if x.get("observed_order_vs_previous") is None else round(x["observed_order_vs_previous"], 2), round(x["solver_s"], 1)) for x in rows])
