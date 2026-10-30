"""Reference discretisation error from the mesh sequence 100/200/400(/800), from per-run records.
Differences d_k between successive meshes; observed order p = log2(d_(k-1)/d_k). If the geometric model d_k = d_0 2^(-p k) holds, the error of the mesh-h solution is e_h = d_(h->h/2)/(1 - 2^-p) and of the finer solution e_(h/2) = e_h 2^-p.
Support test for an assumed order: two successive observed orders both in [0.4, 2.5], differing by <= 0.3, differences decreasing. If not supported, both p = 0.5 and p = 1 are reported as a range and the quantity is flagged `order_not_supported`.
python reference_uncertainty.py -> records/reference_uncertainty.json"""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import analyze_campaign as A
CFG = json.load(open(os.path.join(HERE, "competing_outlets_v1.json"))); I_OF = {s["id"]: s["inflow_m3_s"] for s in CFG["scenarios"]}
KEYS = ("t_act2_s", "t_dom2_s", "sill_lowering1_m", "sill_lowering2_m", "delivered1_m3", "delivered2_m3")


def errors(d, p):
    return d / (1 - 2.0 ** -p), d / (2.0 ** p - 1)      # error of the coarser, error of the finer solution of the last pair


def assess(ds):      # ds: successive differences (len 2 or 3)
    ps = [np.log2(ds[i] / ds[i + 1]) if ds[i + 1] > 0 and ds[i] > 0 else None for i in range(len(ds) - 1)]
    supported = bool(len(ps) >= 2 and all(p is not None and 0.4 <= p <= 2.5 for p in ps) and abs(ps[-1] - ps[-2]) <= 0.3)
    out = dict(differences=ds, observed_orders=ps, order_supported=supported)
    if supported: p = ps[-1]; ec, ef = errors(ds[-1], p); out.update(p_used=p, err_coarser_of_last_pair=ec, err_finer_of_last_pair=ef)
    else: out.update(p_used=None, err_finer_range=[ds[-1] / (2.0 ** q - 1) for q in (1.0, 0.5)], err_coarser_range=[ds[-1] / (1 - 2.0 ** -q) for q in (1.0, 0.5)])
    return out


def run():
    res = {}
    for name, lab in (("frozen", lambda s: A.case_set(s) == "frozen"), ("convergence", lambda s: A.case_set(s) == "convergence" and s["kind"] == "reference" and s.get("overrides", {}).get("cfl", 0.4) == 0.4)):
        runs = A.load(lambda r: r["status"] == "completed" and r["spec"]["kind"] == "reference" and lab(r["spec"]) and not r["spec"].get("overrides"))
        by = {}
        for r in runs: by.setdefault(r["spec"]["scn"], {})[r["spec"]["ncell"]] = r
        for scn, m in by.items():
            ns = sorted(m); I = I_OF[scn]
            if len(ns) < 3: continue
            row = dict(meshes=ns)
            hs = [A.ser_err(A.series(m[ns[i]]), A.series(m[ns[i + 1]]), I)["hydrograph_rmse_over_I"] for i in range(len(ns) - 1)]
            row["hydrograph_rmse_over_I"] = assess(hs)
            for k in KEYS:
                v = [m[n]["metrics"][k] for n in ns]
                row[k] = assess([abs(v[i + 1] - v[i]) for i in range(len(v) - 1)]) if all(x is not None for x in v) else None
            res[f"{name}|{scn}"] = row
    json.dump(res, open(os.path.join(HERE, "records", "reference_uncertainty.json"), "w"), indent=1, default=float); return res


if __name__ == "__main__":
    for k, row in run().items():
        print(k, row["meshes"])
        for q in ("hydrograph_rmse_over_I",) + KEYS:
            a = row.get(q)
            if not a: continue
            print("   %-24s d=%s p=%s supported=%s %s" % (q, [float("%.3g" % x) for x in a["differences"]], [None if p is None else round(p, 2) for p in a["observed_orders"]], a["order_supported"],
                  ("e_fine=%.3g e_coarse=%.3g" % (a["err_finer_of_last_pair"], a["err_coarser_of_last_pair"])) if a["order_supported"] else "range e_fine=%s" % [float("%.3g" % x) for x in a["err_finer_range"]]))
