"""Analysis of the matched-state experiment per FREEZE_TH.md (decision rule frozen). python analyze.py -> analysis.json, figures/th_*.png, RESULTS_TH.md"""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); import baselines as BL
H = 0.01; AL = 1e-6 / 7.0; F = H / AL          # converts raw q*len sums (K m/s * H-lengths) to units alpha/H * H
R = lambda n: json.load(open(os.path.join(HERE, "results", n + ".json"))) if os.path.exists(os.path.join(HERE, "results", n + ".json")) else None
def pick(d, h, dt=10): return d[f"{h}_dt{dt:g}"]
def stats(d, h_pair=("A", "B"), dt=10, plane="-1.0", use_prev=False):
    a, b = pick(d, h_pair[0], dt), pick(d, h_pair[1], dt); s = a["Tb_sec"][plane] / b["Tb_sec"][plane]; out = {"scale": s}
    for X in (10, 30):
        k = f"Q{X}" + ("_prev" if use_prev else ""); qa = a[k] * F; qb = b[k] * F * s; out[f"Q{X}_A"] = qa; out[f"Q{X}_B"] = qb; out[f"delta{X}"] = (qa - qb) / (0.5 * (qa + qb))
    x = np.array(a["x"]); m = (x >= -1e-9) & (x <= 30); qa = np.array(a["q"])[m] * F; qb = np.array(b["q"])[m] * F * s; ln = np.array(a["ln"])[m]; out["e_q"] = float(np.sum(np.abs(qa - qb) * ln) / np.sum(0.5 * (qa + qb) * ln))
    out["Tbout_diff"] = (a["Tb_out"] - b["Tb_out"] * s) / (0.5 * (a["Tb_out"] + b["Tb_out"] * s)); return out
def budget(kind):
    res = {}; base = {n: R(f"{kind}_n{n}_base") for n in (8, 32)}; full = R(f"{kind}_n16_full"); xd = R(f"{kind}_n16_xd28"); xu = R(f"{kind}_n16_xu24")
    if full is None or base[32] is None: return None
    S = {"n8": stats(base[8]) if base[8] else None, "n16": stats(full), "n32": stats(base[32])}; res["mesh"] = S
    for key in ("delta10", "delta30", "e_q"):
        vals = [S[n][key] for n in ("n8", "n16", "n32") if S[n]]; res[f"mesh_seq_{key}"] = vals
    U = {}
    for key in ("delta10", "delta30", "e_q"):
        u_mesh = abs(S["n32"][key] - S["n16"][key]); dts = [stats(full, dt=dt)[key] for dt in (5, 10, 20)]; u_dt = max(dts) - min(dts)
        u_iter = abs(stats(full, use_prev=True)[key] - S["n16"][key]) if key != "e_q" else 0.0
        mp = [stats(full, plane=p)[key] for p in ("-0.5", "-1.0", "-2.0")]; u_match = max(mp) - min(mp)
        u_hist = [stats(full, h_pair=("A", hb))[key] for hb in ("B", "B9", "B13")]; u_locB = max(u_hist) - min(u_hist)
        u_geo = max([abs(stats(v)[key] - S["n16"][key]) for v in (xd, xu) if v] or [0.0])
        U[key] = dict(u_mesh=u_mesh, u_dt=u_dt, u_iter=u_iter, u_match=u_match, u_loc_historyB=u_locB, u_loc_geometry=u_geo, U_frozen=float(np.sqrt(u_mesh ** 2 + u_dt ** 2 + u_iter ** 2 + u_match ** 2 + max(u_locB, u_geo) ** 2)), U_without_B_variants=float(np.sqrt(u_mesh ** 2 + u_dt ** 2 + u_iter ** 2 + u_match ** 2 + u_geo ** 2)),
                      value=S["n32"][key], dts=dts, match_planes=mp, history_B_variants=u_hist)
    res["U"] = U; return res
def main():
    out = {}; st = budget("step"); sr = budget("straight")
    if st is None or sr is None: print("missing results"); return
    out["step"] = st; out["straight"] = sr; dec = {}
    for X in (10, 30):
        k = f"delta{X}"; d1 = st["U"][k]["value"]; U1 = st["U"][k]["U_frozen"]; d2 = sr["U"][k]["value"]; U2 = sr["U"][k]["U_frozen"]; exc = d1 - d2; Ue = float(np.hypot(U1, U2))
        rule = "(1) unresolved/negligible" if abs(d1) <= U1 else ("(2) ordinary developing-flow history effect (straight baseline D reproduces it)" if abs(exc) <= Ue else "(3) step-specific history effect (L and D miss it)")
        dec[f"X{X}H"] = dict(delta_step=d1, U_step=U1, delta_straight=d2, U_straight=U2, excess=exc, U_excess=Ue, decision=rule, L_prediction=0.0, D_prediction=d2, L_error=d1 - 0.0, D_error=d1 - d2)
    out["decision"] = dec
    # Graetz check of the resolved upstream development (n=32 base, A and B at plane -1)
    g = {}
    for kind in ("straight",):
        d = R(f"{kind}_n32_base")
        for h, start in (("A", -3.0), ("B", -11.0)):
            r = pick(d, h); dist = -1.0 - start; Tb_g, qw_g, theta, eta = BL.graetz_section(dist); Tb_r = r["Tb_sec"]["-1.0"]; qw_r = r["qwall_sec"]["-1.0"] * H / AL
            g[h] = dict(distance_H=dist, Tb_resolved=Tb_r, Tb_graetz=Tb_g, rel_Tb=(Tb_r - Tb_g) / Tb_g, qwall_resolved_alpha_over_H_units=qw_r, qwall_graetz=qw_g, rel_q=(qw_r - qw_g) / qw_g)
        a, b = g["A"], g["B"]; g["matched_scale_graetz"] = a["Tb_graetz"] / b["Tb_graetz"]; g["matched_scale_resolved"] = a["Tb_resolved"] / b["Tb_resolved"]
    out["graetz_check"] = g
    # local-closure baseline L (history blind): from the matched bulk temperature at x_c marching 1 H of upstream channel then downstream
    nu, _ = BL.Nu_fd(); L = {}
    for kind in ("step", "straight"):
        d = R(f"{kind}_n32_base"); a = pick(d, "A"); Tb_c = a["Tb_sec"]["-1.0"]; Tb0 = Tb_c * np.exp(-nu / BL.PE * 1.0)
        xs, q, Q = BL.local_closure(Tb0, kind); Qf = lambda X: float(np.interp(X, xs, Q)); res = {}
        for X in (10, 30): qa = pick(d, "A")[f"Q{X}"] * F; b = pick(d, "B"); s = a["Tb_sec"]["-1.0"] / b["Tb_sec"]["-1.0"]; qb = b[f"Q{X}"] * F * s; res[f"Q{X}_L"] = Qf(X); res[f"Q{X}_A"] = qa; res[f"Q{X}_Bstar"] = qb; res[f"rel_err_L_vs_A"] = (Qf(X) - qa) / qa if X == 10 else res.get("rel_err_L_vs_A"); res[f"relerr{X}_L_A"] = (Qf(X) - qa) / qa; res[f"relerr{X}_L_B"] = (Qf(X) - qb) / qb
        L[kind] = res
    out["local_closure"] = L; out["Nu_fd"] = nu
    # cost
    c = {}
    for kind in ("step", "straight"):
        d = R(f"{kind}_n32_base"); c[kind] = dict(cells=pick(d, "A")["cells"], flow_s=pick(d, "A")["tflow"], thermal_s_per_history=pick(d, "A")["wall_s"])
    out["cost"] = c; json.dump(out, open(os.path.join(HERE, "analysis.json"), "w"), indent=1, default=float); print(json.dumps(dec, indent=1, default=float)); print(json.dumps(g, indent=1, default=float))
if __name__ == "__main__": main()
