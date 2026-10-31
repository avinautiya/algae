"""Post-processing of step_vs_ramp cases: wall-heat ledger with face lengths, q(x) comparison, differences vs mesh change. python step_ramp_post.py workdir n1 n2 ..."""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); import step_vs_ramp as S
wd = sys.argv[1]; ns = [int(a) for a in sys.argv[2:]]; out = {}
for n in ns:
    for kind in ("step", "ramp"):
        case = f"{wd}/{kind}_n{n}"; pts, faces, owner, patches = S.mesh_data(case); r = S.analyse(case, n)
        s, c = patches["walls"]; xs = []; qs = []; ln = []
        U = None
        for k in range(s, s + c):
            p = pts[faces[k]]; z0 = [pp[:2] for pp in p if abs(pp[2]) < 1e-12]; ln.append(np.linalg.norm(z0[1] - z0[0]) / S.H); xs.append(0.5 * (z0[0][0] + z0[1][0]) / S.H)
        xs = np.array(xs); o = np.argsort(xs); ln = np.array(ln)[o]; xs = xs[o]; q = np.array(r["q_wall"]); assert len(q) == len(xs)
        # r["q_wall"] was sorted by face-centre x in analyse; same ordering here (centres of straight faces)
        Qw = float(np.sum(q * ln) * S.H * (S.H / n))            # alpha dT/dn integrated over the wall (per unit depth, nondimensional T=1 at inlet)
        Pin = r["Qin"] * 1.0; Pout = r["Qout"] * r["Tb_out"]
        qn = q * S.H / S.ALPHA                       # alpha dT/dn / (alpha dT_in / H)  -> nondimensional wall flux
        m = xs > 0.5
        out[f"{kind}_n{n}"] = dict(cells=r["cells"], Tb_out=r["Tb_out"], reattach=r["reattach_x_over_H"], wall_heat=Qw, Pin=Pin, Pout=Pout, ledger_resid=(Pin - Pout - Qw) / Pin, peak_q=float(qn[m].max()), peak_x=float(xs[m][np.argmax(qn[m])]), mean_q_x0_10=float(np.mean(qn[(xs > 0.5) & (xs < 10)])), x=xs.tolist(), qn=qn.tolist())
        print(kind, n, {k: (round(v, 5) if isinstance(v, float) else v) for k, v in out[f"{kind}_n{n}"].items() if k not in ("x", "qn")})
for n in ns:
    a, b = out[f"step_n{n}"], out[f"ramp_n{n}"]; xg = np.linspace(1, 29, 200); qa = np.interp(xg, a["x"], a["qn"]); qb = np.interp(xg, b["x"], b["qn"])
    out[f"diff_n{n}"] = dict(Tb_out_step_minus_ramp=a["Tb_out"] - b["Tb_out"], wall_heat_step_minus_ramp=a["wall_heat"] - b["wall_heat"], q_rel_L1=float(np.mean(np.abs(qa - qb)) / np.mean(np.abs(qa))), reattach_diff=(a["reattach"] or np.nan) - (b["reattach"] or np.nan)); print("diff", n, out[f"diff_n{n}"])
for kind in ("step", "ramp"):
    for i in range(len(ns) - 1): out[f"meshchange_{kind}_{ns[i]}_{ns[i+1]}"] = dict(dTb=out[f"{kind}_n{ns[i]}"]["Tb_out"] - out[f"{kind}_n{ns[i+1]}"]["Tb_out"], dWall=out[f"{kind}_n{ns[i]}"]["wall_heat"] - out[f"{kind}_n{ns[i+1]}"]["wall_heat"])
json.dump(out, open(os.path.join(HERE, "step_ramp_post.json"), "w"), default=float); print({k: v for k, v in out.items() if k.startswith("meshchange")})
