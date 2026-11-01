"""M2 correction record for the M1 `heat-ice-2d --shield wedge` example. Re-runs the SAME illustrative case (z_f=1.2 cm, theta=62 deg, beta=33 deg, k_deb=0.5, nx=40, ny=32, dt=5 s, 1 h, Ste=0.00631, eps=0.3)
with (a) the M1 thickness formula (deprecated, reproduces the stored M1 behaviour) and (b) the corrected formula. Writes cfd/results/wedge_correction.json. The M1 results are not edited.
Caveat: k/d remains a 1-D normal effective-resistance approximation; the uniform-conductance verification tests do not validate wedge geometry."""
import json, math, os, sys
from .. import geometry
def rows(fn, ny=32, Ly=0.05, zf=0.012, th=62.0, be=33.0, k=0.5):
    out = []
    for j in range(ny):
        z = ((j + 0.5) / ny * Ly) * math.sin(math.radians(th)); d = fn(z, zf, th, be); out.append((z, d, 1e30 if d <= 0 else k / max(d, 1e-4)))
    return out
def main(workdir):
    from . import stefan2d as S
    old, new = rows(geometry.debris_normal_thickness_v1_projected_horizontal), rows(geometry.debris_normal_thickness); res = dict(profile=dict(z=[r[0] for r in new], d_old=[r[1] for r in old], d_new=[r[1] for r in new], h_old=[r[2] for r in old], h_new=[r[2] for r in new]), runs={})
    for name, rr in (("m1_formula", old), ("m2_corrected", new)):
        r = S.run(os.path.join(workdir, name), 0.00631, 40, 32, 5.0, 3600.0, eps=0.3, hShield=[x[2] for x in rr]); L = r["L"]
        res["runs"][name] = dict(s_mean_m=float(L["s_mean"][-1]), s_min_m=float(L["s_min"][-1]), s_max_m=float(L["s_max"][-1]), discrete_residual_over_supply=float(L["resid_disc"][-1] / (L["E_hot_cond"][-1] + L["E_hot_adv"][-1])))
    res["note"] = "M1 formula = horizontal width x sin(theta) (wrong); M2 = (z_f - z) tan(theta - beta) / sin(theta). Illustrative geometry; k_deb assumed."
    p = os.path.join(os.path.dirname(__file__), "results", "wedge_correction.json"); json.dump(res, open(p, "w"), indent=1); return res
if __name__ == "__main__":
    r = main(sys.argv[1]); print(json.dumps(r["runs"], indent=1)); print("d_old/d_new at row 0:", r["profile"]["d_old"][0], r["profile"]["d_new"][0])
