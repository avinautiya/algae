"""python -m ice_cliff_survival heat-ice-2d --workdir DIR [--Ste 0.00631] [--nx 40 --ny 32 --dt 5 --tend 3600 --eps 0.3 --shield none|wedge --out summary.json]
Component run of the verified 2-D heat + moving-ice solver (rigid fixed liquid layer, no flow, no grains, no free surface). Output is labelled COMPONENT, never environmental."""
import argparse, json, os, shutil, sys
def main(argv):
    p = argparse.ArgumentParser(); p.add_argument("--workdir", required=True); p.add_argument("--Ste", type=float, default=0.00631); p.add_argument("--nx", type=int, default=40); p.add_argument("--ny", type=int, default=32); p.add_argument("--dt", type=float, default=5.0)
    p.add_argument("--tend", type=float, default=3600.0); p.add_argument("--eps", type=float, default=0.0); p.add_argument("--shield", default="none", choices=["none", "wedge"]); p.add_argument("--out", default=None); a = p.parse_args(argv)
    if shutil.which("stefan2DFoam") is None and not os.path.exists("/tmp/claude-0/ofbuild/bin/stefan2DFoam"): sys.stderr.write("stefan2DFoam is not built (see ice_network_v2/foam/BUILD.md; cd ice_cliff_survival/cfd/stefan2DFoam && wmake)\n"); return 2
    from . import stefan2d as S
    from ..coupling import static_shield_rows
    h = None if a.shield == "none" else static_shield_rows(a.ny, 0.05, 0.012, 62.0, 33.0, 0.5)       # illustrative wedge: height 1.2 cm on a 5 cm face slice, k_deb = 0.5 W/m/K (assumed)
    r = S.run(os.path.join(a.workdir, "case"), a.Ste, a.nx, a.ny, a.dt, a.tend, eps=a.eps, hShield=h); L = r["L"]
    out = dict(label="COMPONENT run (rigid liquid layer; no flow, grains or free surface) - not an environmental prediction", params=r["params"], s_mean_m=float(L["s_mean"][-1]), s_min_m=float(L["s_min"][-1]), s_max_m=float(L["s_max"][-1]), swept_volume_m3=float(L["Vsw_geom"][-1]),
               ledger=dict(discrete_residual_over_supply=float(L["resid_disc"][-1] / (L["E_hot_cond"][-1] + L["E_hot_adv"][-1])), sensible_plus_latent_residual_over_latent=float(L["resid_phys"][-1] / L["E_latent"][-1]), gcl_volume_residual=float(L["gcl_resid"][-1])), wall_s=r["wall_s"], shield=a.shield)
    if a.out: json.dump(out, open(a.out, "w"), indent=1)
    print(json.dumps(out, indent=1)); return 0
