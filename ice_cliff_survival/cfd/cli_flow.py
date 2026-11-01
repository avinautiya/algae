"""python -m ice_cliff_survival flow-stefan --workdir DIR [--rho-i 917] [--Tin 2] [--Ubar 0.005] [--tend 1800] [--dt 5] [--dx 0.004] [--ny 24] [--shield] [--out summary.json]
Executable first flow-heat-moving-ice coupling (single-phase rigid-lid channel, one melting bare/shielded ice floor). COMPONENT run: assumed inputs, laminar, no grains/free surface/atmosphere; never an environmental prediction.
--rho-i 1000 = equal-density stage A; --rho-i 917 = stage B with the LOCAL interface mass condition. Aborts (explicit error) if the floor slope limit is exceeded."""
import argparse, json, os, sys
def main(argv):
    p = argparse.ArgumentParser(prog="flow-stefan"); p.add_argument("--workdir", required=True); p.add_argument("--rho-i", type=float, default=917.0); p.add_argument("--Tin", type=float, default=2.0); p.add_argument("--Ubar", type=float, default=0.005)
    p.add_argument("--tend", type=float, default=1800.0); p.add_argument("--dt", type=float, default=5.0); p.add_argument("--dx", type=float, default=0.004); p.add_argument("--ny", type=int, default=24); p.add_argument("--shield", action="store_true", help="static corrected-wedge shield on the middle half of the ice section (equal density only)")
    p.add_argument("--out", default=None); a = p.parse_args(argv)
    if not os.path.exists("/tmp/claude-0/ofbuild/bin/flowStefanFoam"): sys.stderr.write("flowStefanFoam is not built (see ice_network_v2/foam/BUILD.md; cd ice_cliff_survival/cfd/flowStefanFoam && wmake)\n"); return 2
    from . import flowstefan as F
    kw = dict(H=0.02, dz=0.01, Lu=0.08, Li=0.2, Ld=0.06, dx=a.dx, ny=a.ny, Ubar=a.Ubar, Tin=a.Tin, dt=a.dt, tend=a.tend, rho_i=a.rho_i)
    if a.shield:
        if a.rho_i != F.RHO_W: sys.stderr.write("UnsupportedRun: shielded face with a density jump would need the enthalpy of the liquid drawn through the shield; use --rho-i 1000\n"); return 2
        from ..geometry import debris_normal_thickness
        nx = round(0.34 / a.dx); Dx = 0.34 / nx; h = 0.5 / debris_normal_thickness(0.003, 0.012, 62.0, 33.0); kw["hShield"] = [h if (0.08 + 0.05 <= (i + 0.5) * Dx <= 0.08 + 0.15) else 1e30 for i in range(nx)]
    R = F.run_case(os.path.join(a.workdir, "case"), **kw); L = R["L"]
    out = dict(label="COMPONENT run (assumed inputs; laminar rigid-lid channel; no grains/free surface/atmosphere) - not an environmental prediction", params={k: v for k, v in R["P"].items() if k != "hShield"}, cells=R["cells"], steps=len(L), ice_volume_melted_m3=float(L["Vsw"][-1]),
               ledger=dict(discrete_residual_over_terms=float(L["R_disc"][-1] / (L["Ein"][-1] + abs(L["Efl"][-1]) + L["Eif"][-1])), physical_residual_over_terms=float(L["R_phys"][-1] / (L["Ein"][-1] + abs(L["Efl"][-1]) + L["Eif"][-1])), interface_heat_J=float(L["Eif"][-1]), latent_J=float(L["Elat"][-1]),
                           local_interface_mass_residual=float(L["massLocal"].max()), gcl_residual_m3=float(abs(L["gcl"]).max()), coupling_defect_max_m=float(L["coupRes"].max()), mean_coupling_iters=float(L["nCoup"].mean())), wall_s=R["wall_s"], cpu_s=R["cpu_s"], maxrss_MB=R["maxrss_MB"])
    if a.out: json.dump(out, open(a.out, "w"), indent=1)
    print(json.dumps(out, indent=1)); return 0
