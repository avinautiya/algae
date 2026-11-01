"""Verification suite for stefan2DFoam (needs the built OpenFOAM app). python verify_cfd.py WORKDIR [subset]  -> results/verify_cfd.json
Subsets: onedlimit, mass, convergence2d, remap, shield, profile"""
import json, os, sys, time, resource
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); import stefan2d as S, remap as RM
OUT = os.path.join(HERE, "results"); os.makedirs(OUT, exist_ok=True)
def save(res): 
    p = os.path.join(OUT, "verify_cfd.json"); old = json.load(open(p)) if os.path.exists(p) else {}; old.update(res); json.dump(old, open(p, "w"), indent=1, default=float)
def ledgers(r):
    L = r["L"]; hot = L["E_hot_cond"][-1] + L["E_hot_adv"][-1]; return dict(resid_disc_rel=float(L["resid_disc"][-1] / hot), resid_phys_rel_latent=float(L["resid_phys"][-1] / L["E_latent"][-1]), gcl_abs=float(L["gcl_resid"][-1]), gcl_rel=float(L["gcl_resid"][-1] / max(abs(L["Vsw_geom"][-1]), 1e-300)))
def onedlimit(wd):
    out = {}
    for Ste in (0.00631, 0.2):
        for tag, rho_i, jump in (("jump917", S.RI0, True), ("equal", S.RW, False)):
            for nx, dt in ((20, 2.5), (40, 2.5), (80, 2.5), (160, 2.5), (40, 5.0), (40, 10.0)):
                r = S.run(f"{wd}/od_{Ste}_{tag}_{nx}_{dt}", Ste, nx, 4, dt, 3600.0, rho_i=rho_i, jump=jump); L = r["L"]; sex = float(S.neumann_s(r, 3600.0)[0])
                out[f"Ste{Ste}_{tag}_nx{nx}_dt{dt}"] = dict(s_rel_err=float((L["s_mean"][-1] - sex) / sex), row_spread=float(L["s_max"][-1] - L["s_min"][-1]), wall_s=r["wall_s"], **ledgers(r))
    save({"onedlimit": out}); return out
def mass(wd):
    """interface mass balance: (a) liquid mass: rho_w dV_L = rho_w V_in + rho_i V_melt with V_in = sum uw A dt; (b) energy-side swept volume sum(q A dt)/(rho_i L) vs geometry; (c) meshPhi swept volume vs geometry (GCL); (d) face relation u_w.n = (rho_i/rho_w - 1) m_n evaluated from the run's m_n."""
    out = {}
    for tag, rho_i, jump in (("jump917", S.RI0, True), ("equal", S.RW, False)):
        r = S.run(f"{wd}/mass_{tag}", 0.00631, 40, 8, 5.0, 3600.0, rho_i=rho_i, jump=jump, eps=0.3); L = r["L"]; t = np.concatenate([[0], L["t"]]); dt = np.diff(t); Ahot = r["params"]["Ly"] * r["params"]["dz"]
        Vin = float(np.sum(L["uw"] * Ahot * dt)); Vsw = L["Vsw_geom"][-1]; rho_ratio = rho_i / S.RW
        mass_resid = (S.RW * Vsw - (S.RW * Vin + rho_i * Vsw)) / (S.RW * Vsw)          # liquid mass balance: rho_w dV_L - rho_w V_in - rho_i V_melt (V_melt = ice volume lost = swept volume)
        Vsw_energy = L["E_if_cond"][-1] / (rho_i * S.LF)
        mean_mn = Vsw / (Ahot * 3600.0); uwn_expected = (rho_ratio - 1) * mean_mn; uw_mean_run = Vin / (Ahot * 3600.0)
        out[tag] = dict(liquid_mass_resid_rel=float(mass_resid), Vin_over_Vsw=float(Vin / Vsw), expected_Vin_over_Vsw=1 - rho_ratio, Vsw_energy_over_geom=float(Vsw_energy / Vsw), Vsw_meshphi_over_geom=float(L["Vsw_meshphi"][-1] / Vsw), mean_uw_x=float(uw_mean_run), expected_uw_n_times_minus1=float(-uwn_expected), **ledgers(r))
    save({"mass": out}); return out
def convergence2d(wd):
    grid = [(10, 8), (20, 16), (40, 32), (80, 64)]; out = {"space": {}, "time": {}}; prev = None
    for nx, ny in grid:
        r = S.run(f"{wd}/c2_{nx}_{ny}", 0.00631, nx, ny, 2.5, 3600.0, jump=True, eps=0.3); L = r["L"]; prof = r["interface"]; yc = 0.5 * (prof[1:, 0] + prof[:-1, 0]); s_row = 0.5 * (prof[1:, 1] + prof[:-1, 1])
        out["space"][f"{nx}x{ny}"] = dict(cells=r["cells"], Vsw=float(L["Vsw_geom"][-1]), s_mean=float(L["s_mean"][-1]), s_min=float(L["s_min"][-1]), s_max=float(L["s_max"][-1]), wall_s=r["wall_s"], s_at_quarter=float(np.interp(0.0125, prof[:, 0], prof[:, 1])), s_at_mid=float(np.interp(0.025, prof[:, 0], prof[:, 1])), s_at_3q=float(np.interp(0.0375, prof[:, 0], prof[:, 1])), **ledgers(r))
        s_n1d = np.interp(prof[:, 0], r["params"]["Ly"] * (np.arange(ny) + 0.5) / ny, S.neumann_s(r, 3600.0)); out["space"][f"{nx}x{ny}"]["max_dev_from_local_1D"] = float(np.max(np.abs(prof[:, 1] - s_n1d)))
    for dt in (2.5, 5.0, 10.0, 20.0):
        r = S.run(f"{wd}/c2t_{dt}", 0.00631, 40, 32, dt, 3600.0, jump=True, eps=0.3); L = r["L"]; prof = r["interface"]
        out["time"][f"dt{dt}"] = dict(Vsw=float(L["Vsw_geom"][-1]), s_at_mid=float(np.interp(0.025, prof[:, 0], prof[:, 1])), s_max=float(L["s_max"][-1]), s_min=float(L["s_min"][-1]), wall_s=r["wall_s"], **ledgers(r))
    save({"convergence2d": out}); return out
def remap(wd):
    out = {}; r = S.run(f"{wd}/rm_src", 0.00631, 20, 16, 5.0, 3600.0, jump=True, eps=0.3); case = f"{wd}/rm_src"
    txt = open(f"{case}/3600/T").read(); i = txt.index("internalField"); j = txt.index("(", i); n = int(txt[i:j].split()[-1]); body = txt[j + 1:].split("\n)")[0].strip().splitlines(); T = np.array([float(x) for x in body[:n]]).reshape(16, 20)
    prof = r["interface"]; s_node = prof[:, 1]; dy = r["params"]["Ly"] / 16; dz = r["params"]["dz"]; E0 = RM.energy(T, s_node, dy * dz, 20)
    for nxn in (40, 80, 7):
        T2 = RM.remap_rows(T, nxn); E1 = RM.energy(T2, s_node, dy * dz, nxn)
        # accuracy against an independent finer solution of the same problem is not available at identical time; compare with the exact local 1-D profile on the new mesh
        sbar = 0.5 * (s_node[:-1] + s_node[1:]); xi = (np.arange(nxn) + 0.5) / nxn; yc = (np.arange(16) + 0.5) * dy; err = []
        from scipy.special import erf
        for jr in range(16):
            lam, b, tr, Tw = r["lam"][jr], r["b"][jr], r["tRow"][jr] + 3600.0, r["Tw"][jr]; x = xi * sbar[jr]; eta = x / (2 * np.sqrt(S.KAP * tr)); I = lambda e: 0.5 * np.sqrt(np.pi) * (erf(e - b) + erf(b)); Tex = Tw - Tw / I(lam) * I(eta); err.append(np.max(np.abs(T2[jr] - Tex)) / Tw)
        out[f"nx{nxn}"] = dict(energy_rel_change=float((E1 - E0) / E0), max_rel_dev_from_local_1D_profile=float(max(err)))
    # non-conservative linear interpolation for contrast
    for nxn in (7, 13, 33):      # point-sampling (non-conservative) contrast on non-nested meshes
        T3 = np.array([np.interp((np.arange(nxn) + 0.5) / nxn, (np.arange(20) + 0.5) / 20, T[jr]) for jr in range(16)]); out[f"linear_interp_contrast_nx{nxn}"] = dict(energy_rel_change=float((RM.energy(T3, s_node, dy * dz, nxn) - E0) / E0))
    save({"remap": out}); return out
def shield(wd):
    out = {}
    # (a) uniform static shield vs quasi-steady analytic (small Ste so storage is negligible), convergence in nx
    for Ste in (0.0005, 0.00631):
        for h in (30.0, 100.0, 500.0):
            for nx in (20, 40, 80):
                r = S.run(f"{wd}/sh_{Ste}_{h}_{nx}", Ste, nx, 4, 5.0, 3600.0, jump=False, rho_i=S.RW, hShield=[h] * 4); L = r["L"]; Tw = r["Tw"][0]; sex = S.shielded_s(Tw, h, r["params"]["s0"], 3600.0, rho_i=S.RW)
                out[f"Ste{Ste}_h{h}_nx{nx}"] = dict(s_mean=float(L["s_mean"][-1]), s_quasi_steady=sex, rel_dev=float((L["s_mean"][-1] - sex) / sex), **ledgers(r))
    # (b) limits: bare (h huge) equals the unshielded run; h -> 0 gives no melt
    rb = S.run(f"{wd}/sh_bare", 0.00631, 40, 4, 5.0, 3600.0, jump=False, rho_i=S.RW); rh = S.run(f"{wd}/sh_huge", 0.00631, 40, 4, 5.0, 3600.0, jump=False, rho_i=S.RW, hShield=[1e12] * 4); r0 = S.run(f"{wd}/sh_tiny", 0.00631, 40, 4, 5.0, 3600.0, jump=False, rho_i=S.RW, hShield=[1e-6] * 4)
    out["limits"] = dict(bare_vs_h1e12_rel=float((rh["L"]["s_mean"][-1] - rb["L"]["s_mean"][-1]) / rb["L"]["s_mean"][-1]), tiny_h_total_advance_m=float(r0["L"]["s_mean"][-1] - r0["params"]["s0"]), bare_advance_m=float(rb["L"]["s_mean"][-1] - rb["params"]["s0"]))
    # (c) 2-D partial cover: rows y > Ly/2 shielded by a thick debris layer (h = 30), exposed rows bare, same hot wall
    ny = 32; h = np.where((np.arange(ny) + 0.5) / ny > 0.5, 30.0, 1e30)
    for nx, nyy in ((20, 16), (40, 32), (80, 64)):
        hh = np.where((np.arange(nyy) + 0.5) / nyy > 0.5, 30.0, 1e30); r = S.run(f"{wd}/shp_{nx}", 0.00631, nx, nyy, 5.0, 3600.0, jump=True, hShield=hh); L = r["L"]; prof = r["interface"]
        out[f"partial_{nx}x{nyy}"] = dict(s_exposed_at_y0p0125=float(np.interp(0.0125, prof[:, 0], prof[:, 1])), s_covered_at_y0p0375=float(np.interp(0.0375, prof[:, 0], prof[:, 1])), Vsw=float(L["Vsw_geom"][-1]), exposed_advance_ratio_covered_over_exposed=float((np.interp(0.0375, prof[:, 0], prof[:, 1]) - 1e-3) / (np.interp(0.0125, prof[:, 0], prof[:, 1]) - 1e-3)), **ledgers(r))
    save({"shield": out}); return out
def profile(wd):
    out = {}
    for nx, ny, dt in ((40, 32, 5.0), (80, 64, 5.0), (160, 128, 5.0), (80, 64, 2.5)):
        ru0 = resource.getrusage(resource.RUSAGE_CHILDREN); r = S.run(f"{wd}/pf_{nx}_{ny}_{dt}", 0.00631, nx, ny, dt, 3600.0, jump=True, eps=0.3, timeout=7200); ru1 = resource.getrusage(resource.RUSAGE_CHILDREN)
        d = os.popen(f"du -sk {wd}/pf_{nx}_{ny}_{dt}").read().split()[0]; steps = int(3600.0 / dt)
        out[f"{nx}x{ny}_dt{dt}"] = dict(cells=r["cells"], steps=steps, wall_s=r["wall_s"], s_per_cell_step=r["wall_s"] / (r["cells"] * steps), peak_rss_MB_children_max=ru1.ru_maxrss / 1024.0, disk_kB=int(d))
    save({"profile": out}); return out
if __name__ == "__main__":
    wd = sys.argv[1]; subs = sys.argv[2:] or ["onedlimit", "mass", "convergence2d", "remap", "shield", "profile"]
    for s in subs: t = time.time(); globals()[s](wd); print(s, "done", round(time.time() - t, 1), "s", flush=True)
