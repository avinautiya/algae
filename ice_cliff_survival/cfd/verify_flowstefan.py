"""Verification ladder for flowStefanFoam (M2). `python -m ice_cliff_survival.cfd.verify_flowstefan WORKDIR [rung ...]` runs rungs (default all) and writes cfd/results/m2/<rung>.json.
rungs: r1_graetz  r2_stefan  r3_frozen  r4_coupled  r5_localmass  r6_refine  r7_domain  r8_shield   (all inputs are ASSUMED test parameters, not measurements)"""
import json, math, os, sys, time
import numpy as np
from scipy.linalg import eigh
from . import flowstefan as F, flowstefan_post as Pp
from .stefan2d import solve_lambda
RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", "m2"); KAP = F.K_W / (F.RHO_W * F.CP_W)
BASE = dict(H=0.02, dz=0.01, Lu=0.08, Li=0.2, Ld=0.06, dx=0.004, ny=24, Ubar=0.005, Tin=2.0, dt=5.0, tend=1800.0)
def save(name, d):
    os.makedirs(RES, exist_ok=True); json.dump(d, open(os.path.join(RES, name + ".json"), "w"), indent=1, default=float)
def lastdump(case, R):
    st = Pp.read_dump(f"{case}/dump.csv", R["nx"], R["P"]["ny"]); return st
# ---------------------------------------------------------------- Graetz-Poiseuille reference (parallel plates, both walls and lid T=0, inlet T=T_in, u = 6 Ub eta (1-eta)); axial conduction O(1/Pe^2) neglected
def graetz_modes(N=2000):
    h = 1.0 / N; eta = np.arange(1, N) * h; W = 6 * eta * (1 - eta)
    A = (2 * np.eye(N - 1) - np.eye(N - 1, k=1) - np.eye(N - 1, k=-1)) / h ** 2
    lam2, V = eigh(A, np.diag(W)); V = V / np.sign(V[np.argmax(np.abs(V), axis=0), np.arange(V.shape[1])])
    nm = 60; lam2 = lam2[:nm]; V = V[:, :nm]; C = (W[:, None] * V).sum(0) * h / ((W[:, None] * V ** 2).sum(0) * h); a = (W[:, None] * V).sum(0) * h; dphi0 = (4 * V[0] - V[1]) / (2 * h)
    return lam2, C, a, dphi0
def graetz_ref(x, Ub, H, Tin, modes):
    lam2, C, a, d0 = modes; xi = x / (H * Ub * H / KAP); e = np.exp(-np.outer(xi, lam2))
    return dict(Tb=Tin * (e * (C * a)).sum(1), qwall=F.K_W * Tin / H * (e * (C * d0)).sum(1))      # qwall: per-wall flux magnitude
def r1_graetz(wd, nys=(12, 24, 48)):
    modes = graetz_modes(); out = []; Ub = 0.001
    for ny in nys:
        dx = 0.004 * 12 / ny * 2 if False else 0.28 / (70 * ny // 24); nx = round(0.28 / dx); dx = 0.28 / nx
        R = F.run_case(f"{wd}/r1_ny{ny}", H=0.02, dz=0.01, Lu=0.0, Ld=0.0, Li=0.28, dx=dx, ny=ny, Ubar=Ub, Tin=1.0, Tlid_fixed=0.0, melt=False, dt=10.0, tend=2400.0, dump=True, dumpEvery=10 ** 6)
        st = lastdump(f"{wd}/r1_ny{ny}", R)[-1]; Dx = 0.28 / R["nx"]; xc = (np.arange(R["nx"]) + 0.5) * Dx; Tb_n = (st["Ux"] * st["T"]).sum(1) / st["Ux"].sum(1)
        q = np.zeros(R["nx"]); q[st["F"][:, 0].astype(int)] = st["F"][:, 5]; sel = (xc > 0.1) & (xc < 0.26); ref = graetz_ref(xc, Ub, 0.02, 1.0, modes)
        eta = (np.arange(ny) + 0.5) / ny; i9 = int(0.9 * R["nx"]); u = st["Ux"][i9]; ub = u.mean(); eU = float(np.sqrt(np.mean((u - 6 * ub * eta * (1 - eta)) ** 2)) / ub)
        out.append(dict(ny=ny, nx=R["nx"], err_Tb_rel_max=float(np.max(np.abs(Tb_n[sel] / ref["Tb"][sel] - 1))), err_qwall_rel_max=float(np.max(np.abs(q[sel] / ref["qwall"][sel] - 1))), err_qwall_rel_at_0p9L=float(abs(q[i9] / ref["qwall"][i9] - 1)),
                        err_u_rms_over_Ub=eU, steady_dEsens_over_Ein_rate=float(abs(R['L']['Esens'][-1] - R['L']['Esens'][-2]) / (R['L']['Ein'][-1] - R['L']['Ein'][-2])), wall_s=R["wall_s"]))
    d = dict(rung="r1 no melting: laminar Graetz-Poiseuille (Pe_H=%.0f, axial conduction neglected in reference)" % (Ub * 0.02 / KAP), results=out); save("r1_graetz", d); return d
# ---------------------------------------------------------------- no imposed flow: Stefan limit
def r2_one(wd, rho_i, ny, dt, tend=3600.0, Ste=0.00631, s0=1e-3, tag=""):
    r = rho_i / F.RHO_W; Tw = Ste * F.LF / F.CP_W; lam, b = solve_lambda(Ste, r, True); t0 = (s0 / (2 * lam)) ** 2 / KAP
    R = F.run_case(f"{wd}/r2_{tag}", H=s0, dz=0.01, Lu=0.0, Ld=0.0, Li=0.02, dx=0.004, ny=ny, Ubar=0.0, Tin=0.0, closedEnds=True, lidClosure=(r != 1.0), Tlid_fixed=Tw, rho_i=rho_i, dt=dt, tend=tend, neumann=(Tw, lam, b, t0))
    L = R["L"]; sn = s0 + L["Vsw"] / (0.02 * 0.01); sa = 2 * lam * np.sqrt(KAP * (L["t"] + t0)); e = (sn - sa) / sa
    return dict(rho_i=rho_i, ny=ny, dt=dt, s_num=float(sn[-1]), s_exact=float(sa[-1]), rel_err_end=float(e[-1]), rel_err_max=float(np.max(np.abs(e))), R_disc=float(L["R_disc"][-1]), R_phys=float(L["R_phys"][-1]), massLocal=float(L["massLocal"].max()), gcl=float(np.abs(L["gcl"]).max()), umax=float(L["umax"].max()), lam=lam, b=b)
def r2_stefan(wd):
    res = []
    for rho in (1000.0, 917.0):
        for ny, dt in ((20, 40.0), (40, 40.0), (80, 40.0), (40, 20.0), (40, 10.0)): res.append(r2_one(wd, rho, ny, dt, tag=f"{int(rho)}_{ny}_{int(dt)}"))
    d = dict(rung="r2 no imposed flow: Stefan/Neumann limit (1-D, equal density and density jump with LOCAL interface BC; lid inflow closure is a TEST-ONLY lid condition)", results=res); save("r2_stefan", d); return d
# ---------------------------------------------------------------- vanishing latent motion -> fixed geometry
def r3_frozen(wd):
    kw = dict(H=0.02, dz=0.01, Lu=0.0, Ld=0.0, Li=0.28, dx=0.004, ny=24, Ubar=0.001, Tin=1.0, Tlid_fixed=0.0, dt=10.0, tend=1200.0, dump=True, dumpEvery=10 ** 6)
    A = F.run_case(f"{wd}/r3_fixed", melt=False, **kw); B = F.run_case(f"{wd}/r3_slow", melt=True, Lmult=1e6, rho_i=F.RHO_W, **kw)
    sa, sb = lastdump(f"{wd}/r3_fixed", A)[-1], lastdump(f"{wd}/r3_slow", B)[-1]; qa = sa["F"][:, 5]; qb = sb["F"][:, 5]
    d = dict(rung="r3 vanishing latent motion (L x 1e6) vs fixed geometry", max_rel_diff_wall_flux=float(np.max(np.abs(qb / qa - 1))), max_rel_diff_T=float(np.max(np.abs(sb["T"] - sa["T"])) / np.max(sa["T"])), max_floor_displacement_m=float(np.abs(sb["b"]).max()),
             Eif_fixed=float(A["L"]["Eif"][-1]), Eif_slow=float(B["L"]["Eif"][-1]), R_disc_slow=float(B["L"]["R_disc"][-1])); save("r3_frozen", d); return d
# ---------------------------------------------------------------- coupled throughflow: ledgers (solver) + independent reconstruction from saved fields
def budget_stats(case, R, rho):
    st = lastdump(case, R); B = Pp.budgets(st, R["P"], R["nx"], R["Ltot"], R["P"]["Lu"], R["P"]["Li"], rho)
    Rp = np.array([b["R_phys"] for b in B]); sc = np.array([b["scale"] for b in B]); ml = np.array([b["mloc"] for b in B]); mg = np.array([b["mglob"] for b in B])
    return dict(indep_energy_step_resid_max_rel=float(np.max(np.abs(Rp) / sc)), indep_energy_cum_resid_rel=float(abs(Rp.sum()) / sc.sum()), indep_local_mass_resid_max_rel=float(ml.max()), indep_global_mass_resid_rel_max=float(np.abs(mg).max()), n_steps=len(B))
def summarize(R, rho):
    L = R["L"]; Eref = float(L["Ein"][-1] + abs(L["Efl"][-1]) + L["Eif"][-1] + 1e-300); lat = float(L["Elat"][-1])
    return dict(rho_i=rho, cells=R["cells"], steps=len(L), mean_coupling_iters=float(L["nCoup"].mean()), max_coupling_iters=int(L["nCoup"].max()), final_coupling_defect_m=float(L["coupRes"][-1]), max_coupling_defect_m=float(L["coupRes"].max()),
                Vsw=float(L["Vsw"][-1]), E_if=float(L["Eif"][-1]), E_latent=lat, R_disc_over_scale=float(L["R_disc"][-1] / Eref), R_phys_over_scale=float(L["R_phys"][-1] / Eref), interface_heat_vs_latent_rel=float((L["Eif"][-1] - lat) / lat),
                gcl_max_over_Vsw=float(np.abs(L["gcl"]).max() / max(abs(L["Vsw"][-1]), 1e-300)), massLocal_solver=float(L["massLocal"].max()), mass_net_Qout_minus_Qin=float(L["Qout"][-1] - L["Qin"][-1]), wall_s=R["wall_s"], cpu_s=R["cpu_s"], maxrss_MB=R["maxrss_MB"])
def r4_coupled(wd):
    res = []
    for rho in (1000.0, 917.0):
        R = F.run_case(f"{wd}/r4_{int(rho)}", rho_i=rho, dump=True, **BASE); s = summarize(R, rho); s.update(budget_stats(f"{wd}/r4_{int(rho)}", R, rho)); res.append(s)
        st = lastdump(f"{wd}/r4_{int(rho)}", R)[-1]; os.makedirs(RES, exist_ok=True)
        np.savez(os.path.join(RES, f"fields_r4_{int(rho)}.npz"), T=st["T"], Ux=st["Ux"], Uy=st["Uy"], b=st["b"], q=st["F"][:, 5], t=st["t"], Dx=R["Ltot"] / R["nx"], H=R["P"]["H"], t_series=R["L"]["t"], Vsw_series=R["L"]["Vsw"], Eif_series=R["L"]["Eif"], Elat_series=R["L"]["Elat"])
        s["floor_displacement_mm"] = [float(v * 1e3) for v in st["b"][::5]]; s["wall_flux_W_m2_along_ice"] = [float(v) for v in st["F"][::5, 5]]; s["max_slope"] = float(np.max(np.abs(np.diff(st["b"]))) / (BASE["Lu"] + BASE["Li"] + BASE["Ld"]) * R["nx"])
    d = dict(rung="r4 coupled throughflow flow-heat-moving-ice (bare ice), equal density and density jump", results=res); save("r4_coupled", d); return d
def r5_localmass(wd):
    R = F.run_case(f"{wd}/r5", rho_i=917.0, dump=True, dumpEvery=1, **{**BASE, "tend": 900.0}); st = lastdump(f"{wd}/r5", R); r = 917.0 / F.RHO_W; out = dict(rung="r5 local interface mass balance with non-uniform melting")
    a, b = st[-2], st[-1]; dt = b["t"] - a["t"]; Dx = R["Ltot"] / R["nx"]; dz = R["P"]["dz"]; F_ = b["F"]; col = F_[:, 0].astype(int); db = b["b"][col + 1] - b["b"][col]; phib = dz * (F_[:, 1] * db - F_[:, 2] * Dx)
    Vf = dz * Dx * 0.5 * ((a["b"][col] - b["b"][col]) + (a["b"][col + 1] - b["b"][col + 1])); tgt = (1 - r) * Vf / dt
    out.update(melt_speed_ratio_max_over_min=float(np.max(Vf) / np.min(Vf)), local_resid_over_max_target=float(np.max(np.abs(phib - tgt)) / np.max(np.abs(tgt))), uniform_closure_local_error_over_max_target=float(np.max(np.abs(tgt.mean() - tgt)) / np.max(np.abs(tgt))),
               note="dt between saved states = %g s; uniform_closure_* is what a single global uniform interface velocity would have mis-assigned (not used)" % dt); save("r5_localmass", out); return out
def qoi(R, st):
    s = st[-1]; Dx = R["Ltot"] / R["nx"]; P = R["P"]; x = (np.arange(R["nx"] + 1)) * Dx; o = s["O"]; Uo = o[:, 0]; To = o[:, 2]
    f = lambda xx: float(np.interp(xx, x, s["b"]))
    return dict(Vsw=float(R["L"]["Vsw"][-1]), b_at=[f(P["Lu"] + a * P["Li"]) for a in (0.1, 0.25, 0.5, 0.75, 0.95)], T_out=float((Uo * To).sum() / Uo.sum()), E_if=float(R["L"]["Eif"][-1]))
def order(f0, f1, f2): return float(math.log2(abs(f0 - f1) / abs(f1 - f2))) if abs(f1 - f2) > 0 and abs(f0 - f1) > 0 else float("nan")
def r6_refine(wd):
    out = dict(rung="r6 separate refinement: flow/thermal mesh, physical timestep, coupling tolerance")
    meshes = ((0.008, 12), (0.004, 24), (0.002, 48)); m = []
    for dx, ny in meshes:
        R = F.run_case(f"{wd}/r6_m{ny}", rho_i=917.0, dump=True, dumpEvery=10 ** 6, **{**BASE, "dx": dx, "ny": ny, "dt": 2.5}); q = qoi(R, lastdump(f"{wd}/r6_m{ny}", R)); q.update(dx=dx, ny=ny, cells=R["cells"], wall_s=R["wall_s"], cpu_s=R["cpu_s"], R_phys=float(R["L"]["R_phys"][-1]), R_disc=float(R["L"]["R_disc"][-1])); m.append(q)
    out["space"] = m; out["space_order_Vsw"] = order(m[0]["Vsw"], m[1]["Vsw"], m[2]["Vsw"]); out["space_order_b_mid"] = order(m[0]["b_at"][2], m[1]["b_at"][2], m[2]["b_at"][2]); out["space_order_T_out"] = order(m[0]["T_out"], m[1]["T_out"], m[2]["T_out"])
    tt = []
    for dt in (20.0, 10.0, 5.0, 2.5):
        R = F.run_case(f"{wd}/r6_dt{dt}", rho_i=917.0, dump=True, dumpEvery=10 ** 6, **{**BASE, "dt": dt}); q = qoi(R, lastdump(f"{wd}/r6_dt{dt}", R)); q.update(dt=dt, wall_s=R["wall_s"], R_phys=float(R["L"]["R_phys"][-1])); tt.append(q)
    out["time"] = tt; out["time_order_Vsw"] = order(tt[0]["Vsw"], tt[1]["Vsw"], tt[2]["Vsw"]); out["time_order_Vsw_fine"] = order(tt[1]["Vsw"], tt[2]["Vsw"], tt[3]["Vsw"]); out["time_order_b_mid"] = order(tt[1]["b_at"][2], tt[2]["b_at"][2], tt[3]["b_at"][2])
    cc = []
    for mc, tol in ((1, 1e-6), (12, 1e-2), (12, 1e-4), (12, 1e-6), (12, 1e-9)):
        R = F.run_case(f"{wd}/r6_c{mc}_{tol}", rho_i=917.0, dump=True, dumpEvery=10 ** 6, maxCoupling=mc, tolCoupling=tol, **BASE); q = qoi(R, lastdump(f"{wd}/r6_c{mc}_{tol}", R)); L = R["L"]
        q.update(maxCoupling=mc, tol=tol, mean_iters=float(L["nCoup"].mean()), R_phys_rel=float(L["R_phys"][-1] / L["Eif"][-1]), coupling_defect_max_m=float(L["coupRes"].max()), wall_s=R["wall_s"]); cc.append(q)
    out["coupling"] = cc; save("r6_refine", out); return out
def r7_domain(wd):
    out = dict(rung="r7 domain length / outlet influence (QoI on the ice section)"); res = []
    for Lu, Ld in ((0.08, 0.06), (0.16, 0.06), (0.08, 0.12), (0.04, 0.03), (0.08, 0.02)):
        R = F.run_case(f"{wd}/r7_{Lu}_{Ld}", rho_i=917.0, dump=True, dumpEvery=10 ** 6, **{**BASE, "Lu": Lu, "Ld": Ld}); q = qoi(R, lastdump(f"{wd}/r7_{Lu}_{Ld}", R)); q.update(Lu=Lu, Ld=Ld); res.append(q)
    out["cases"] = res; save("r7_domain", out); return out
def r8_shield(wd):
    from ..geometry import debris_normal_thickness
    P = {**BASE}; nx = round((P["Lu"] + P["Li"] + P["Ld"]) / P["dx"]); Dx = (P["Lu"] + P["Li"] + P["Ld"]) / nx; d = debris_normal_thickness(0.003, 0.012, 62.0, 33.0); h = 0.5 / d
    hs = [h if (P["Lu"] + 0.25 * P["Li"] <= (i + 0.5) * Dx <= P["Lu"] + 0.75 * P["Li"]) else 1e30 for i in range(nx)]
    A = F.run_case(f"{wd}/r8_bare", rho_i=F.RHO_W, dump=True, dumpEvery=10 ** 6, **P); B = F.run_case(f"{wd}/r8_shield", rho_i=F.RHO_W, dump=True, dumpEvery=10 ** 6, hShield=hs, **P)
    sa, sb = summarize(A, F.RHO_W), summarize(B, F.RHO_W); sb["h_shield_W_m2K"] = h; sb["d_normal_m"] = d
    out = dict(rung="r8 static debris shield (corrected thickness, middle half of ice section, equal density), vs bare", bare=sa, shield=sb, melt_volume_ratio=sb["Vsw"] / sa["Vsw"], note="uniform h over the covered stretch; massless debris; k_deb assumed 0.5; not a wedge-geometry validation"); save("r8_shield", out); return out
def r9_profile(wd):
    """measured cost of the first coupled case, run ALONE on an otherwise idle machine (record: cells, steps, coupling iterations, flow/pressure solves, wall+CPU time, peak RSS, output size)."""
    import subprocess
    res = []
    for rho, dump in ((917.0, False), (917.0, True)):
        R = F.run_case(f"{wd}/r9_{int(dump)}", rho_i=rho, dump=dump, writeEvery=BASE["tend"], **BASE); L = R["L"]; iters = int(L["nCoup"].sum()); size = int(subprocess.run(["du", "-sk", f"{wd}/r9_{int(dump)}"], capture_output=True, text=True).stdout.split()[0])
        res.append(dict(dump=dump, cells=R["cells"], steps=len(L), coupling_iterations_total=iters, momentum_solves=0, pressure_solves=3 * iters, temperature_solves=iters, wall_s=R["wall_s"], cpu_s=R["cpu_s"], maxrss_MB=R["maxrss_MB"], output_kB=size,
                        us_per_cell_per_coupling_iteration=R["wall_s"] / (R["cells"] * iters) * 1e6, us_per_cell_per_step=R["wall_s"] / (R["cells"] * len(L)) * 1e6))
    d = dict(rung="r9 measured profile of the first coupled case (single process, machine otherwise idle)", results=res, base=BASE); save("r9_profile", d); return d
RUNGS = dict(r1_graetz=r1_graetz, r2_stefan=r2_stefan, r3_frozen=r3_frozen, r4_coupled=r4_coupled, r5_localmass=r5_localmass, r6_refine=r6_refine, r7_domain=r7_domain, r8_shield=r8_shield, r9_profile=r9_profile)
if __name__ == "__main__":
    wd = sys.argv[1]; names = sys.argv[2:] or list(RUNGS)
    for n in names:
        t = time.time()
        try: RUNGS[n](wd); print(n, "ok", round(time.time() - t, 1), "s", flush=True)
        except Exception as e:
            import traceback; save(n + "_FAILED", dict(rung=n, error=repr(e), trace=traceback.format_exc())); print(n, "FAILED", repr(e), flush=True)
