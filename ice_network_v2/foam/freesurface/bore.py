"""B. Hydraulic jump as a translating BORE (Galilean image of the stationary jump), 2-D, stock interFoam, laminar (nu_w = 1e-6), sigma = 0.07.
WHY A BORE: the stationary-jump version (fixed outlet water level via a per-face p_rgh profile, kept in jump.py / jump_stationary_outlet_failed.log) FAILED: a fixed-value p_rgh profile that does not follow the actual free surface drives
tens-to-hundreds of m/s of air through cells that momentarily change phase at the outlet (max |U| 114-155 m/s in the air, dt collapsing to 1e-5 s). Preserved as a failed case. The bore needs no outlet level condition.
Setup: ahead of the bore water depth h0 = 0.05 m at rest; behind it depth h_b = h0 (sqrt(1+8 Fr^2) - 1)/2 with Fr = c/sqrt(g h0) = 2.5 (c = 1.751 m/s), moving at u_p = c (1 - h0/h_b) = 1.181 m/s (mass + momentum jump conditions in the bore frame).
Bore frame = stationary jump with upstream Froude number Fr1 = Fr = 2.5 (relative velocity c over depth h0), i.e. Re_h0 = c h0/nu = 8.8e4, We = 2.2e3 (same dimensionless numbers as the stationary case).
BCs: left inlet U = u_p over the water depth (variableHeightFlowRateInletVelocity, flow rate u_p h_b dz) with alpha prescribed (water for y < h_b); right end closed wall; top atmosphere totalPressure p_rgh 0; bed noSlip. L = 4 m, height 0.3 m, bore starts at x = 0.5 m; t_end = 1.5 s.
Measured: bore speed from the depth-profile front, depth behind (h_b,meas), mass and momentum jump conditions in the bore frame using measured c, h0, h_b, u; head loss vs (h_b-h0)^3/(4 h_b h0); resolved viscous dissipation fraction. dy = h0/n, n = 4, 6, 8, 12, and a time-step check (maxCo halved) at n = 8."""
import json, os, shutil, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); import common as C, jump as J
G = 9.81; H0 = 0.05; FR = 2.5; CB = FR * np.sqrt(G * H0); HB = 0.5 * H0 * (np.sqrt(1 + 8 * FR ** 2) - 1); UP = CB * (1 - H0 / HB); L = 4.0; HT = 0.3; X0 = 0.5
def make(case, n, tend=1.5, maxco=0.3):
    shutil.rmtree(case, ignore_errors=True); dy = H0 / n; nx = int(round(L / dy)); ny = int(round(HT / dy)); dz = dy
    C.write(f"{case}/system/blockMeshDict", "dictionary", "blockMeshDict", f"convertToMeters 1; vertices ((0 0 0) ({L} 0 0) ({L} {HT} 0) (0 {HT} 0) (0 0 {dz}) ({L} 0 {dz}) ({L} {HT} {dz}) (0 {HT} {dz})); blocks (hex (0 1 2 3 4 5 6 7) ({nx} {ny} 1) simpleGrading (1 1 1)); boundary (inlet {{type patch; faces ((0 4 7 3));}} wall {{type wall; faces ((1 2 6 5) (0 1 5 4));}} atmosphere {{type patch; faces ((3 7 6 2));}} frontAndBack {{type empty; faces ((0 3 2 1) (4 5 6 7));}});\n")
    C.physics(case); C.control(case, tend, 1e-4, maxco, 0.25, maxdt=0.005)
    yc = (np.arange(ny) + 0.5) * dy; xc = (np.arange(nx) + 0.5) * dy
    fmt = lambda v: "nonuniform List<scalar> %d(%s)" % (len(v), " ".join("%.10g" % x for x in v))
    C.write(f"{case}/0/alpha.water.orig", "volScalarField", "alpha.water", f"dimensions [0 0 0 0 0 0 0]; internalField uniform 0; boundaryField {{ inlet {{type fixedValue; value {fmt((yc < HB).astype(float))};}} wall {{type zeroGradient;}} atmosphere {{type inletOutlet; inletValue uniform 0; value uniform 0;}} frontAndBack {{type empty;}} }}\n")
    C.write(f"{case}/0/p_rgh", "volScalarField", "p_rgh", "dimensions [1 -1 -2 0 0 0 0]; internalField uniform 0; boundaryField { inlet {type fixedFluxPressure; value uniform 0;} wall {type fixedFluxPressure; value uniform 0;} atmosphere {type totalPressure; p0 uniform 0; value uniform 0;} frontAndBack {type empty;} }\n")
    C.write(f"{case}/0/U", "volVectorField", "U", f"dimensions [0 1 -1 0 0 0 0]; internalField uniform (0 0 0); boundaryField {{ inlet {{type variableHeightFlowRateInletVelocity; flowRate constant {UP * HB * dz}; alpha alpha.water; value uniform (0 0 0);}} wall {{type noSlip;}} atmosphere {{type pressureInletOutletVelocity; value uniform (0 0 0);}} frontAndBack {{type empty;}} }}\n")
    a = np.zeros((ny, nx)); u = np.zeros((ny, nx))
    for j in range(ny):
        for i in range(nx):
            h = HB if xc[i] < X0 else H0; a[j, i] = 1.0 if yc[j] < h else 0.0; u[j, i] = UP if (xc[i] < X0 and a[j, i] > 0.5) else 0.0
    return dict(dy=dy, nx=nx, ny=ny, a=a, u=u)
def analyse(case, g, times, tend):
    dy, nx, ny = g["dy"], g["nx"], g["ny"]; xs = (np.arange(nx) + 0.5) * dy; fr = []
    for t in times:
        a = J.read(f"{case}/{t:g}/alpha.water", "sca").reshape(ny, nx); u = J.read(f"{case}/{t:g}/U", "vec"); ux = u[:, 0].reshape(ny, nx); assert np.all(np.isfinite(ux))
        h = a.sum(0) * dy; q = (a * ux).sum(0) * dy; thr = 0.5 * (HB + H0); k = np.where(h > thr)[0]; xf = xs[k.max()] if len(k) else np.nan
        # sub-cell front position: linear interpolation of h through thr
        i = k.max() if len(k) else 0; xf = xs[i] + (h[i] - thr) / max(h[i] - h[min(i + 1, nx - 1)], 1e-12) * dy if len(k) else np.nan
        fr.append((t, xf, h, q, ux, a))
    t = np.array([f[0] for f in fr]); xf = np.array([f[1] for f in fr]); m = t >= 0.5; c_meas = float(np.polyfit(t[m], xf[m], 1)[0])
    t_last, xf_last, h, q, ux, a = fr[-1]; sel = lambda lo, hi: (xs >= lo) & (xs <= hi)
    h0m = float(h[sel(xf_last + 0.3, xf_last + 0.9)].mean()); hbm = float(h[sel(xf_last - 0.9, xf_last - 0.3)].mean()); qb = float(q[sel(xf_last - 0.9, xf_last - 0.3)].mean()); ub = qb / hbm; q0 = float(q[sel(xf_last + 0.3, xf_last + 0.9)].mean())
    Frm = c_meas / np.sqrt(G * h0m); hb_th = 0.5 * h0m * (np.sqrt(1 + 8 * Frm ** 2) - 1)
    # bore-frame fluxes: mass  (c-u) h conserved ; momentum  (c-u)^2 h + g h^2/2 conserved
    mass_b = (c_meas - ub) * hbm; mass_0 = c_meas * h0m; mom_b = (c_meas - ub) ** 2 * hbm + G * hbm ** 2 / 2; mom_0 = c_meas ** 2 * h0m + G * h0m ** 2 / 2
    E0 = h0m + c_meas ** 2 / (2 * G); E1 = hbm + (c_meas - ub) ** 2 / (2 * G); dE = E0 - E1; dEth = (hbm - h0m) ** 3 / (4 * hbm * h0m)
    nu = 1e-6; UX = ux; dudy, dudx = np.gradient(UX, dy, axis=0), np.gradient(UX, dy, axis=1); eps = float(np.sum(a * 1000 * nu * (2 * dudx ** 2 + dudy ** 2)) * dy * dy)
    return dict(c_meas=c_meas, c_design=CB, Fr_meas=Frm, h0_meas=h0m, h_behind_meas=hbm, h_behind_conjugate_of_meas=hb_th, h_rel_err=(hbm - hb_th) / hb_th, u_behind=ub, u_p_design=UP, mass_flux_rel_diff=(mass_b - mass_0) / mass_0, momentum_flux_rel_diff=(mom_b - mom_0) / mom_0,
                headloss=dE, headloss_theory=dEth, headloss_rel_err=(dE - dEth) / dEth, front_x_final=float(xf_last), resolved_viscous_dissipation_fraction_of_theory=float(eps / (1000 * G * mass_0 * dEth)), q_ahead=q0, depth_profile=[float(v) for v in h[::max(1, nx // 160)]])
if __name__ == "__main__":
    wd = sys.argv[1]; out = {}
    for n, co in ((4, 0.3), (6, 0.3), (8, 0.3), (12, 0.3), (8, 0.15)):
        case = f"{wd}/bore_n{n}_co{co}"; g = make(case, n, maxco=co); C.run(["blockMesh"], case, "bm.log"); shutil.copy(f"{case}/0/alpha.water.orig", f"{case}/0/alpha.water"); J.write_internal(case, g); C.run(["checkMesh"], case, "cm.log")
        try: C.run(["interFoam"], case, "run.log")
        except Exception as e: out[f"n{n}_co{co}"] = dict(failed=str(e)); print(n, "FAILED", flush=True); json.dump(out, open(os.path.join(HERE, "bore_results.json"), "w"), indent=1); continue
        out[f"n{n}_co{co}"] = analyse(case, g, [0.5, 0.75, 1.0, 1.25, 1.5], 1.5); out[f"n{n}_co{co}"]["dy"] = g["dy"]; print(n, co, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in out[f"n{n}_co{co}"].items() if k != "depth_profile"}, flush=True)
        json.dump(out, open(os.path.join(HERE, "bore_results.json"), "w"), indent=1)
