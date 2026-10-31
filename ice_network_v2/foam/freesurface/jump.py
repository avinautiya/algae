"""B. Rectangular-channel hydraulic jump, 2-D, stock interFoam (laminar, nu_w = 1e-6, sigma = 0.07), frictionless-ish noSlip bed.
Case: h1 = 0.05 m, Fr1 = U1/sqrt(g h1) = 2.5 (U1 = 1.751 m/s, Re_h1 = 8.8e4, We_h1 = rho U1^2 h1/sigma = 2.2e3, Bo = rho g h1^2/sigma = 350). Conjugate depth (Belanger): h2/h1 = (sqrt(1+8 Fr1^2) - 1)/2 = 3.0711, h2 = 0.1536 m.
BCs: inlet U by variableHeightFlowRateInletVelocity (flow rate U1 h1 dz) with alpha prescribed per face (water for y < h1), p_rgh fixedFluxPressure; outlet p_rgh = hydrostatic profile of water level h2 (per-face fixedValue: (rho_w-rho_a) g h2 in water, 0 in air consistent with the atmosphere BC; a first version with an inconsistent air value drove 114 m/s air flow and was abandoned), U zeroGradient, alpha inletOutlet; top atmosphere totalPressure; bed noSlip.
Initial: depth h1 and U1 for x < 0.8 m, depth h2 and U1 h1/h2 for x > 0.8 m. L = 2 m, domain height 0.3 m. Three resolutions dy = h1/4, h1/6, h1/8 (dx = dy); t_end = 8 s; time-average of 0.5 s snapshots over t in [6, 8] s.
APPLICABILITY: Belanger's relation holds for a (quasi-)frictionless horizontal rectangular channel with hydrostatic pressure upstream and downstream and uniform velocity; here the 2-D laminar run at Re ~ 1e5 has no turbulence model, so the dissipation in the jump is carried by resolved eddies and numerical diffusion; bed friction (noSlip, laminar) is small."""
import json, os, shutil, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); import common as C
G = 9.81; H1 = 0.05; FR = 2.5; U1 = FR * np.sqrt(G * H1); H2 = 0.5 * H1 * (np.sqrt(1 + 8 * FR ** 2) - 1); U2 = U1 * H1 / H2; L = 2.0; HT = 0.3; XJ = 0.8
def make(case, n, tend=8.0, maxco=0.3, dt_ratio=1.0):
    shutil.rmtree(case, ignore_errors=True); dy = H1 / n; nx = int(round(L / dy)); ny = int(round(HT / dy)); dz = dy
    C.write(f"{case}/system/blockMeshDict", "dictionary", "blockMeshDict", f"convertToMeters 1; vertices ((0 0 0) ({L} 0 0) ({L} {HT} 0) (0 {HT} 0) (0 0 {dz}) ({L} 0 {dz}) ({L} {HT} {dz}) (0 {HT} {dz})); blocks (hex (0 1 2 3 4 5 6 7) ({nx} {ny} 1) simpleGrading (1 1 1)); boundary (inlet {{type patch; faces ((0 4 7 3));}} outlet {{type patch; faces ((1 2 6 5));}} bed {{type wall; faces ((0 1 5 4));}} atmosphere {{type patch; faces ((3 7 6 2));}} frontAndBack {{type empty; faces ((0 3 2 1) (4 5 6 7));}});\n")
    C.physics(case); C.control(case, tend, 1e-4, maxco, 0.5, maxdt=0.005 * dt_ratio)
    yc = (np.arange(ny) + 0.5) * dy; xc = (np.arange(nx) + 0.5) * dy
    a_in = (yc < H1).astype(float); ph_out = np.where(yc < H2, 999.0 * G * H2, 0.0)      # p_rgh = p - rho gh: air column at p_rgh = 0 (matches the atmosphere totalPressure p0 = 0), water (rho_w - rho_a) g h2
    fmt = lambda v: "nonuniform List<scalar> %d(%s)" % (len(v), " ".join("%.10g" % x for x in v))
    C.write(f"{case}/0/alpha.water.orig", "volScalarField", "alpha.water", f"dimensions [0 0 0 0 0 0 0]; internalField uniform 0; boundaryField {{ inlet {{type fixedValue; value {fmt(a_in)};}} outlet {{type inletOutlet; inletValue uniform 0; value uniform 0;}} bed {{type zeroGradient;}} atmosphere {{type inletOutlet; inletValue uniform 0; value uniform 0;}} frontAndBack {{type empty;}} }}\n")
    C.write(f"{case}/0/p_rgh", "volScalarField", "p_rgh", f"dimensions [1 -1 -2 0 0 0 0]; internalField uniform 0; boundaryField {{ inlet {{type fixedFluxPressure; value uniform 0;}} outlet {{type fixedValue; value {fmt(ph_out)};}} bed {{type fixedFluxPressure; value uniform 0;}} atmosphere {{type totalPressure; p0 uniform 0; value uniform 0;}} frontAndBack {{type empty;}} }}\n")
    C.write(f"{case}/0/U", "volVectorField", "U", f"dimensions [0 1 -1 0 0 0 0]; internalField uniform (0 0 0); boundaryField {{ inlet {{type variableHeightFlowRateInletVelocity; flowRate constant {U1 * H1 * dz}; alpha alpha.water; value uniform (0 0 0);}} outlet {{type zeroGradient;}} bed {{type noSlip;}} atmosphere {{type pressureInletOutletVelocity; value uniform (0 0 0);}} frontAndBack {{type empty;}} }}\n")
    # initial fields written directly (setFields-free): alpha, U for the whole mesh, in x-fastest cell order
    a = np.zeros((ny, nx)); u = np.zeros((ny, nx))
    for j in range(ny):
        for i in range(nx):
            h = H1 if xc[i] < XJ else H2
            a[j, i] = 1.0 if yc[j] < h else 0.0
            u[j, i] = (U1 if xc[i] < XJ else U2) if a[j, i] > 0.5 else 0.0
    return dict(dy=dy, nx=nx, ny=ny, a=a, u=u)
def write_internal(case, g):
    def put(name, body_internal):
        p = f"{case}/0/{name}"; t = open(p).read(); i = t.index("internalField"); j = t.index(";", i); open(p, "w").write(t[:i] + body_internal + t[j:])
    put("alpha.water", "internalField nonuniform List<scalar> %d(%s)" % (g["a"].size, " ".join("%g" % x for x in g["a"].ravel())))
    put("U", "internalField nonuniform List<vector> %d(%s)" % (g["u"].size, " ".join("(%g 0 0)" % x for x in g["u"].ravel())))
def read(path, kind):
    t = open(path).read(); i = t.index("internalField"); j = t.index("(", i); n = int(t[i:j].split()[-1]); body = t[j + 1:].split("\n)")[0].strip().splitlines()
    return np.array([[float(v) for v in l.strip("()").split()] for l in body[:n]]) if kind == "vec" else np.array([float(l) for l in body[:n]])
def analyse(case, g, times):
    dy, nx, ny = g["dy"], g["nx"], g["ny"]; res = dict(times=times); Hs = []; Qs = []; Ms = []; hu = []; UU = []; AA = []
    for t in times:
        a = read(f"{case}/{t:g}/alpha.water", "sca").reshape(ny, nx); u = read(f"{case}/{t:g}/U", "vec"); ux = u[:, 0].reshape(ny, nx); uy = u[:, 1].reshape(ny, nx)
        assert np.all(np.isfinite(a)) and np.all(np.isfinite(ux))
        AA.append(a); UU.append((ux, uy))
    A = np.mean(AA, axis=0); UX = np.mean([x[0] for x in UU], axis=0); UY = np.mean([x[1] for x in UU], axis=0)
    h = A.sum(0) * dy; Q = (A * UX).sum(0) * dy; yc = (np.arange(ny) + 0.5) * dy
    # momentum flux per unit width (hydrostatic pressure, using the time-mean depth): M = Q^2/h + g h^2/2 (section-integrated alpha*u^2 is used instead of Q^2/h to include profile shape)
    Mflux = (A * UX ** 2).sum(0) * dy + G * h ** 2 / 2
    def sec(x0, x1): m = (np.arange(nx) + 0.5) * dy; k = (m >= x0) & (m <= x1); return k
    up = sec(0.15, 0.45); dn = sec(1.5, 1.9); h1m = float(h[up].mean()); h2m = float(h[dn].mean()); Q1 = float(Q[up].mean()); Q2 = float(Q[dn].mean()); U1m = Q1 / h1m; U2m = Q2 / h2m; Fr1m = U1m / np.sqrt(G * h1m)
    h2_th = 0.5 * h1m * (np.sqrt(1 + 8 * Fr1m ** 2) - 1)
    E1 = h1m + Q1 ** 2 / (2 * G * h1m ** 2); E2 = h2m + Q2 ** 2 / (2 * G * h2m ** 2); dE = E1 - E2; dE_th = (h2m - h1m) ** 3 / (4 * h1m * h2m)
    M1 = float(Mflux[up].mean()); M2 = float(Mflux[dn].mean())
    # toe: first x where mean depth exceeds h1 + 0.5 (h2 - h1)
    xs = (np.arange(nx) + 0.5) * dy; thr = h1m + 0.5 * (h2m - h1m); toe = float(xs[np.argmax(h > thr)])
    # resolved viscous dissipation (water only, time-mean fields -> underestimates fluctuating dissipation): rho nu 2 S:S integrated
    dudy, dudx = np.gradient(UX, dy, axis=0), np.gradient(UX, dy, axis=1); dvdy, dvdx = np.gradient(UY, dy, axis=0), np.gradient(UY, dy, axis=1); S2 = 2 * (dudx ** 2 + dvdy ** 2) + (dudy + dvdx) ** 2
    eps = float(np.sum(A * 1000 * 1e-6 * S2) * dy * dy); P_theory = 1000 * G * Q1 * dE_th      # W per unit width
    res.update(h1=h1m, h2=h2m, Fr1=Fr1m, h2_conjugate=h2_th, h2_rel_err=(h2m - h2_th) / h2_th, Q_in=Q1, Q_out=Q2, Q_rel_diff=(Q2 - Q1) / Q1, M_up=M1, M_down=M2, M_rel_diff=(M2 - M1) / M1, headloss=dE, headloss_theory=dE_th, headloss_rel_err=(dE - dE_th) / dE_th,
               toe_x=toe, resolved_viscous_dissipation_W_per_m=eps, theory_loss_W_per_m=P_theory, resolved_fraction=eps / P_theory, depth_profile=[float(x) for x in h[::max(1, nx // 100)]])
    return res
if __name__ == "__main__":
    wd = sys.argv[1]; out = {}; plan = [(4, 1.0), (6, 1.0), (8, 1.0)]
    for n, r in plan:
        case = f"{wd}/jump_n{n}"; g = make(case, n); C.run(["blockMesh"], case, "bm.log"); shutil.copy(f"{case}/0/alpha.water.orig", f"{case}/0/alpha.water"); write_internal(case, g); C.run(["checkMesh"], case, "cm.log")
        try: C.run(["interFoam"], case, "run.log")
        except Exception as e: out[f"n{n}"] = dict(failed=str(e)); print(n, "FAILED", e, flush=True); json.dump(out, open(os.path.join(HERE, "jump_results.json"), "w"), indent=1); continue
        out[f"n{n}"] = analyse(case, g, [6.0, 6.5, 7.0, 7.5, 8.0]); out[f"n{n}"]["dy"] = g["dy"]; print(n, {k: v for k, v in out[f"n{n}"].items() if k != "depth_profile"}, flush=True)
        json.dump(out, open(os.path.join(HERE, "jump_results.json"), "w"), indent=1)
