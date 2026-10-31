"""Free-surface thermal verification of interThermalFoam (v2 custom solver: interFoam + sensible-heat equation in BOTH phases; passive temperature, no buoyancy, no phase change).
Thermal model: d(Cv T)/dt + div(T [Cw alphaPhi + Ca (phi - alphaPhi)]) - div(k_alpha grad T) = 0 over the whole domain, Cv = alpha Cw + (1-alpha) Ca [J/m3/K], k_alpha = alpha k_w + (1-alpha) k_a (arithmetic mean), alphaPhi = MULES volume flux of water.
Free-surface thermal condition: NONE imposed explicitly: temperature and the conductive flux are continuous through the (diffuse, ~1-2 cell) interface because one equation is solved with alpha-weighted properties; air has effusivity ratio sqrt(k C)_a/sqrt(k C)_w = 3.3e-3 so the surface is nearly adiabatic for water.
No evaporation, radiation or convective air-side coefficient: the air-side boundary is whatever the air column does (declared limitation).
Cases:  T1  vertical diffusion eigenmode in water at rest (adiabatic bottom, water layer h = 0.05 m, air above): T = T_m + A cos(pi y/h) exp(-alpha (pi/h)^2 t), alpha = k_w/Cw = 1.33e-7 m2/s; air initialised at the interface value. Three resolutions + time-step check. Energy ledger: closed domain (all walls adiabatic).
        T2  pure advection of a Gaussian temperature pulse by a uniform flow U = 0.1 m/s through a cyclic two-layer (water/air) domain, one transit; exact solution = translated Gaussian (plus negligible diffusion). Three resolutions + time-step check. Local (L2 profile error) and global (energy) checks."""
import json, os, shutil, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); import common as C
C.ENV["PATH"] = "/tmp/claude-0/ofbuild/bin:" + C.ENV.get("PATH", ""); C.ENV["LD_LIBRARY_PATH"] = "/usr/lib:/usr/lib/x86_64-linux-gnu/openmpi/lib:/tmp/claude-0/ofbuild/lib:" + C.ENV.get("LD_LIBRARY_PATH", "")
HW = 0.05; HA = 0.05; ALPHA_W = 0.561 / 4.217e6
def thermal_props(case):
    C.write(f"{case}/constant/thermalProperties", "dictionary", "thermalProperties", "water {rho 1000; cp 4217; k 0.561;} air {rho 1; cp 1005; k 0.0257;}\n")
    s = open(f"{case}/system/fvSchemes").read().replace("laplacianSchemes", "div(CvPhi,T) Gauss linearUpwind grad(T); laplacianSchemes") if False else None
    C.write(f"{case}/system/fvSchemes", "dictionary", "fvSchemes", "ddtSchemes {default Euler;} gradSchemes {default Gauss linear;} divSchemes {div(rhoPhi,U) Gauss linearUpwind grad(U); div(phi,alpha) Gauss vanLeer; div(phirb,alpha) Gauss linear; div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear; div(CvPhi,T) Gauss linearUpwind grad(T);} laplacianSchemes {default Gauss linear corrected;} interpolationSchemes {default linear;} snGradSchemes {default corrected;}\n")
    t = open(f"{case}/system/fvSolution").read().replace("solvers {", 'solvers { T {solver smoothSolver; smoother symGaussSeidel; tolerance 1e-14; relTol 0;} ', 1); open(f"{case}/system/fvSolution", "w").write(t)
def blocks(case, nx, ny_w, ny_a, Lx, cyclic, wall_sides):
    dz = Lx / nx; yl = [0, HW, HW + HA]; V = [(x, y) for y in yl for x in (0, Lx)]       # index = 2*iy + ix
    verts = " ".join(f"({x} {y} 0)" for x, y in V) + " " + " ".join(f"({x} {y} {dz})" for x, y in V); nv = len(V)
    hx = lambda i0, i1, i2, i3, nxx, nyy: f"hex ({i0} {i1} {i3} {i2} {i0+nv} {i1+nv} {i3+nv} {i2+nv}) ({nxx} {nyy} 1) simpleGrading (1 1 1)"
    # block w: bottom row (0,1) top row (2,3); block a: bottom row (2,3) top row (4,5)
    bl = hx(0, 1, 2, 3, nx, ny_w) + " " + hx(2, 3, 4, 5, nx, ny_a)
    bottom = "(0 %d %d 1)" % (nv, nv + 1)   # (0 1 1' 0') ordering: (v0 v1 v1' v0')
    bottom = f"(0 1 {1+nv} {nv})"; top = f"(4 {4+nv} {5+nv} 5)"
    leftw = f"(0 {nv} {2+nv} 2)"; lefta = f"(2 {2+nv} {4+nv} 4)"; rightw = f"(1 3 {3+nv} {1+nv})"; righta = f"(3 5 {5+nv} {3+nv})"
    fb = f"(0 2 3 1) ({nv} {1+nv} {3+nv} {2+nv}) (2 4 5 3) ({2+nv} {3+nv} {5+nv} {4+nv})"
    side = (f"left {{type cyclic; neighbourPatch right; faces ({leftw} {lefta});}} right {{type cyclic; neighbourPatch left; faces ({rightw} {righta});}}" if cyclic else f"walls {{type wall; faces ({leftw} {lefta} {rightw} {righta});}}")
    C.write(f"{case}/system/blockMeshDict", "dictionary", "blockMeshDict", f"convertToMeters 1; vertices ({verts}); blocks ({bl}); boundary (bottom {{type wall; faces ({bottom});}} atmosphere {{type patch; faces ({top});}} {side} frontAndBack {{type empty; faces ({fb});}});\n")
    return dz
def fields(case, cyclic, Tfun, Ufun, nx, ny_w, ny_a, Lx):
    sides = "left {type cyclic;} right {type cyclic;}" if cyclic else "walls {type wall;}"
    wall_alpha = "" if cyclic else "walls {type zeroGradient;}"
    dx = Lx / nx; xs = (np.arange(nx) + 0.5) * dx; yw = (np.arange(ny_w) + 0.5) * HW / ny_w; ya = HW + (np.arange(ny_a) + 0.5) * HA / ny_a
    T = np.concatenate([np.array([[Tfun(x, y) for x in xs] for y in yw]).ravel(), np.array([[Tfun(x, y) for x in xs] for y in ya]).ravel()])
    al = np.concatenate([np.ones(nx * ny_w), np.zeros(nx * ny_a)])
    u = Ufun
    sidebc = lambda bc: (f"left {{type cyclic;}} right {{type cyclic;}}" if cyclic else f"walls {{{bc}}}")
    C.write(f"{case}/0/alpha.water", "volScalarField", "alpha.water", f"dimensions [0 0 0 0 0 0 0]; internalField nonuniform List<scalar> {al.size}({' '.join('%g' % x for x in al)}); boundaryField {{ bottom {{type zeroGradient;}} atmosphere {{type inletOutlet; inletValue uniform 0; value uniform 0;}} {sidebc('type zeroGradient;')} frontAndBack {{type empty;}} }}\n")
    C.write(f"{case}/0/p_rgh", "volScalarField", "p_rgh", f"dimensions [1 -1 -2 0 0 0 0]; internalField uniform 0; boundaryField {{ bottom {{type fixedFluxPressure; value uniform 0;}} atmosphere {{type totalPressure; p0 uniform 0; value uniform 0;}} {sidebc('type fixedFluxPressure; value uniform 0;')} frontAndBack {{type empty;}} }}\n")
    C.write(f"{case}/0/U", "volVectorField", "U", f"dimensions [0 1 -1 0 0 0 0]; internalField uniform ({u} 0 0); boundaryField {{ bottom {{type {'slip' if u else 'noSlip'};}} atmosphere {{type pressureInletOutletVelocity; value uniform ({u} 0 0);}} {sidebc('type slip;' if u else 'type noSlip;')} frontAndBack {{type empty;}} }}\n")
    C.write(f"{case}/0/T", "volScalarField", "T", f"dimensions [0 0 0 1 0 0 0]; internalField nonuniform List<scalar> {T.size}({' '.join('%.12g' % x for x in T)}); boundaryField {{ bottom {{type zeroGradient;}} atmosphere {{type zeroGradient;}} {sidebc('type zeroGradient;')} frontAndBack {{type empty;}} }}\n")
    return T
def run_case(kind, n, cof, wd, tend):
    case = f"{wd}/{kind}_n{n}_co{cof}"; shutil.rmtree(case, ignore_errors=True)
    if kind == "T1":
        nx = 4; ny_w = n; ny_a = max(2, n // 2); Lx = 0.02; cyc = False; U = 0.0; A = 1.0; Tm = 10.0
        Tf = lambda x, y: (Tm + A * np.cos(np.pi * y / HW)) if y < HW else Tm - A
    else:
        nx = n; ny_w = 4; ny_a = 4; Lx = 0.5; cyc = True; U = 0.1; Tf = lambda x, y: 1.0 + np.exp(-(((x - 0.15) / 0.03) ** 2))
    # capillary time-step limit (Brackbill): dt < sqrt(rho_avg dx^3/(2 pi sigma)); half of it is used
    blocks(case, nx, ny_w, ny_a, Lx, cyc, not cyc); C.physics(case)
    if kind == 'T1': t = open(f'{case}/system/fvSolution').read().replace('PIMPLE {', 'PIMPLE {frozenFlow yes; '); open(f'{case}/system/fvSolution', 'w').write(t)      # water kept at rest: spurious interface currents (3e-4 m/s) would otherwise dominate the diffusion test (Pe ~ 100)
    dmin = min(Lx / nx, HW / ny_w); capdt = 0.5 * np.sqrt(500.0 * dmin ** 3 / (2 * np.pi * 0.07)); C.control(case, tend, min(1e-4, capdt), cof, tend / 2, maxdt=capdt); thermal_props(case); T0 = fields(case, cyc, Tf, U, nx, ny_w, ny_a, Lx)
    C.run(["blockMesh"], case, "bm.log"); C.run(["checkMesh"], case, "cm.log"); C.run(["interThermalFoam"], case, "run.log")
    return case, dict(nx=nx, ny_w=ny_w, ny_a=ny_a, Lx=Lx, U=U)
def read(path, kind="sca"):
    t = open(path).read(); i = t.index("internalField"); j = t.index("(", i); n = int(t[i:j].split()[-1]); body = t[j + 1:].split("\n)")[0].strip().splitlines(); return np.array([float(l) for l in body[:n]])
def ledger(case):
    E = [float(l.split("=")[1].split()[0]) for l in open(f"{case}/run.log") if l.startswith("ThermalEnergy")]; return dict(E_first=E[0], E_last=E[-1], rel_drift=(E[-1] - E[0]) / E[0], abs_drift_J_per_m=E[-1] - E[0])
def analyse(kind, case, g, tend):
    tl = f"{tend:g}"; T = read(f"{case}/{tl}/T"); a = read(f"{case}/{tl}/alpha.water"); nx, ny_w, ny_a, Lx = g["nx"], g["ny_w"], g["ny_a"], g["Lx"]; out = dict(ledger=ledger(case), Tmin=float(T.min()), Tmax=float(T.max()))
    if kind == "T1":
        yw = (np.arange(ny_w) + 0.5) * HW / ny_w; Tw = T[: nx * ny_w].reshape(ny_w, nx).mean(1); c = np.cos(np.pi * yw / HW); Am = float(np.sum((Tw - 10.0) * c) / np.sum(c * c)); ex = np.exp(-ALPHA_W * (np.pi / HW) ** 2 * tend)
        # exact discrete water mean is 10; amplitude of the cosine mode; also profile error
        prof_err = float(np.sqrt(np.mean((Tw - (10.0 + ex * c)) ** 2)) / 1.0); out.update(amp=Am, amp_exact=float(ex), amp_rel_err=(Am - ex) / ex, profile_rms_err=prof_err, T_air_mean=float(T[nx * ny_w:].mean()), mean_water_T=float(Tw.mean()))
    else:
        dx = Lx / nx; xs = (np.arange(nx) + 0.5) * dx; Tw = T[: nx * ny_w].reshape(ny_w, nx).mean(0); Ta = T[nx * ny_w:].reshape(ny_a, nx).mean(0); sh = g["U"] * tend
        ex = 1.0 + np.exp(-(((((xs - 0.15 - sh) + 0.25) % 0.5 - 0.25) / 0.03) ** 2)); out.update(L2_err_water=float(np.sqrt(np.mean((Tw - ex) ** 2))), L2_err_air=float(np.sqrt(np.mean((Ta - ex) ** 2))), peak=float(Tw.max() - 1), peak_exact=1.0, dx=dx, y_variation=float(np.abs(T[: nx * ny_w].reshape(ny_w, nx) - Tw[None, :]).max()))
    return out
if __name__ == "__main__":
    wd = sys.argv[1]; res = {}
    for kind, tend, ns in (("T1", 300.0, (10, 20, 40)), ("T2", 2.5, (50, 100, 200))):
        for n in ns:
            for cof in ((0.2,) if not (n == ns[1]) else (0.2, 0.1)):
                try: case, g = run_case(kind, n, cof, wd, tend); r = analyse(kind, case, g, tend)
                except Exception as e: r = dict(failed=str(e))
                res[f"{kind}_n{n}_co{cof}"] = r; print(kind, n, cof, json.dumps(r)[:400], flush=True); json.dump(res, open(os.path.join(HERE, "thermal_fs_results.json"), "w"), indent=1)
