"""Matched-state upstream-thermal-history experiment (see FREEZE_TH.md).  python th.py WORKDIR KIND N [extras]  -> th_<kind>_n<N>.npz/json in ./results
KIND: step | straight. One flow solve (simpleFoam) per mesh, then unit-inlet temperature fields (scalarTransportFoam, implicit pseudo-time) for each thermal history. Adiabatic upstream wall segments remove the inlet/wall corner."""
import json, os, re, shutil, subprocess, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "cfd")); import poiseuille_graetz as P, step_vs_ramp as S
H = 0.01; NU = 1e-6; UB = 5e-3; PR = 7.0; ALPHA = NU / PR; RES = os.path.join(HERE, "results"); os.makedirs(RES, exist_ok=True)
def edges(XU): return sorted(set([-XU, -16, -13, -11, -9, -3, 0]))
def make(kind, n, case, XU=16, XD=40):
    shutil.rmtree(case, ignore_errors=True); xe = [e for e in edges(XU) if e >= -XU]; m = len(xe) - 1; V = []; idx = {}
    def add(name, x, y): idx[name] = len(V); V.append((x, y))
    for i, x in enumerate(xe): add(f"ub{i}", x, 1)
    for i, x in enumerate(xe): add(f"ut{i}", x, 2)
    if kind == "step": add("b0", 0, 0); add("b1", XD, 0); add("b2", XD, 1); add("b3", XD, 2)
    else: add("c1", XD, 1); add("c2", XD, 2)
    blk = [(f"ub{i}", f"ub{i+1}", f"ut{i+1}", f"ut{i}", int(round((xe[i + 1] - xe[i]) * n)), n) for i in range(m)]
    if kind == "step": blk += [("b0", "b1", "b2", f"ub{m}", XD * n, n), (f"ub{m}", "b2", "b3", f"ut{m}", XD * n, n)]
    else: blk += [(f"ub{m}", "c1", "c2", f"ut{m}", XD * n, n)]
    dz = H / n; nv = len(V); verts = " ".join(f"({x*H} {y*H} 0)" for x, y in V) + " " + " ".join(f"({x*H} {y*H} {dz})" for x, y in V)
    bl = " ".join(f"hex ({idx[a]} {idx[b]} {idx[c]} {idx[d]} {idx[a]+nv} {idx[b]+nv} {idx[c]+nv} {idx[d]+nv}) ({nx} {ny} 1) simpleGrading (1 1 1)" for a, b, c, d, nx, ny in blk)
    face = lambda a, b: f"({idx[a]} {idx[b]} {idx[b]+nv} {idx[a]+nv})"
    inlet = f"({idx['ub0']} {idx['ub0']+nv} {idx['ut0']+nv} {idx['ut0']})"
    outlet = (f"({idx['b1']} {idx['b2']} {idx['b2']+nv} {idx['b1']+nv}) ({idx['b2']} {idx['b3']} {idx['b3']+nv} {idx['b2']+nv})" if kind == "step" else f"({idx['c1']} {idx['c2']} {idx['c2']+nv} {idx['c1']+nv})")
    wallU = {f"wallU{i+1}": face(f"ub{i}", f"ub{i+1}") for i in range(m)}
    wallD = (face("b0", "b1") + " " + face(f"ub{m}", "b0")) if kind == "step" else face(f"ub{m}", "c1")
    lid = " ".join(face(f"ut{i}", f"ut{i+1}") for i in range(m)) + " " + (face(f"ut{m}", "b3") if kind == "step" else face(f"ut{m}", "c2"))
    fb = " ".join(f"({idx[a]} {idx[d]} {idx[c]} {idx[b]}) ({idx[a]+nv} {idx[b]+nv} {idx[c]+nv} {idx[d]+nv})" for a, b, c, d, _, _ in blk)
    bnd = f"inlet {{type patch; faces ({inlet});}} outlet {{type patch; faces ({outlet});}} " + " ".join(f"{k} {{type wall; faces ({v});}}" for k, v in wallU.items()) + f" wallD {{type wall; faces ({wallD});}} lid {{type symmetryPlane; faces ({lid});}} frontAndBack {{type empty; faces ({fb});}}"
    P.write(f"{case}/system/blockMeshDict", "dictionary", "blockMeshDict", f"convertToMeters 1; vertices ({verts}); blocks ({bl}); boundary ({bnd});\n")
    P.write(f"{case}/system/controlDict", "dictionary", "controlDict", "application simpleFoam; startFrom startTime; startTime 0; stopAt endTime; endTime 4000; deltaT 1; writeControl timeStep; writeInterval 4000; purgeWrite 1; writeFormat ascii; writePrecision 12; timeFormat general; timePrecision 8; runTimeModifiable no;\n")
    P.write(f"{case}/system/fvSchemes", "dictionary", "fvSchemes", "ddtSchemes {default steadyState;} gradSchemes {default Gauss linear;} divSchemes {default none; div(phi,U) bounded Gauss linearUpwind grad(U); div(phi,T) bounded Gauss linearUpwind grad(T); div((nuEff*dev2(T(grad(U))))) Gauss linear;} laplacianSchemes {default Gauss linear corrected;} interpolationSchemes {default linear;} snGradSchemes {default corrected;}\n")
    P.write(f"{case}/system/fvSolution", "dictionary", "fvSolution", "solvers { p {solver GAMG; smoother GaussSeidel; tolerance 1e-10; relTol 0.01;} U {solver smoothSolver; smoother symGaussSeidel; tolerance 1e-10; relTol 0.01;} T {solver smoothSolver; smoother symGaussSeidel; tolerance 1e-14; relTol 0; maxIter 5000;} } SIMPLE {nNonOrthogonalCorrectors 0; consistent yes; residualControl {p 1e-8; U 1e-9;}} relaxationFactors {equations {U 0.9; T 0.95;} fields {p 0.3;}}\n")
    P.write(f"{case}/constant/transportProperties", "dictionary", "transportProperties", f"transportModel Newtonian; nu [0 2 -1 0 0 0 0] {NU}; DT [0 2 -1 0 0 0 0] {ALPHA};\n")
    P.write(f"{case}/constant/turbulenceProperties", "dictionary", "turbulenceProperties", "simulationType laminar;\n")
    walls = " ".join(list(wallU) + ["wallD"]); P.write(f"{case}/0/U", "volVectorField", "U", f"dimensions [0 1 -1 0 0 0 0]; internalField uniform ({UB} 0 0); boundaryField {{ inlet {{type fixedValue; value uniform ({UB} 0 0);}} outlet {{type zeroGradient;}} " + " ".join(f"{w} {{type noSlip;}}" for w in walls.split()) + " lid {type symmetryPlane;} frontAndBack {type empty;} }\n")
    P.write(f"{case}/0/p", "volScalarField", "p", "dimensions [0 2 -2 0 0 0 0]; internalField uniform 0; boundaryField { inlet {type zeroGradient;} outlet {type fixedValue; value uniform 0;} " + " ".join(f"{w} {{type zeroGradient;}}" for w in walls.split()) + " lid {type symmetryPlane;} frontAndBack {type empty;} }\n")
    return list(wallU), xe
def flow(case):
    P.run(["blockMesh"], case, "bm.log"); P.run(["checkMesh"], case, "cm.log"); P.run(["simpleFoam"], case, "flow.log")
    return sorted([d for d in os.listdir(case) if d.replace(".", "").isdigit() and float(d) > 0], key=float)[-1]
def thermal(case, tf, wallU, xe, cold_start, dt, tend, out):
    shutil.rmtree(out, ignore_errors=True); os.makedirs(out); shutil.copytree(f"{case}/constant", f"{out}/constant"); shutil.copytree(f"{case}/system", f"{out}/system"); os.makedirs(f"{out}/{tf}")
    for f in ("U", "phi"): shutil.copy(f"{case}/{tf}/{f}", f"{out}/{tf}/{f}")
    cold = [w for i, w in enumerate(wallU) if xe[i] >= cold_start - 1e-9]
    bc = lambda w: ("{type fixedValue; value uniform 0;}" if (w in cold or w == "wallD") else "{type zeroGradient;}")
    P.write(f"{out}/{tf}/T", "volScalarField", "T", "dimensions [0 0 0 1 0 0 0]; internalField uniform 1; boundaryField { inlet {type fixedValue; value uniform 1;} outlet {type zeroGradient;} " + " ".join(f"{w} {bc(w)}" for w in wallU + ["wallD"]) + " lid {type symmetryPlane;} frontAndBack {type empty;} }\n")
    t0 = float(tf); txt = open(f"{out}/system/controlDict").read().replace("application simpleFoam;", "application scalarTransportFoam;").replace("startTime 0;", f"startTime {t0};").replace("deltaT 1;", f"deltaT {dt};").replace("endTime 4000;", f"endTime {t0 + tend};").replace("writeControl timeStep;", "writeControl runTime;").replace("writeInterval 4000;", f"writeInterval {tend/10};").replace("purgeWrite 1;", "purgeWrite 3;")
    open(f"{out}/system/controlDict", "w").write(txt); P.run(["scalarTransportFoam"], out, "T.log")
    return sorted([d for d in os.listdir(out) if d.replace(".", "").isdigit() and float(d) > t0], key=float)
class Mesh:
    def __init__(self, case, n):
        self.n = n; self.pts, self.faces, self.owner, self.patches = S.mesh_data(case); nn, b = S.read_list(f"{case}/constant/polyMesh/neighbour"); self.nbr = np.array([int(l) for l in b.strip().splitlines()[:nn]])
        ncell = int(np.max(self.owner)) + 1; self.C, self.fc = S.cell_centres(self.pts, self.faces, self.owner, self.nbr, ncell); self.ncell = ncell; self.dz = H / n
def read_T(out, tdir): return S.read_cell_field(f"{out}/{tdir}/T", "sca")
def observables(M, T, U, xc_list=(-1.0,)):
    D = ALPHA; s, c = M.patches["wallD"]; x = []; q = []; ln = []
    for k in range(s, s + c):
        p = M.pts[M.faces[k]]; z0 = [pp[:2] for pp in p if abs(pp[2]) < 1e-12]; e = z0[1] - z0[0]; L = np.linalg.norm(e); nrm = np.array([e[1], -e[0]]) / L; P_ = M.owner[k]; d = M.fc[k][:2] - M.C[P_][:2]
        if nrm @ d < 0: nrm = -nrm
        delta = 1.0 / max(nrm @ d, 0.05 * np.linalg.norm(d)); x.append(M.fc[k][0] / H); q.append(D * T[P_] * delta); ln.append(L / H)
    x = np.array(x); q = np.array(q); ln = np.array(ln); o = np.argsort(x, kind="stable"); res = dict(x=x[o], q=q[o], ln=ln[o], Tb_sec={}, prof={}, qwall_sec={})
    ux = U[:, 0]; Cx = M.C[:, 0] / H; Cy = M.C[:, 1] / H
    for xc in xc_list:
        dx = 1.0 / M.n; col = np.where((np.abs(Cx - xc) < 0.5 * dx + 1e-9) & (Cy > 1.0))[0]; col = col[np.argsort(Cy[col])]; w = ux[col]; Tb = float(np.sum(w * T[col]) / np.sum(w))
        res["Tb_sec"][xc] = Tb; res["prof"][xc] = (Cy[col], T[col], w)
    # wall flux on the upstream cold wall at the section column: first cell above the bottom (y near 1)
    for xc in xc_list:
        col = np.where((np.abs(Cx - xc) < 0.5 / M.n + 1e-9) & (Cy > 1.0))[0]; k = col[np.argmin(Cy[col])]; res["qwall_sec"][xc] = float(D * T[k] / (0.5 * (1.0 / M.n) * H))
    # flow-weighted bulk temperature along the channel (stations every 0.5 H, x = -1 ... 30): used for the conductance-shape test (a(x) = q / T_b)
    st = np.arange(-1.0, 30.01, 0.5); dx = 1.0 / M.n; tbx = []
    for xs_ in st:
        col = np.where(np.abs(Cx - xs_) < 0.5 * dx + 1e-9)[0]; w = ux[col]; tbx.append(float(np.sum(w * T[col]) / np.sum(w)))
    res["Tb_x"] = (st.tolist(), tbx)
    # outlet flow-weighted temperature
    s2, c2 = M.patches["outlet"]; num = den = 0.0
    for k in range(s2, s2 + c2):
        p = M.pts[M.faces[k]]; z0 = [pp[:2] for pp in p if abs(pp[2]) < 1e-12]; ar = np.linalg.norm(z0[1] - z0[0]) * M.dz; u = U[M.owner[k], 0]; num += u * ar * T[M.owner[k]]; den += u * ar
    res["Tb_out"] = num / den; return res
def Q_of(res, X): m = (res["x"] >= -1e-9) & (res["x"] <= X); return float(np.sum(res["q"][m] * res["ln"][m]))
def run_pipeline(wd, kind, n, tag="base", XU=16, XD=40, histories=(("A", -3.0), ("B", -11.0)), dts=(10.0,), tend=1500.0):
    case = f"{wd}/{kind}_n{n}_{tag}"; wallU, xe = make(kind, n, case, XU, XD); t0 = time.time(); tf = flow(case); tflow = time.time() - t0; M = Mesh(case, n); out = {}
    U = S.read_cell_field(f"{case}/{tf}/U", "vec")
    for hname, cs in histories:
        for dt in dts:
            v = f"{case}_{hname}_dt{dt:g}"; t1 = time.time(); tdirs = thermal(case, tf, wallU, xe, cs, dt, tend, v); wall = time.time() - t1
            T = read_T(v, tdirs[-1]); Tp = read_T(v, tdirs[-2]); res = observables(M, T, U, (-0.5, -1.0, -2.0)); resp = observables(M, Tp, U, (-1.0,))
            out[f"{hname}_dt{dt:g}"] = dict(res=res, Q10=Q_of(res, 10), Q30=Q_of(res, 30), Q10_prev=Q_of(resp, 10), Q30_prev=Q_of(resp, 30), wall_s=wall, cells=M.ncell, tflow=tflow, tdirs=tdirs[-2:])
            shutil.rmtree(v, ignore_errors=True)
    return out
def save(out, name):
    ser = {}
    for k, v in out.items():
        r = v["res"]; ser[k] = dict(x=r["x"].tolist(), q=r["q"].tolist(), ln=r["ln"].tolist(), Tb_sec={str(a): b for a, b in r["Tb_sec"].items()}, qwall_sec={str(a): b for a, b in r["qwall_sec"].items()}, Tb_out=r["Tb_out"], Tb_x=r["Tb_x"], prof={str(a): [list(map(float, z)) for z in b] for a, b in r["prof"].items()},
                      Q10=v["Q10"], Q30=v["Q30"], Q10_prev=v["Q10_prev"], Q30_prev=v["Q30_prev"], wall_s=v["wall_s"], cells=v["cells"], tflow=v["tflow"])
    json.dump(ser, open(os.path.join(RES, name + ".json"), "w"))
if __name__ == "__main__":
    wd, kind, n = sys.argv[1], sys.argv[2], int(sys.argv[3]); extra = sys.argv[4] if len(sys.argv) > 4 else "base"
    if extra == "base": save(run_pipeline(wd, kind, n), f"{kind}_n{n}_base")
    elif extra == "n16full":
        save(run_pipeline(wd, kind, n, dts=(5.0, 10.0, 20.0), histories=(("A", -3.0), ("B", -11.0), ("B9", -9.0), ("B13", -13.0))), f"{kind}_n{n}_full")
        save(run_pipeline(wd, kind, n, tag="xd28", XD=28), f"{kind}_n{n}_xd28"); save(run_pipeline(wd, kind, n, tag="xu24", XU=24), f"{kind}_n{n}_xu24")
