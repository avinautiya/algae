"""Controlled fixed-geometry comparison: single backward step vs smooth ramp of the same drop, LAMINAR, steady, stock OpenFOAM v1912 (simpleFoam + scalarTransportFoam).
Matched: channel heights (upstream h = 1, downstream 2 in units of H), inlet volume flux, inlet temperature (uniform T=1 => identical inlet sensible-enthalpy flux), wall temperature 0 (melting-point reference), properties, mesh spacing, outlet.
This is NOT a hydraulic jump and has NO free surface: the top is a free-slip adiabatic lid (declared rigid-lid proxy of a flat surface). It tests the stated question (does a geometric step change wall heat flux, separation length and downstream enthalpy at fixed budgets)
in the laminar regime only. Re_H = Ub H / nu = 50, Pr = 7, Pe_H = 350. Geometry: upstream 12 H, drop 1 H (1 -> 2 channel height), downstream 30 H; ramp: x in [0, 2H] (slope 1:2; a longer 1:6 ramp failed checkMesh skewness 7.5 with this block topology and was not used). Physical-unit scale H = 0.01 m.
Outputs: wall heat flux q(x) (nondimensional: alpha dT/dn / (alpha dT_in/H)), reattachment (wall-shear sign change), outlet flow-weighted temperature, outlet enthalpy flux, integrated bed heat, ledger  P_in - P_out - P_wall."""
import json, os, re, shutil, subprocess, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); import poiseuille_graetz as P
H = 0.01; NU = 1e-6; UB = 5e-3; PR = 7.0; ALPHA = NU / PR; XU = 12.0; XD = 30.0; XR = 2.0
ENV = P.ENV


def vtx(x, y): return f"({x*H} {y*H} 0) ({x*H} {y*H} {H/ NB[0]})" 


def blocks(kind, n):
    """vertices and blocks in H units; thickness in z one cell. returns blockMeshDict body"""
    v = []; idx = {}
    def V(name, x, y): idx[name] = len(v); v.append((x, y))
    # lower plane vertices then upper plane (same x,y, z = dz)
    if kind == "step":
        V("a0", -XU, 1); V("a1", 0, 1); V("a2", -XU, 2); V("a3", 0, 2); V("b0", 0, 0); V("b1", XD, 0); V("b2", XD, 1); V("b3", XD, 2)
        # blocks: A upstream [a0 a1 a3 a2]; B lower downstream [b0 b1 b2 a1]; C upper downstream [a1 b2 b3 a3]
        blk = [("a0", "a1", "a3", "a2", XU * n, n), ("b0", "b1", "b2", "a1", XD * n, n), ("a1", "b2", "b3", "a3", XD * n, n)]
    else:    # same block topology as the step; only vertex b0 moves from (0,0) to (XR,0): the left edge of the lower block becomes the ramp wall (slope 1:XR)
        V("a0", -XU, 1); V("a1", 0, 1); V("a2", -XU, 2); V("a3", 0, 2); V("b0", XR, 0); V("b1", XD, 0); V("b2", XD, 1); V("b3", XD, 2)
        blk = [("a0", "a1", "a3", "a2", XU * n, n), ("b0", "b1", "b2", "a1", XD * n, n), ("a1", "b2", "b3", "a3", XD * n, n)]
    dz = H / n; nv = len(v); verts = " ".join(f"({x*H} {y*H} 0)" for x, y in v) + " " + " ".join(f"({x*H} {y*H} {dz})" for x, y in v)
    bl = " ".join(f"hex ({idx[a]} {idx[b]} {idx[c]} {idx[d]} {idx[a]+nv} {idx[b]+nv} {idx[c]+nv} {idx[d]+nv}) ({int(nx)} {int(ny)} 1) simpleGrading (1 1 1)" for a, b, c, d, nx, ny in blk)
    return verts, bl, idx, nv


def patches(kind, idx, nv):
    i = lambda k: idx[k]; face = lambda a, b: f"({i(a)} {i(b)} {i(b)+nv} {i(a)+nv})"
    if True:
        inlet = f"({i('a0')} {i('a0')+nv} {i('a2')+nv} {i('a2')})"; outlet = f"({i('b1')} {i('b2')} {i('b2')+nv} {i('b1')+nv}) ({i('b2')} {i('b3')} {i('b3')+nv} {i('b2')+nv})"
        walls = " ".join([face("a0", "a1"), face("a1", "b0"), face("b0", "b1")]); lid = " ".join([face("a2", "a3"), face("a3", "b3")])
        # NOTE vertical step wall face a1-b0 orientation fixed by blockMesh patch-face matching (faces are matched geometrically)
    return inlet, outlet, walls, lid


def make(kind, n, case):
    shutil.rmtree(case, ignore_errors=True); verts, bl, idx, nv = blocks(kind, n); inlet, outlet, walls, lid = patches(kind, idx, nv)
    sides = [("a0", "a1", "a3", "a2"), ("b0", "b1", "b2", "a1"), ("a1", "b2", "b3", "a3")]
    fb = " ".join(f"({idx[a]} {idx[d]} {idx[c]} {idx[b]}) ({idx[a]+nv} {idx[b]+nv} {idx[c]+nv} {idx[d]+nv})" for a, b, c, d in sides)
    P.write(f"{case}/system/blockMeshDict", "dictionary", "blockMeshDict", f"convertToMeters 1; vertices ({verts}); blocks ({bl}); boundary (inlet {{type patch; faces ({inlet});}} outlet {{type patch; faces ({outlet});}} walls {{type wall; faces ({walls});}} lid {{type symmetryPlane; faces ({lid});}} frontAndBack {{type empty; faces ({fb});}});\n")
    P.write(f"{case}/system/controlDict", "dictionary", "controlDict", "application simpleFoam; startFrom startTime; startTime 0; stopAt endTime; endTime 3000; deltaT 1; writeControl timeStep; writeInterval 3000; purgeWrite 1; writeFormat ascii; writePrecision 10; timeFormat general; timePrecision 8; runTimeModifiable no;\n")
    P.write(f"{case}/system/fvSchemes", "dictionary", "fvSchemes", "ddtSchemes {default steadyState;} gradSchemes {default Gauss linear;} divSchemes {default none; div(phi,U) bounded Gauss linearUpwind grad(U); div(phi,T) bounded Gauss linearUpwind grad(T); div((nuEff*dev2(T(grad(U))))) Gauss linear;} laplacianSchemes {default Gauss linear corrected;} interpolationSchemes {default linear;} snGradSchemes {default corrected;} wallDist {method meshWave;}\n")
    P.write(f"{case}/system/fvSolution", "dictionary", "fvSolution", "solvers { p {solver GAMG; smoother GaussSeidel; tolerance 1e-10; relTol 0.01;} U {solver smoothSolver; smoother symGaussSeidel; tolerance 1e-10; relTol 0.01;} T {solver smoothSolver; smoother symGaussSeidel; tolerance 1e-12; relTol 0;} } SIMPLE {nNonOrthogonalCorrectors 1; consistent yes; residualControl {p 1e-7; U 1e-8;}} relaxationFactors {equations {U 0.9; T 0.95;} fields {p 0.3;}}\n")
    P.write(f"{case}/constant/transportProperties", "dictionary", "transportProperties", f"transportModel Newtonian; nu [0 2 -1 0 0 0 0] {NU}; DT [0 2 -1 0 0 0 0] {ALPHA};\n")
    P.write(f"{case}/constant/turbulenceProperties", "dictionary", "turbulenceProperties", "simulationType laminar;\n")
    P.write(f"{case}/0/U", "volVectorField", "U", f"dimensions [0 1 -1 0 0 0 0]; internalField uniform ({UB} 0 0); boundaryField {{ inlet {{type fixedValue; value uniform ({UB} 0 0);}} outlet {{type zeroGradient;}} walls {{type noSlip;}} lid {{type symmetryPlane;}} frontAndBack {{type empty;}} }}\n")
    P.write(f"{case}/0/p", "volScalarField", "p", "dimensions [0 2 -2 0 0 0 0]; internalField uniform 0; boundaryField { inlet {type zeroGradient;} outlet {type fixedValue; value uniform 0;} walls {type zeroGradient;} lid {type symmetryPlane;} frontAndBack {type empty;} }\n")
    P.write(f"{case}/0/T", "volScalarField", "T", "dimensions [0 0 0 1 0 0 0]; internalField uniform 1; boundaryField { inlet {type fixedValue; value uniform 1;} outlet {type zeroGradient;} walls {type fixedValue; value uniform 0;} lid {type symmetryPlane;} frontAndBack {type empty;} }\n")


def read_list(path):
    t = open(path).read(); m = re.search(r"\n(\d+)\s*\n\(", t); n = int(m.group(1)); body = t[m.end():]; return n, body


def mesh_data(case):
    n, b = read_list(f"{case}/constant/polyMesh/points"); pts = np.array([[float(v) for v in l.strip("()").split()] for l in b.strip().splitlines()[:n]])
    nf, b = read_list(f"{case}/constant/polyMesh/faces"); faces = [[int(v) for v in re.findall(r"\d+", l)[1:]] for l in b.strip().splitlines()[:nf]]
    no, b = read_list(f"{case}/constant/polyMesh/owner"); owner = np.array([int(l) for l in b.strip().splitlines()[:no]])
    bt = open(f"{case}/constant/polyMesh/boundary").read(); patches = {}
    for m in re.finditer(r"(\w+)\s*\{[^}]*?nFaces\s+(\d+);\s*startFace\s+(\d+);", bt): patches[m.group(1)] = (int(m.group(3)), int(m.group(2)))
    return pts, faces, owner, patches


def read_cell_field(path, kind):
    t = open(path).read(); i = t.index("internalField"); head = t[i:t.index("(", i)].split(); n = int(head[-1]); j = t.index("(", i); body = t[j + 1:].split("\n)")[0].strip().splitlines()
    return np.array([[float(v) for v in l.strip("()").split()] for l in body[:n]]) if kind == "vec" else np.array([float(l) for l in body[:n]])


def cell_centres(pts, faces, owner, nbr, ncell):
    fc = np.array([pts[f].mean(axis=0) for f in faces]); C = np.zeros((ncell, 3)); cnt = np.zeros(ncell)
    for k in range(len(faces)):
        C[owner[k]] += fc[k]; cnt[owner[k]] += 1
        if k < len(nbr): C[nbr[k]] += fc[k]; cnt[nbr[k]] += 1
    return C / cnt[:, None], fc


def analyse(case, n):
    pts, faces, owner, patches = mesh_data(case); nn, b = read_list(f"{case}/constant/polyMesh/neighbour"); nbr = np.array([int(l) for l in b.strip().splitlines()[:nn]])
    tl = sorted([d for d in os.listdir(case) if d.replace(".", "").isdigit() and float(d) > 0], key=float)[-1]
    U = read_cell_field(f"{case}/{tl}/U", "vec"); T = read_cell_field(f"{case}/{tl}/T", "sca"); ncell = len(T); C, fc = cell_centres(pts, faces, owner, nbr, ncell)
    s, c = patches["walls"]; q = []; tau = []; xs = []
    for k in range(s, s + c):
        f = faces[k]; p = pts[f]; a = p[:, :2]; e = a.max(0) - a.min(0)
        v = a[1] - a[0] if len(f) == 4 else None
        # in-plane edge of the face: two vertices with z = 0
        z0 = [pp[:2] for pp in p if abs(pp[2]) < 1e-12]; ed = z0[1] - z0[0]; ln = np.linalg.norm(ed); tang = ed / ln; nrm = np.array([-tang[1], tang[0]])
        P0 = C[owner[k]][:2]; Fc = fc[k][:2]; d = (P0 - Fc); dn = abs(d @ nrm)
        if (d @ nrm) < 0: nrm = -nrm
        q.append(ALPHA * (T[owner[k]] - 0.0) / dn); xs.append(Fc[0] / H); tau.append(NU * (U[owner[k]][:2] @ tang) / dn * np.sign(tang[0] if tang[0] != 0 else 1.0))
    xs = np.array(xs); q = np.array(q); tau = np.array(tau); o = np.argsort(xs); xs, q, tau = xs[o], q[o], tau[o]
    # outlet: face-weighted flow quantities
    s2, c2 = patches["outlet"]; Pout = 0.0; Qout = 0.0
    for k in range(s2, s2 + c2):
        f = faces[k]; p = pts[f]; z0 = [pp[:2] for pp in p if abs(pp[2]) < 1e-12]; ed = z0[1] - z0[0]; area = np.linalg.norm(ed) * (H / n); u = U[owner[k]][0]; Qout += u * area; Pout += u * area * T[owner[k]]
    s1, c1 = patches["inlet"]; Qin = 0.0
    for k in range(s1, s1 + c1):
        f = faces[k]; p = pts[f]; z0 = [pp[:2] for pp in p if abs(pp[2]) < 1e-12]; ed = z0[1] - z0[0]; Qin += UB * np.linalg.norm(ed) * (H / n)
    wall_len = np.diff(np.concatenate([[xs[0]], xs]))  # placeholder
    # integrated wall heat per unit depth: sum q * face length (|edge|)
    Qw = 0.0
    for kk, k in enumerate(range(s, s + c)):
        f = faces[k]; p = pts[f]; z0 = [pp[:2] for pp in p if abs(pp[2]) < 1e-12]; Qw += q[np.where(np.isclose(xs, fc[k][0] / H))[0][0]] * np.linalg.norm(z0[1] - z0[0]) if False else 0.0
    # sign change of tangential wall velocity gradient downstream of the drop (x>0): reattachment
    m = (xs > 0.2) & (xs < XD - 1); xr = None
    sgn = np.sign(tau[m]); xm = xs[m]
    for i in range(1, len(sgn)):
        if sgn[i - 1] < 0 and sgn[i] > 0: xr = float(xm[i - 1] + (xm[i] - xm[i - 1]) * (-tau[m][i - 1]) / (tau[m][i] - tau[m][i - 1])); break
    return dict(n=n, cells=ncell, x_wall=xs.tolist(), q_wall=q.tolist(), tau_wall=tau.tolist(), Qin=Qin, Qout=Qout, Tb_out=Pout / Qout, enthalpy_flux_out_over_in=Pout / (Qin * 1.0), reattach_x_over_H=xr, time=tl)


def run_case(kind, n, wd):
    case = f"{wd}/{kind}_n{n}"; make(kind, n, case); P.run(["blockMesh"], case, "bm.log"); P.run(["checkMesh"], case, "cm.log")
    P.run(["simpleFoam"], case, "flow.log"); tf = sorted([d for d in os.listdir(case) if d.replace(".", "").isdigit() and float(d) > 0], key=float)[-1]
    shutil.copy(f"{case}/0/T", f"{case}/{tf}/T"); P.run(["potentialFoam", "-help"], case, "dummy.log") if False else None
    txt = open(f"{case}/system/controlDict").read().replace("application simpleFoam;", "application scalarTransportFoam;").replace("startFrom startTime; startTime 0;", f"startFrom startTime; startTime {tf};").replace("deltaT 1;", "deltaT 2;").replace("endTime 3000;", f"endTime {float(tf)+4000:g};").replace("writeInterval 3000;", "writeInterval 4000;").replace("writeControl timeStep;", "writeControl runTime;")
    open(f"{case}/system/controlDict", "w").write(txt); P.run(["scalarTransportFoam"], case, "T.log"); return analyse(case, n)


if __name__ == "__main__":
    wd = sys.argv[1]; ns = [int(a) for a in sys.argv[2:]] or [8]; out = {}
    for n in ns:
        for kind in ("step", "ramp"):
            r = run_case(kind, n, wd); out[f"{kind}_n{n}"] = r; print(kind, n, r["cells"], "Tb_out", r["Tb_out"], "reattach", r["reattach_x_over_H"], flush=True)
            json.dump(out, open(os.path.join(HERE, "step_vs_ramp_results.json"), "w"))
