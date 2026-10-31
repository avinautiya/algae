"""Fixed-geometry laminar verification with stock OpenFOAM v1912 (ESI, Ubuntu package 1912.200626-2build3): 2-D plane channel, height H, length L = 20 H.
 (1) icoFoam: developed plane Poiseuille profile u(y) = 6 Ub (y/H)(1 - y/H) against the exact parabola, three resolutions.
 (2) scalarTransportFoam on the converged velocity: Graetz problem with wall temperature 0 (melting-point reference), inlet T = 1: developed Nu_Dh -> 7.5407 (isothermal parallel plates, Dh = 2H).
This is a verification of the standard solvers on a fixed wall. It has no free surface, no phase change, no ice. Run:  python poiseuille_graetz.py [workdir]"""
import os, subprocess, sys, json, shutil
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ENV = dict(os.environ, WM_PROJECT_DIR="/usr/share/openfoam", FOAM_ETC="/usr/share/openfoam/etc")
H = 0.01; L = 20 * H; NU = 1e-6; UB = 1e-3; PR = 7.0; DT_T = NU / PR          # Re_H = 10, Pe_H = 70
HDR = 'FoamFile { version 2.0; format ascii; class %s; %sobject %s; }\n'


def write(path, cls, obj, body, loc=""):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    open(path, "w").write(HDR % (cls, ('location "%s"; ' % loc) if loc else "", obj) + body)


def make(case, ny, dt, tend, solver):
    nx = 20 * ny
    write(f"{case}/system/blockMeshDict", "dictionary", "blockMeshDict", f"""convertToMeters 1;
vertices ((0 0 0) ({L} 0 0) ({L} {H} 0) (0 {H} 0) (0 0 {H/ny}) ({L} 0 {H/ny}) ({L} {H} {H/ny}) (0 {H} {H/ny}));
blocks (hex (0 1 2 3 4 5 6 7) ({nx} {ny} 1) simpleGrading (1 1 1));
boundary (inlet {{type patch; faces ((0 4 7 3));}} outlet {{type patch; faces ((1 2 6 5));}} walls {{type wall; faces ((0 1 5 4) (3 7 6 2));}} frontAndBack {{type empty; faces ((0 3 2 1) (4 5 6 7));}});
""")
    write(f"{case}/system/controlDict", "dictionary", "controlDict", f"""application {solver}; startFrom latestTime; startTime 0; stopAt endTime; endTime {tend}; deltaT {dt};
writeControl runTime; writeInterval {tend}; purgeWrite 0; writeFormat ascii; writePrecision 10; timeFormat general; timePrecision 8; runTimeModifiable no;
""")
    write(f"{case}/system/fvSchemes", "dictionary", "fvSchemes", """ddtSchemes {default Euler;} gradSchemes {default Gauss linear;} divSchemes {default none; div(phi,U) Gauss linear; div(phi,T) Gauss linear;}
laplacianSchemes {default Gauss linear corrected;} interpolationSchemes {default linear;} snGradSchemes {default corrected;}
""")
    write(f"{case}/system/fvSolution", "dictionary", "fvSolution", """solvers { p {solver PCG; preconditioner DIC; tolerance 1e-12; relTol 0.01;} pFinal {$p; relTol 0;} U {solver smoothSolver; smoother symGaussSeidel; tolerance 1e-12; relTol 0;} T {solver smoothSolver; smoother symGaussSeidel; tolerance 1e-14; relTol 0;} }
PISO {nCorrectors 2; nNonOrthogonalCorrectors 0; pRefCell 0; pRefValue 0;} SIMPLE {nNonOrthogonalCorrectors 0;}
""")
    write(f"{case}/constant/transportProperties", "dictionary", "transportProperties", f"nu [0 2 -1 0 0 0 0] {NU}; DT [0 2 -1 0 0 0 0] {DT_T};\n")


def fields_u(case):
    write(f"{case}/0/U", "volVectorField", "U", f"dimensions [0 1 -1 0 0 0 0]; internalField uniform (0 0 0); boundaryField {{ inlet {{type fixedValue; value uniform ({UB} 0 0);}} outlet {{type zeroGradient;}} walls {{type noSlip;}} frontAndBack {{type empty;}} }}\n")
    write(f"{case}/0/p", "volScalarField", "p", "dimensions [0 2 -2 0 0 0 0]; internalField uniform 0; boundaryField { inlet {type zeroGradient;} outlet {type fixedValue; value uniform 0;} walls {type zeroGradient;} frontAndBack {type empty;} }\n")


def fields_T(case):
    write(f"{case}/0/T", "volScalarField", "T", "dimensions [0 0 0 1 0 0 0]; internalField uniform 0; boundaryField { inlet {type fixedValue; value uniform 1;} outlet {type zeroGradient;} walls {type fixedValue; value uniform 0;} frontAndBack {type empty;} }\n")


def run(cmd, case, log):
    with open(os.path.join(case, log), "w") as f: r = subprocess.run(cmd, cwd=case, env=ENV, stdout=f, stderr=subprocess.STDOUT)
    if r.returncode: raise RuntimeError(f"{cmd} failed in {case}: see {log}")


def read_field(path, kind):
    txt = open(path).read(); i = txt.index("internalField"); j = txt.index("(", i); n = int(txt[txt.rindex("\n", 0, j) - 12:j].split()[-1]) if False else None
    head = txt[i:j].split(); n = int(head[-1]); body = txt[j + 1:].split("\n)")[0].strip().splitlines()
    return np.array([[float(v) for v in l.strip("()").split()] for l in body[:n]]) if kind == "vec" else np.array([float(l) for l in body[:n]])


def analyse(case, ny):
    nx = 20 * ny; tl = sorted([d for d in os.listdir(case) if d.replace(".", "").isdigit() and float(d) > 0], key=float)[-1]
    # uniform structured mesh, cell index = i + nx*j (x fastest): centres analytic (postProcess writeCellCentres crashes in this Debian build: sha1 IO error, logged cc.log in early attempt)
    dxm = L / nx; dym = H / ny; X = np.tile((np.arange(nx) + 0.5) * dxm, (ny, 1)); Y = np.tile(((np.arange(ny) + 0.5) * dym)[:, None], (1, nx))
    U = read_field(f"{case}/{tl}/U", "vec")[:, 0]; T = read_field(f"{case}/{tl}/T", "sca"); Uu = U.reshape(ny, nx); TT = T.reshape(ny, nx)
    i = int(0.9 * nx); y = Y[:, i]; u = Uu[:, i]; tt = TT[:, i]; ub = u.mean()
    ex = 6 * ub * (y / H) * (1 - y / H); err_u = float(np.sqrt(np.mean((u - ex) ** 2)) / ub); mass = float(np.mean(Uu[:, -1]) / UB)
    # developed Nusselt: T_b = flow-weighted mean; q_w from a 2nd-order one-sided wall gradient (cell values at y1,y2 with wall value 0)
    Tb = float((u * tt).sum() / u.sum()); dy = H / ny; y1, y2 = 0.5 * dy, 1.5 * dy; t1, t2 = tt[0], tt[1]
    dTdy = (t1 * y2 ** 2 - t2 * y1 ** 2) / (y1 * y2 * (y2 - y1))        # derivative at 0 of the quadratic through (0,0),(y1,t1),(y2,t2)
    Nu = float(2 * H * dTdy / (0 - Tb)) * -1.0
    return dict(ny=ny, nx=nx, cells=nx * ny, dy=dy, err_u_rel=err_u, outlet_mass_ratio=mass, Tb=Tb, Nu_Dh=abs(Nu), Nu_exact=7.5407)


if __name__ == "__main__":
    wd = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "work"); shutil.rmtree(wd, ignore_errors=True); out = []
    for ny, dt in ((10, 0.2), (20, 0.1), (40, 0.05)):
        case = f"{wd}/ny{ny}"; make(case, ny, dt, 300, "icoFoam"); fields_u(case); run(["blockMesh"], case, "bm.log"); run(["icoFoam"], case, "ico.log")
        # thermal: restart from developed U, solve T to steady state
        make(case, ny, dt, 300, "scalarTransportFoam"); 
        tdir = sorted([d for d in os.listdir(case) if d.replace(".", "").isdigit()], key=float)[-1]; fields_T(case)
        shutil.copy(f"{case}/0/T", f"{case}/{tdir}/T"); txt = open(f"{case}/system/controlDict").read().replace("startFrom latestTime;", "startFrom latestTime;").replace("endTime 300;", "endTime 600;"); open(f"{case}/system/controlDict", "w").write(txt)
        run(["scalarTransportFoam"], case, "T.log"); r = analyse(case, ny); r["dt"] = dt; out.append(r); print(r, flush=True)
    for a, b in zip(out[:-1], out[1:]):
        b["order_u"] = float(np.log2(a["err_u_rel"] / b["err_u_rel"])) if b["err_u_rel"] > 0 else None; b["order_Nu"] = float(np.log2(abs(a["Nu_Dh"] - 7.5407) / abs(b["Nu_Dh"] - 7.5407)))
    json.dump(dict(openfoam="v1912 patch 200626 (Ubuntu openfoam 1912.200626-2build3)", Re_H=UB * H / NU, Pe_H=UB * H / DT_T, Pr=PR, results=out), open(os.path.join(HERE, "poiseuille_graetz_results.json"), "w"), indent=1); print(json.dumps(out, indent=1))
