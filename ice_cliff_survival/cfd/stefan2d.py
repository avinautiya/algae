"""Driver, analytic references and ledger analysis for stefan2DFoam (see cfd/VERIFICATION.md). Properties are physical (rho_w 1000, c_w 4217, k_w 0.561, rho_i 917 or equal-density test, L 3.34e5); Ste enters only through T_w."""
import json, math, os, shutil, subprocess, sys, time
import numpy as np
from scipy.special import erf
from scipy.optimize import brentq
HERE = os.path.dirname(os.path.abspath(__file__))
ENV = {k: v for k, v in os.environ.items() if k != "FOAM_SIGFPE"}; ENV.update(WM_PROJECT_DIR="/usr/share/openfoam", FOAM_ETC="/usr/share/openfoam/etc", PATH="/tmp/claude-0/ofbuild/bin:" + os.environ.get("PATH", ""), LD_LIBRARY_PATH="/usr/lib:/usr/lib/x86_64-linux-gnu/openmpi/lib:/tmp/claude-0/ofbuild/lib:" + os.environ.get("LD_LIBRARY_PATH", ""))
RW, CW, KW, RI0, LF = 1000.0, 4217.0, 0.561, 917.0, 3.34e5; KAP = KW / (RW * CW)
HDR = 'FoamFile { version 2.0; format ascii; class %s; object %s; }\n'
def W(path, cls, obj, body): os.makedirs(os.path.dirname(path), exist_ok=True); open(path, "w").write(HDR % (cls, obj) + body)
def solve_lambda(Ste, r, jump):
    def F(l):
        b = (1 - r) * l if jump else 0.0; I = 0.5 * math.sqrt(math.pi) * (erf(l - b) + erf(b)); return 2 * l * r / Ste - math.exp(-(l - b) ** 2) / I
    l = brentq(F, 1e-8, 5.0, xtol=1e-15, rtol=1e-15); return l, ((1 - r) * l if jump else 0.0)
def fmt(a): a = list(map(float, a)); return "%d(%s)" % (len(a), " ".join("%.15g" % x for x in a))
def run(case, Ste, nx, ny, dt, tend, rho_i=RI0, jump=True, eps=0.0, Ly=0.05, s0=1.0e-3, hShield=None, dz=1e-3, timeout=3600, mode="cos"):
    shutil.rmtree(case, ignore_errors=True); r = rho_i / RW; Tw0 = Ste * LF / CW
    dyv = Ly / ny; yc = (np.arange(ny) + 0.5) * dyv
    Tw = Tw0 * (1.0 + eps * np.cos(math.pi * yc / Ly)) if mode == "cos" else np.full(ny, Tw0)
    lam = np.zeros(ny); b = np.zeros(ny)
    for j in range(ny): lam[j], b[j] = solve_lambda(Tw[j] * CW / LF, r, jump)
    tRow = (s0 / (2 * lam)) ** 2 / KAP; hS = np.full(ny, 1e30) if hShield is None else np.array(hShield, float)
    W(f"{case}/system/blockMeshDict", "dictionary", "blockMeshDict", f"convertToMeters 1; vertices ((0 0 0) ({s0} 0 0) ({s0} {Ly} 0) (0 {Ly} 0) (0 0 {dz}) ({s0} 0 {dz}) ({s0} {Ly} {dz}) (0 {Ly} {dz})); blocks (hex (0 1 2 3 4 5 6 7) ({nx} {ny} 1) simpleGrading (1 1 1)); boundary (hot {{type patch; faces ((0 4 7 3));}} ice {{type patch; faces ((1 2 6 5));}} sides {{type patch; faces ((0 1 5 4) (3 7 6 2));}} frontAndBack {{type empty; faces ((0 3 2 1) (4 5 6 7));}});\n")
    W(f"{case}/system/controlDict", "dictionary", "controlDict", f"application stefan2DFoam; startFrom startTime; startTime 0; stopAt endTime; endTime {tend}; deltaT {dt}; writeControl runTime; writeInterval {tend}; writeFormat ascii; writePrecision 12; timeFormat general; timePrecision 8; runTimeModifiable no;\n")
    W(f"{case}/system/fvSchemes", "dictionary", "fvSchemes", "ddtSchemes {default Euler;} gradSchemes {default Gauss linear;} divSchemes {default Gauss upwind;} laplacianSchemes {default Gauss linear corrected;} interpolationSchemes {default linear;} snGradSchemes {default corrected;}\n")
    W(f"{case}/system/fvSolution", "dictionary", "fvSolution", "solvers { T {solver PBiCGStab; preconditioner DILU; tolerance 1e-15; relTol 0; maxIter 20000;} }\n")
    W(f"{case}/constant/stefanProperties", "dictionary", "stefanProperties", f"rhoW {RW}; cpW {CW}; kW {KW}; rhoI {rho_i}; L {LF}; s0 {s0}; t0 0; Ly {Ly}; dz {dz}; ny {ny}; densityJump {'true' if jump else 'false'};\n tRow {fmt(tRow)}; TwRow {fmt(Tw)}; lambdaRow {fmt(lam)}; bRow {fmt(b)}; hShield {fmt(hS)};\n")
    W(f"{case}/0/T", "volScalarField", "T", "dimensions [0 0 0 1 0 0 0]; internalField uniform 0; boundaryField { hot {type fixedValue; value uniform 0;} ice {type mixed; refValue uniform 0; refGradient uniform 0; valueFraction uniform 1; value uniform 0;} sides {type zeroGradient;} frontAndBack {type empty;} }\n")
    for cmd, log in ((["blockMesh"], "bm.log"), (["stefan2DFoam"], "run.log")):
        t = time.time()
        with open(f"{case}/{log}", "w") as f: rc = subprocess.run(cmd, cwd=case, env=ENV, stdout=f, stderr=subprocess.STDOUT, timeout=timeout).returncode
        if rc: raise RuntimeError(f"{cmd} failed in {case}")
        if cmd[0] == "stefan2DFoam": wall = time.time() - t
    L = np.genfromtxt(f"{case}/ledger.csv", delimiter=",", names=True); pr = open(f"{case}/{tend:g}/interface.csv").read().strip().splitlines()[1:]; prof = np.array([[float(x) for x in l.split(",")] for l in pr])
    return dict(L=L, interface=prof, lam=lam, b=b, tRow=tRow, Tw=Tw, wall_s=wall, cells=nx * ny, params=dict(Ste=Ste, nx=nx, ny=ny, dt=dt, tend=tend, rho_i=rho_i, jump=jump, eps=eps, Ly=Ly, s0=s0, dz=dz))
def neumann_s(res, t): return 2 * res["lam"] * np.sqrt(KAP * (t + res["tRow"]))
def shielded_s(Tw, h, s0, t, k=KW, rho_i=RI0, L=LF):
    """quasi-steady (Ste -> 0) Stefan with a static conductive shield: q = Tw/(1/h + s/k); s/h + s^2/(2k) = s0/h + s0^2/(2k) + Tw t/(rho_i L)"""
    rhs = (s0 / h + s0 ** 2 / (2 * k)) + Tw * t / (rho_i * L); return -k / h + math.sqrt(k * k / h / h + 2 * k * rhs)
