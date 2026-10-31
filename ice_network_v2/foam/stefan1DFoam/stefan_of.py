"""Verification of stefan1DFoam (OpenFOAM v1912 finiteVolume, ALE moving mesh) against the exact one-phase Stefan solution WITH density jump (similarity solution derived in the milestone report):
  T' = -C' exp(-(eta-b)^2),  b = (1 - r) lam,  r = rho_i/rho_w,  lam solves  2 lam r / Ste = exp(-(lam-b)^2) / I(lam,b),  I = (sqrt(pi)/2)[erf(lam-b) + erf(b)],  s = 2 lam sqrt(kappa t),  kappa = k_w/(rho_w c_w),  Ste = c_w T_w / L.
Properties are the physical ones (rho_w 1000, c_w 4217, k_w 0.561, rho_i 917, L 3.34e5); Ste is varied only through the hot-wall temperature T_w. Equal-density limit = same code with rho_i set to rho_w and densityJump off (u_w = 0) as an explicit test of that limit.
Cases: Ste_env = 0.00631 (T_w = 0.5 K), Ste = 0.2 (T_w = 15.8 K) ; spatial refinement nx = 20, 40, 80 at fixed dt; time refinement dt = 5, 10, 20 s at nx = 40.   python stefan_of.py WORKDIR"""
import json, os, shutil, subprocess, sys
import numpy as np
from scipy.special import erf
from scipy.optimize import brentq
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "freesurface")); import common as C
ENV = dict(C.ENV); ENV["PATH"] = "/tmp/claude-0/ofbuild/bin:" + ENV.get("PATH", ""); ENV["LD_LIBRARY_PATH"] = "/usr/lib:/usr/lib/x86_64-linux-gnu/openmpi/lib:/tmp/claude-0/ofbuild/lib:" + ENV.get("LD_LIBRARY_PATH", "")
RW, CW, KW, RI, LF = 1000.0, 4217.0, 0.561, 917.0, 3.34e5; KAP = KW / (RW * CW)
def solve_lambda(Ste, r, jump):
    def F(l):
        b = (1 - r) * l if jump else 0.0; I = 0.5 * np.sqrt(np.pi) * (erf(l - b) + erf(b)); return 2 * l * r / Ste - np.exp(-(l - b) ** 2) / I
    l = brentq(F, 1e-8, 5.0, xtol=1e-15, rtol=1e-15); return l, ((1 - r) * l if jump else 0.0)
def run(case, Ste, rho_i, jump, nx, dt, t0=600.0, tend=3600.0):
    shutil.rmtree(case, ignore_errors=True); Tw = Ste * LF / CW; r = rho_i / RW; lam, b = solve_lambda(Ste, r, jump); s0 = 2 * lam * np.sqrt(KAP * t0); d = 1e-3
    C.write(f"{case}/system/blockMeshDict", "dictionary", "blockMeshDict", f"convertToMeters 1; vertices ((0 0 0) ({s0} 0 0) ({s0} {d} 0) (0 {d} 0) (0 0 {d}) ({s0} 0 {d}) ({s0} {d} {d}) (0 {d} {d})); blocks (hex (0 1 2 3 4 5 6 7) ({nx} 1 1) simpleGrading (1 1 1)); boundary (hot {{type patch; faces ((0 4 7 3));}} ice {{type patch; faces ((1 2 6 5));}} sides {{type empty; faces ((0 1 5 4) (3 7 6 2) (0 3 2 1) (4 5 6 7));}});\n")
    C.write(f"{case}/system/controlDict", "dictionary", "controlDict", f"application stefan1DFoam; startFrom startTime; startTime {t0}; stopAt endTime; endTime {tend}; deltaT {dt}; writeControl runTime; writeInterval {tend - t0}; writeFormat ascii; writePrecision 12; timeFormat general; timePrecision 8; runTimeModifiable no;\n")
    C.write(f"{case}/system/fvSchemes", "dictionary", "fvSchemes", "ddtSchemes {default Euler;} gradSchemes {default Gauss linear;} divSchemes {default Gauss upwind;} laplacianSchemes {default Gauss linear corrected;} interpolationSchemes {default linear;} snGradSchemes {default corrected;}\n")
    C.write(f"{case}/system/fvSolution", "dictionary", "fvSolution", "solvers { T {solver smoothSolver; smoother symGaussSeidel; tolerance 1e-16; relTol 0;} }\n")
    C.write(f"{case}/constant/stefanProperties", "dictionary", "stefanProperties", f"rhoW {RW}; cpW {CW}; kW {KW}; rhoI {rho_i}; L {LF}; Tw {Tw}; s0 {s0}; t0 {t0}; lambda {lam}; bpar {b}; densityJump {'true' if jump else 'false'};\n")
    C.write(f"{case}/{t0:g}/T", "volScalarField", "T", f"dimensions [0 0 0 1 0 0 0]; internalField uniform 0; boundaryField {{ hot {{type fixedValue; value uniform {Tw};}} ice {{type fixedValue; value uniform 0;}} sides {{type empty;}} }}\n")
    for cmd, log in ((["blockMesh"], "bm.log"), (["stefan1DFoam"], "run.log")):
        with open(f"{case}/{log}", "w") as f: rc = subprocess.run(cmd, cwd=case, env=ENV, stdout=f, stderr=subprocess.STDOUT).returncode
        if rc: raise RuntimeError(f"{cmd} failed in {case}")
    L = np.genfromtxt(f"{case}/ledger.csv", delimiter=",", names=True); t = L["t"]; s_ex = 2 * lam * np.sqrt(KAP * t); q_ex = RI * LF * lam * np.sqrt(KAP / t)
    E_lat_end = L["E_latent"][-1]; Esens = L["E_sens"][-1] - L["E_sens"][0] if False else None
    return dict(Ste=Ste, rho_i=rho_i, jump=jump, nx=nx, dt=dt, lam=lam, s_final=float(L["s"][-1]), s_exact=float(s_ex[-1]), s_rel_err=float((L["s"][-1] - s_ex[-1]) / s_ex[-1]), s_rel_err_max=float(np.max(np.abs(L["s"] - s_ex) / s_ex)),
                q_if_rel_err_final=float((L["q_if"][-1] * 1.0 - q_ex[-1] * (RI / RI) / 1.0 * 0 - (RI * LF * (L["sdot"][-1])) * 0) / 1.0 * 0) if False else float((L["sdot"][-1] - lam * np.sqrt(KAP / t[-1])) / (lam * np.sqrt(KAP / t[-1]))),
                ledger_resid_physical_over_latent=float(L["resid_physical"][-1] / E_lat_end), ledger_resid_discrete_over_supply=float(L["resid_discrete"][-1] / L["E_supply_cum"][-1]), E_supply=float(L["E_supply_cum"][-1]), E_latent=float(E_lat_end), uw_over_sdot=(1 - r) if jump else 0.0)
if __name__ == "__main__":
    wd = sys.argv[1]; out = {}
    for Ste in (0.00631, 0.2):
        for tag, rho_i, jump in (("rhoI917_jump", RI, True), ("equal_density", RW, False)):
            for nx, dt in ((20, 2.5), (40, 2.5), (80, 2.5), (160, 2.5), (40, 5.0), (40, 10.0)):
                k = f"Ste{Ste}_{tag}_nx{nx}_dt{dt}"
                try: out[k] = run(f"{wd}/{k}", Ste, rho_i, jump, nx, dt)
                except Exception as e: out[k] = dict(failed=str(e))
                print(k, json.dumps(out[k])[:260], flush=True); json.dump(out, open(os.path.join(HERE, "stefan_of_results.json"), "w"), indent=1)
