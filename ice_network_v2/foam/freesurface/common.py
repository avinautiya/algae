"""Shared case writers for stock-interFoam verification cases (OpenFOAM ESI v1912 patch 200626). Isothermal: interFoam solves NO temperature equation."""
import os, subprocess, sys
import numpy as np
ENV = {k: v for k, v in os.environ.items() if k != "FOAM_SIGFPE"}; ENV.update(WM_PROJECT_DIR="/usr/share/openfoam", FOAM_ETC="/usr/share/openfoam/etc")
# FOAM_SIGFPE (floating-point trapping) is unset for the interFoam runs: in this Debian-built v1912 the trap aborts inside PCG::scalarSolve (sumProd) at the first steps although the same case runs to completion without it
# (continuity errors, residuals and alpha bounds are checked in the analysis; NaN would still be detected there). Recorded in MILESTONE docs.
HDR = 'FoamFile { version 2.0; format ascii; class %s; object %s; }\n'
def write(path, cls, obj, body):
    os.makedirs(os.path.dirname(path), exist_ok=True); open(path, "w").write(HDR % (cls, obj) + body)
def run(cmd, case, log):
    with open(os.path.join(case, log), "w") as f: r = subprocess.run(cmd, cwd=case, env=ENV, stdout=f, stderr=subprocess.STDOUT)
    if r.returncode: raise RuntimeError(f"{cmd} failed in {case}, see {log}")
def physics(case, sigma=0.07, nu_w=1e-6):
    write(f"{case}/constant/transportProperties", "dictionary", "transportProperties", f"phases (water air); water {{transportModel Newtonian; nu {nu_w}; rho 1000;}} air {{transportModel Newtonian; nu 1.48e-05; rho 1;}} sigma {sigma};\n")
    write(f"{case}/constant/g", "uniformDimensionedVectorField", "g", "dimensions [0 1 -2 0 0 0 0]; value (0 -9.81 0);\n")
    write(f"{case}/constant/turbulenceProperties", "dictionary", "turbulenceProperties", "simulationType laminar;\n")
    write(f"{case}/system/fvSchemes", "dictionary", "fvSchemes", "ddtSchemes {default Euler;} gradSchemes {default Gauss linear;} divSchemes {div(rhoPhi,U) Gauss linearUpwind grad(U); div(phi,alpha) Gauss vanLeer; div(phirb,alpha) Gauss linear; div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;} laplacianSchemes {default Gauss linear corrected;} interpolationSchemes {default linear;} snGradSchemes {default corrected;}\n")
    write(f"{case}/system/fvSolution", "dictionary", "fvSolution", '''solvers { "alpha.water.*" {nAlphaCorr 1; nAlphaSubCycles 1; cAlpha 1; MULESCorr yes; nLimiterIter 3; solver smoothSolver; smoother symGaussSeidel; tolerance 1e-10; relTol 0;}
 "pcorr.*" {solver PCG; preconditioner DIC; tolerance 1e-8; relTol 0;} p_rgh {solver PCG; preconditioner DIC; tolerance 1e-9; relTol 0.01;} p_rghFinal {$p_rgh; relTol 0;} "(U|k|epsilon)" {solver smoothSolver; smoother symGaussSeidel; tolerance 1e-9; relTol 0;} "(U|k|epsilon)Final" {$U; relTol 0;} }
PIMPLE {momentumPredictor no; nOuterCorrectors 1; nCorrectors 3; nNonOrthogonalCorrectors 3;}
''')
def control(case, tend, dt0, maxco, wint, maxdt=0.01):
    write(f"{case}/system/controlDict", "dictionary", "controlDict", f"application interFoam; startFrom startTime; startTime 0; stopAt endTime; endTime {tend}; deltaT {dt0}; writeControl adjustableRunTime; writeInterval {wint}; purgeWrite 0; writeFormat ascii; writePrecision 10; timeFormat general; timePrecision 8; runTimeModifiable no; adjustTimeStep yes; maxCo {maxco}; maxAlphaCo {maxco}; maxDeltaT {maxdt};\n")
