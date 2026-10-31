"""Time-step and linear-solver sensitivity at ny=20 (same case as poiseuille_graetz.py): dt/2, dt/4 and a tighter pressure tolerance. Fields kept in workdir for graetz_exact.py."""
import os, shutil, sys, re
def sub(path, a, b):
    t = open(path).read().replace(a, b); open(path, "w").write(t)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import poiseuille_graetz as P
wd = sys.argv[1]; ny = 20
for tag, dt, tight in (("dt0.05", 0.05, False), ("dt0.025", 0.025, False), ("tol", 0.1, True)):
    case = f"{wd}/ny20_{tag}"; shutil.rmtree(case, ignore_errors=True)
    P.make(case, ny, dt, 300, "icoFoam"); P.fields_u(case)
    if tight: sub(f"{case}/system/fvSolution", "relTol 0.01;", "relTol 0.0001;"); sub(f"{case}/system/fvSolution", "nCorrectors 2;", "nCorrectors 4;")
    P.run(["blockMesh"], case, "bm.log"); P.run(["icoFoam"], case, "ico.log")
    P.make(case, ny, dt, 300, "scalarTransportFoam"); tdir = sorted([d for d in os.listdir(case) if d.replace(".", "").isdigit()], key=float)[-1]; P.fields_T(case)
    shutil.copy(f"{case}/0/T", f"{case}/{tdir}/T"); sub(f"{case}/system/controlDict", "endTime 300;", "endTime 600;")
    P.run(["scalarTransportFoam"], case, "T.log"); print(tag, "done", flush=True)
