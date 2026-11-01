"""Driver for flowStefanFoam (M2): case generation, run, ledger parsing.  `run_case(case, **params)`; see PARAM_DEFAULTS.
Geometry: floor [inert Lu | ice Li | inert Ld], lid at H. Equal density when rho_i = rho_w (default for stage A); stage B uses rho_i = 917 with the LOCAL interface mass condition."""
import math, os, shutil, subprocess, time, resource
import numpy as np
ENV = dict(os.environ)
import shlex
_out = subprocess.run(["bash", "-c", "source /tmp/claude-0/ofenv.sh && env -0"], capture_output=True, text=True).stdout
ENV.update(dict(kv.split("=", 1) for kv in _out.split("\0") if "=" in kv))
RHO_W, CP_W, K_W, LF, MU = 1000.0, 4217.0, 0.561, 3.34e5, 1.0e-3
PARAM_DEFAULTS = dict(H=0.02, dz=0.01, Lu=0.04, Li=0.20, Ld=0.04, dx=0.004, ny=24, ygrade=1.0, Ubar=0.01, Tin=1.0, Tlid_fixed=None, rho_i=RHO_W, melt=True, Lmult=1.0, dt=0.5, tend=100.0, closedEnds=False,
                      lidClosure=False, maxCoupling=12, tolCoupling=1e-6, nPiso=3, dump=False, tdiv="linearUpwind grad(T)", udiv="linearUpwind grad(U)", writeEvery=None, ddtCorr=True, jumpBC=True, slopeMax=1.0, convU=True, momPred=False, outletOpen=False, nuMult=1.0, dumpEvery=1, Tinit=None, neumann=None, hShield=None, floorBC="isothermal")
def W(p, cls, obj, body):
    os.makedirs(os.path.dirname(p), exist_ok=True)
    open(p, "w").write(f"FoamFile {{version 2.0; format ascii; class {cls}; object {obj};}}\n{body}\n")
def make_case(case, **kw):
    P = dict(PARAM_DEFAULTS); P.update(kw); shutil.rmtree(case, ignore_errors=True)
    H, dz, dx = P["H"], P["dz"], P["dx"]; secs = [(n, L) for n, L in (("floorUp", P["Lu"]), ("ice", P["Li"]), ("floorDown", P["Ld"])) if L > 0]
    Ltot = sum(L for _, L in secs); nxs = [max(1, round(L / dx)) for _, L in secs]; nx = sum(nxs); Dx = Ltot / nx
    if abs(sum(nxs) * Dx - Ltot) > 1e-12: raise ValueError
    xs = [0.0]; 
    for (_, L) in secs: xs.append(xs[-1] + L)
    vt = [(x, 0.0, 0.0) for x in xs] + [(x, H, 0.0) for x in xs] + [(x, 0.0, dz) for x in xs] + [(x, H, dz) for x in xs]; nb = len(secs) + 1
    V = lambda i, k: i + k * nb       # k: 0 bottom-front(z0) , 1 top z0, 2 bottom z dz, 3 top z dz
    g = P["ygrade"]; grad = "1" if g == 1.0 else f"((0.5 0.5 {1.0/g:.12g}) (0.5 0.5 {g:.12g}))"       # first half: from floor to mid grows ratio 1/g? see below
    if g != 1.0: grad = f"((0.5 0.5 {g:.12g}) (0.5 0.5 {1.0/g:.12g}))"
    blocks = ""; bpatch = {n: [] for n, _ in secs}
    for b, ((n, L), nxb) in enumerate(zip(secs, nxs)):
        blocks += f"hex ({V(b,0)} {V(b+1,0)} {V(b+1,1)} {V(b,1)} {V(b,2)} {V(b+1,2)} {V(b+1,3)} {V(b,3)}) ({nxb} {P['ny']} 1) simpleGrading (1 {grad} 1)\n"
    faces = lambda b, which: {"bot": f"({V(b,0)} {V(b+1,0)} {V(b+1,2)} {V(b,2)})", "top": f"({V(b,1)} {V(b,3)} {V(b+1,3)} {V(b+1,1)})"}[which]
    pat = "inlet {type patch; faces ((%d %d %d %d));}\noutlet {type patch; faces ((%d %d %d %d));}\n" % (V(0,0), V(0,2), V(0,3), V(0,1), V(nb-1,0), V(nb-1,1), V(nb-1,3), V(nb-1,2))
    pat += "lid {type wall; faces (" + " ".join(faces(b, "top") for b in range(nb - 1)) + ");}\n"
    for b, (n, _) in enumerate(secs): pat += f"{n} {{type wall; faces ({faces(b,'bot')});}}\n"
    fr = " ".join("(%d %d %d %d)" % (V(b,0), V(b,1), V(b+1,1), V(b+1,0)) for b in range(nb - 1)); bk = " ".join("(%d %d %d %d)" % (V(b,2), V(b+1,2), V(b+1,3), V(b,3)) for b in range(nb - 1))
    pat += f"frontAndBack {{type empty; faces ({fr} {bk});}}\n"
    W(f"{case}/system/blockMeshDict", "dictionary", "blockMeshDict", "convertToMeters 1;\nvertices (" + " ".join("(%.15g %.15g %.15g)" % v for v in vt) + f");\nblocks ({blocks});\nedges ();\nboundary ({pat});\nmergePatchPairs ();")
    wi = P["writeEvery"] or P["tend"]
    W(f"{case}/system/controlDict", "dictionary", "controlDict", f"application flowStefanFoam; startFrom startTime; startTime 0; stopAt endTime; endTime {P['tend']}; deltaT {P['dt']}; writeControl runTime; writeInterval {wi}; writeFormat ascii; writePrecision 10; timeFormat general; timePrecision 8; runTimeModifiable no;")
    W(f"{case}/system/fvSchemes", "dictionary", "fvSchemes", f"ddtSchemes {{default Euler;}} gradSchemes {{default Gauss linear;}} divSchemes {{default none; div(phiRel,U) Gauss {P['udiv']}; div(phiEq,T) Gauss {P['tdiv']};}} laplacianSchemes {{default Gauss linear corrected;}} interpolationSchemes {{default linear;}} snGradSchemes {{default corrected;}} fluxRequired {{default no; p;}}")
    W(f"{case}/system/fvSolution", "dictionary", "fvSolution", "solvers { p {solver PCG; preconditioner DIC; tolerance 1e-13; relTol 0; maxIter 2000;} U {solver smoothSolver; smoother symGaussSeidel; tolerance 1e-12; relTol 1e-3; maxIter 300;} T {solver PBiCGStab; preconditioner DILU; tolerance 1e-14; relTol 0;} }")
    hs = "" if P["hShield"] is None else "hShield %d(%s);" % (len(P["hShield"]), " ".join("%.15g" % x for x in P["hShield"]))
    nm = "" if P["neumann"] is None else "neumannTw %.15g; neumannLam %.15g; neumannB %.15g; neumannT0 %.15g;" % tuple(P["neumann"])
    rho_i = P["rho_i"]; Lm = LF * P["Lmult"]
    W(f"{case}/constant/flowProperties", "dictionary", "flowProperties", f"rhoW {RHO_W}; cpW {CP_W}; kW {K_W}; nu {MU/RHO_W*P['nuMult']}; rhoI {rho_i}; L {Lm}; H {H}; dz {dz}; Lu {P['Lu']}; Li {P['Li']}; Ltot {Ltot}; Ubar {P['Ubar']}; Tin {P['Tin']}; melt {'true' if P['melt'] else 'false'}; closedEnds {'true' if P['closedEnds'] else 'false'}; lidClosure {'true' if P['lidClosure'] else 'false'}; dump {'true' if P['dump'] else 'false'}; nx {nx}; ny {P['ny']}; maxCoupling {P['maxCoupling']}; tolCoupling {P['tolCoupling']}; nPiso {P['nPiso']}; ddtCorr {'true' if P['ddtCorr'] else 'false'}; jumpBC {'true' if P['jumpBC'] else 'false'}; slopeMax {P['slopeMax']}; convU {'true' if P['convU'] else 'false'}; momPred {'true' if P['momPred'] else 'false'}; dumpEvery {P['dumpEvery']}; {hs} {nm}")
    W(f"{case}/constant/transportProperties", "dictionary", "transportProperties", "")
    wall0 = "type fixedValue; value uniform (0 0 0);"; zg = "type zeroGradient;"
    fl = {n: wall0 for n in ("floorUp", "floorDown", "lid")}
    ub = "inlet {type fixedValue; value uniform (0 0 0);} outlet {type zeroGradient;} lid {%s} ice {type fixedValue; value uniform (0 0 0);} floorUp {%s} floorDown {%s} frontAndBack {type empty;}" % (wall0, wall0, wall0)
    if P["closedEnds"] and not P["outletOpen"]: ub = ub.replace("outlet {type zeroGradient;}", "outlet {type fixedValue; value uniform (0 0 0);}")
    W(f"{case}/0/U", "volVectorField", "U", f"dimensions [0 1 -1 0 0 0 0]; internalField uniform (0 0 0); boundaryField {{ {ub} }}")
    pb = "inlet {type zeroGradient;} outlet {type fixedValue; value uniform 0;} lid {type zeroGradient;} ice {type zeroGradient;} floorUp {type zeroGradient;} floorDown {type zeroGradient;} frontAndBack {type empty;}"
    if P["closedEnds"] and not P["outletOpen"]: pb = pb.replace("outlet {type fixedValue; value uniform 0;}", "outlet {type zeroGradient;}")
    W(f"{case}/0/p", "volScalarField", "p", f"dimensions [0 2 -2 0 0 0 0]; internalField uniform 0; boundaryField {{ {pb} }}")
    lidT = f"type fixedValue; value uniform {P['Tlid_fixed']};" if P["Tlid_fixed"] is not None else zg
    flT = "type fixedValue; value uniform 0;" if P["floorBC"] == "isothermal" else zg
    outT = zg; inT = f"type fixedValue; value uniform {P['Tin']};"
    if P["closedEnds"]: inT = zg
    W(f"{case}/0/T", "volScalarField", "T", f"dimensions [0 0 0 1 0 0 0]; internalField uniform {P['Tin'] if P['Tinit'] is None else P['Tinit']}; boundaryField {{ inlet {{{inT}}} outlet {{{outT}}} lid {{{lidT}}} ice {{type mixed; refValue uniform 0; refGradient uniform 0; valueFraction uniform 1; value uniform 0;}} floorUp {{{flT}}} floorDown {{{flT}}} frontAndBack {{type empty;}} }}")
    return P, nx, Ltot
def run_case(case, timeout=7200, **kw):
    P, nx, Ltot = make_case(case, **kw); t0 = time.time(); r0 = resource.getrusage(resource.RUSAGE_CHILDREN)
    for cmd, log in ((["blockMesh"], "bm.log"), (["flowStefanFoam"], "run.log")):
        with open(f"{case}/{log}", "w") as f: rc = subprocess.run(cmd, cwd=case, env=ENV, stdout=f, stderr=subprocess.STDOUT, timeout=timeout).returncode
        if rc: raise RuntimeError(f"{cmd} failed in {case}: see {case}/{log}")
    wall = time.time() - t0; r1 = resource.getrusage(resource.RUSAGE_CHILDREN)
    L = np.genfromtxt(f"{case}/ledger.csv", delimiter=",", names=True)
    return dict(L=L, P=P, nx=nx, Ltot=Ltot, cells=nx * P["ny"], wall_s=wall, cpu_s=(r1.ru_utime + r1.ru_stime) - (r0.ru_utime + r0.ru_stime), maxrss_MB=r1.ru_maxrss / 1024.0)
