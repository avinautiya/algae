"""Independent reconstruction of mass and energy budgets from the SAVED per-step fields of flowStefanFoam (dump.csv), using its own geometry (cell volumes, face areas and normals from the floor node heights b_i),
not the solver's accumulated ledger. Structured-column geometry: column i spans x in [i Dx, (i+1) Dx], floor nodes b_i, b_{i+1}, lid y=H, `ny` rows uniform in eta (use ygrade=1 for budgeted cases)."""
import numpy as np
from .flowstefan import RHO_W, CP_W, K_W, LF
def read_dump(path, nx, ny):
    steps = []; cur = None
    for ln in open(path):
        if ln.startswith("#"): continue
        t = ln.split()
        if t[0] == "S": cur = dict(n=int(t[1]), t=float(t[2]), T=np.zeros((nx, ny)), Ux=np.zeros((nx, ny)), Uy=np.zeros((nx, ny)), I=[], O=[], F=[]); steps.append(cur)
        elif t[0] == "B": cur["b"] = np.array(t[1:], float)
        elif t[0] == "C": i, j = int(t[1]), int(t[2]); cur["T"][i, j] = float(t[3]); cur["Ux"][i, j] = float(t[4]); cur["Uy"][i, j] = float(t[5])
        elif t[0] == "I": cur["I"].append([float(x) for x in t[2:]])
        elif t[0] == "O": cur["O"].append([float(x) for x in t[2:]])
        elif t[0] == "F": cur["F"].append([float(x) for x in t[1:]])
    for s in steps:
        for k in "IOF": s[k] = np.array(s[k])
    return steps
def budgets(steps, P, nx, Ltot, Lu, Li, rho_i, T_floor_inert=0.0):
    """returns per-step arrays: independent energy residual (physical: sensible + rho_i L swept volume vs boundary heat), mass residuals (global, local interface), all with their normalising scales."""
    H, dz, ny = P["H"], P["dz"], P["ny"]; Dx = Ltot / nx; r = rho_i / RHO_W; rc = RHO_W * CP_W; out = []
    icols = [i for i in range(nx) if Lu - 1e-12 <= (i + 0.5) * Dx <= Lu + Li + 1e-12]
    def Esens(s):
        h = (H - 0.5 * (s["b"][:-1] + s["b"][1:])) / ny; return rc * dz * Dx * np.sum(s["T"] * h[:, None])
    for a, b in zip(steps[:-1], steps[1:]):
        dt = b["t"] - a["t"]; Vsw = dz * Dx * np.sum(0.5 * ((a["b"][:-1] - b["b"][:-1]) + (a["b"][1:] - b["b"][1:])))      # swept volume (floor moves down)
        # inlet / outlet faces (rows j): area dz (H-b)/ny ; advective + diffusive heat (diffusive: first-order, face to cell centre Dx/2)
        Ai = dz * (H - b["b"][0]) / ny; Ao = dz * (H - b["b"][-1]) / ny
        Ui, Ti = b["I"][:, 0], b["I"][:, 2]; Uo, To = b["O"][:, 0], b["O"][:, 2]
        Qin = np.sum(Ui) * Ai; Qout = np.sum(Uo) * Ao; ein = RHO_W * CP_W * np.sum(Ui * Ti) * Ai
        ein_d = K_W * np.sum((Ti - b["T"][0, :]) / (0.5 * Dx)) * Ai      # heat INTO the domain through the inlet face: k (T_f - T_P)/(Dx/2)
        eout = RHO_W * CP_W * np.sum(Uo * To) * Ao; eout_d = -K_W * np.sum((To - b["T"][-1, :]) / (0.5 * Dx)) * Ao
        # inert isothermal floor conduction (first row cells, vertical half-height): heat into domain k (T_f - T_P) / (h/2)
        efl = 0.0
        for i in range(nx):
            if i in icols: continue
            hh = (H - 0.5 * (b["b"][i] + b["b"][i + 1])) / ny; efl += K_W * (T_floor_inert - b["T"][i, 0]) / (0.5 * hh) * Dx * dz
        heat_in = (ein + ein_d - eout - eout_d + efl) * dt; dEs = Esens(b) - Esens(a); Elat = rho_i * LF * Vsw
        scale = abs(ein) * dt + abs(eout) * dt + abs(efl) * dt + abs(dEs) + Elat + 1e-300
        # local interface mass: phi_b = dz (Ux db - Uy Dx)  vs (1-r) V_f / dt, V_f from consecutive geometry
        F = b["F"]; col = F[:, 0].astype(int); db = b["b"][col + 1] - b["b"][col]; phib = dz * (F[:, 1] * db - F[:, 2] * Dx)
        Vf = dz * Dx * 0.5 * ((a["b"][col] - b["b"][col]) + (a["b"][col + 1] - b["b"][col + 1])); target = (1 - r) * Vf / dt if r != 1.0 else np.zeros_like(Vf)
        mloc = np.max(np.abs(phib - target)) / max(np.max(np.abs(Vf / dt)), 1e-300)
        # global: Qout - Qin = -(net floor outflow) -> residual of continuity  Qin - Qout - sum phi_b (floor) - (inert junction faces ~ included in F only for ice patch)
        mglob = (Qin - Qout - np.sum(phib)) / max(abs(Qin), 1e-300)
        out.append(dict(t=b["t"], dt=dt, R_phys=dEs + Elat - heat_in, scale=scale, Vsw=Vsw, mloc=mloc, mglob=mglob, Qin=Qin, Qout=Qout, dEs=dEs, Elat=Elat, heat_in=heat_in))
    return out
