"""Plan-view depth-averaged free-surface solver with advected sensible heat and Stefan-type ice-floor lowering (M3). NOT 3-D VOF CFD.
Why this and not the OpenFOAM stack: the repository's free-surface solver (interThermalFoam) has no phase change/ice-motion coupling, and the verified flow-heat-moving-ice solver (flowStefanFoam) is single-phase rigid-lid. A verified VOF+Stefan+moving-mesh solver does not exist here; this depth-averaged model is the supported alternative.
Physics: shallow-water equations (hydrostatic reconstruction -> exactly well balanced, wetting/drying, HLL), Manning friction, optional Ergun porous drag in debris cells, conservative depth-weighted temperature, bed-exchange heat sink h_eff*T on wet ice cells that is the SAME flux used to lower the ice floor (rho_i L dz/dt = h_eff T, quasi-static ice).
Hydraulics are advanced to quasi-steady between morphological (bed) updates; the water temperature field too. Bed update uses the time-mean sink accumulated in the last window.
Idealizations (all labelled): solid debris = raised bed (non-porous limit); debris at T_m so its export carries no enthalpy relative to T_m; debris cells do not receive water heat when dry; no seepage-heat; steady-flow approximation between bed updates (error controlled by the bed-update interval series)."""
import math
import numpy as np
from numba import njit
G = 9.81; RHO_W = 1000.0; CP = 4217.0; RHO_I = 917.0; LF = 3.34e5; HDRY = 1e-4

@njit(cache=True)
def _hll(hL, unL, utL, hR, unR, utR, g):
    if hL <= 0.0 and hR <= 0.0: return 0.0, 0.0, 0.0
    cL = math.sqrt(g * hL); cR = math.sqrt(g * hR)
    if hL <= 0.0: sL = unR - 2.0 * cR
    else: sL = min(unL - cL, unR - cR) if hR > 0.0 else unL - cL
    if hR <= 0.0: sR = unL + 2.0 * cL
    else: sR = max(unL + cL, unR + cR) if hL > 0.0 else unR + cR
    FL0 = hL * unL; FL1 = hL * unL * unL + 0.5 * g * hL * hL; FL2 = hL * unL * utL
    FR0 = hR * unR; FR1 = hR * unR * unR + 0.5 * g * hR * hR; FR2 = hR * unR * utR
    if sL >= 0.0: return FL0, FL1, FL2
    if sR <= 0.0: return FR0, FR1, FR2
    d = sR - sL
    return ((sR * FL0 - sL * FR0 + sL * sR * (hR - hL)) / d, (sR * FL1 - sL * FR1 + sL * sR * (hR * unR - hL * unL)) / d, (sR * FL2 - sL * FR2 + sL * sR * (hR * utR - hL * utL)) / d)

@njit(cache=True)
def _bc(h, hu, hv, hT, nx, ny, qin_w, Tin, closed_out):
    for j in range(ny + 2):
        h[0, j] = h[1, j]; hu[0, j] = qin_w; hv[0, j] = 0.0; hT[0, j] = h[1, j] * Tin       # inlet ghost: prescribed discharge per width, T_in
        h[nx + 1, j] = h[nx, j]; hu[nx + 1, j] = hu[nx, j]; hv[nx + 1, j] = hv[nx, j]; hT[nx + 1, j] = hT[nx, j]      # outlet: zero gradient (free outfall; verified supercritical by the outlet-length study)
        if closed_out: hu[nx + 1, j] = 0.0; hu[nx, j] = hu[nx, j]
    for i in range(nx + 2):
        h[i, 0] = h[i, 1]; hu[i, 0] = hu[i, 1]; hv[i, 0] = -hv[i, 1]; hT[i, 0] = hT[i, 1]
        h[i, ny + 1] = h[i, ny]; hu[i, ny + 1] = hu[i, ny]; hv[i, ny + 1] = -hv[i, ny]; hT[i, ny + 1] = hT[i, ny]

@njit(cache=True)
def advance(h, hu, hv, hT, zs, zdp, heff, iceMask, manning, qin_w, Tin, dx, nx, ny, tmax, cfl, porous_k, porous_phi, porous_D, nu, cool, qacc, flux_acc, closed_out):
    """advance for tmax seconds of hydraulic time (or less steps); returns (t, steps, max wave speed used). qacc accumulates sink heat flux*dt per cell [J/m2]; flux_acc[0..3] accumulates inlet water volume, outlet water volume, inlet heat, outlet heat [m3, m3, J, J]."""
    t = 0.0; steps = 0; smax_all = 0.0
    nh = np.zeros_like(h); nhu = np.zeros_like(h); nhv = np.zeros_like(h); nhT = np.zeros_like(h)
    while t < tmax:
        _bc(h, hu, hv, hT, nx, ny, qin_w, Tin, closed_out)
        smax = 1e-9
        for i in range(1, nx + 1):
            for j in range(1, ny + 1):
                if h[i, j] > HDRY:
                    c = math.sqrt(G * h[i, j]); smax = max(smax, (abs(hu[i, j]) + abs(hv[i, j])) / h[i, j] + 2 * c)
        dt = min(cfl * dx / smax, tmax - t)
        for i in range(nx + 2):
            for j in range(ny + 2):
                nh[i, j] = h[i, j]; nhu[i, j] = hu[i, j]; nhv[i, j] = hv[i, j]; nhT[i, j] = hT[i, j]
        lam = dt / dx
        for i in range(1, nx + 2):               # x faces between (i-1) and i
            for j in range(1, ny + 1):
                a = i - 1; b = i
                zf = max(zs[a, j], zs[b, j])
                hL = max(0.0, h[a, j] + zs[a, j] - zf); hR = max(0.0, h[b, j] + zs[b, j] - zf)
                uL = hu[a, j] / h[a, j] if h[a, j] > HDRY else 0.0; vL = hv[a, j] / h[a, j] if h[a, j] > HDRY else 0.0
                uR = hu[b, j] / h[b, j] if h[b, j] > HDRY else 0.0; vR = hv[b, j] / h[b, j] if h[b, j] > HDRY else 0.0
                f0, f1, f2 = _hll(hL, uL, vL, hR, uR, vR, G)
                TL = hT[a, j] / h[a, j] if h[a, j] > HDRY else 0.0; TR = hT[b, j] / h[b, j] if h[b, j] > HDRY else 0.0
                fT = f0 * (TL if f0 > 0 else TR)
                if a >= 1:
                    nh[a, j] -= lam * f0; nhu[a, j] -= lam * (f1 + 0.5 * G * (h[a, j] ** 2 - hL * hL)); nhv[a, j] -= lam * f2; nhT[a, j] -= lam * fT
                if b <= nx:
                    nh[b, j] += lam * f0; nhu[b, j] += lam * (f1 + 0.5 * G * (h[b, j] ** 2 - hR * hR)); nhv[b, j] += lam * f2; nhT[b, j] += lam * fT
                if i == 1: flux_acc[0] += f0 * dx * dt; flux_acc[2] += RHO_W * CP * fT * dx * dt
                if i == nx + 1: flux_acc[1] += f0 * dx * dt; flux_acc[3] += RHO_W * CP * fT * dx * dt
        for i in range(1, nx + 1):               # y faces between (j-1) and j
            for j in range(1, ny + 2):
                a = j - 1; b = j
                zf = max(zs[i, a], zs[i, b])
                hL = max(0.0, h[i, a] + zs[i, a] - zf); hR = max(0.0, h[i, b] + zs[i, b] - zf)
                vL = hv[i, a] / h[i, a] if h[i, a] > HDRY else 0.0; uL = hu[i, a] / h[i, a] if h[i, a] > HDRY else 0.0
                vR = hv[i, b] / h[i, b] if h[i, b] > HDRY else 0.0; uR = hu[i, b] / h[i, b] if h[i, b] > HDRY else 0.0
                f0, f1, f2 = _hll(hL, vL, uL, hR, vR, uR, G)       # f1: normal (y) momentum flux, f2: tangential (x)
                TL = hT[i, a] / h[i, a] if h[i, a] > HDRY else 0.0; TR = hT[i, b] / h[i, b] if h[i, b] > HDRY else 0.0
                fT = f0 * (TL if f0 > 0 else TR)
                if a >= 1:
                    nh[i, a] -= lam * f0; nhv[i, a] -= lam * (f1 + 0.5 * G * (h[i, a] ** 2 - hL * hL)); nhu[i, a] -= lam * f2; nhT[i, a] -= lam * fT
                if b <= ny:
                    nh[i, b] += lam * f0; nhv[i, b] += lam * (f1 + 0.5 * G * (h[i, b] ** 2 - hR * hR)); nhu[i, b] += lam * f2; nhT[i, b] += lam * fT
        for i in range(1, nx + 1):
            for j in range(1, ny + 1):
                hn = nh[i, j]
                if hn <= HDRY:
                    flux_acc[4] += (max(hn, 0.0) - hn) * dx * dx; h[i, j] = max(hn, 0.0); hu[i, j] = 0.0; hv[i, j] = 0.0; hT[i, j] = 0.0 if hn <= 0 else nhT[i, j]
                    continue
                ux = nhu[i, j] / hn; uy = nhv[i, j] / hn; sp = math.sqrt(ux * ux + uy * uy)
                S = G * manning * manning * sp / hn ** (4.0 / 3.0)
                if zdp[i, j] > 0.0 and porous_k > 0:       # Ergun drag within the debris layer fraction
                    frac = min(hn, zdp[i, j]) / hn; Sp = 150.0 * nu * (1 - porous_phi) ** 2 / (porous_phi ** 3 * porous_D ** 2) + 1.75 * (1 - porous_phi) * sp / (porous_phi ** 3 * porous_D)
                    S += Sp * frac
                fac = 1.0 / (1.0 + dt * S)
                h[i, j] = hn; hu[i, j] = nhu[i, j] * fac; hv[i, j] = nhv[i, j] * fac
                Tn = nhT[i, j] / hn
                if iceMask[i, j] and cool:
                    q = heff[i, j] * max(Tn, 0.0); dT = dt * q / (RHO_W * CP * hn)
                    if dT > Tn: dT = Tn; q = dT * RHO_W * CP * hn / dt
                    Tn -= dT; qacc[i, j] += q * dt
                hT[i, j] = Tn * hn
        t += dt; steps += 1; smax_all = max(smax_all, smax)
    return t, steps, smax_all

class Sim:
    """cfg keys: dx, Lpool, Lr, Ltail, Wl (lane width), Q, Tw, z_toe, z_byp, deposit: array spec via set_deposit(), heff, k_deb, manning."""
    def __init__(self, dx=0.05, Lpool=2.0, Lr=2.0, Ltail=1.0, Wl=0.5, Q=0.07, Tw=0.5, z_toe=-0.08, z_byp=0.0, heff=200.0, k_deb=2.0, manning=0.02, cfl=0.4, dz=1.0, porous=None, closed_out=False, wall_h=5.0, tail_slope=0.5, ramp=0.0, pool_floor=-0.7, tail_drop=1.5, ramp_ratio=3.0):
        self.dx = dx; self.nx = int(round((Lpool + Lr + Ltail) / dx)); self.nw = int(round(Wl / dx)); self.ny = 2 * self.nw + 1
        self.Lpool, self.Lr, self.Wl, self.Q, self.Tw, self.heff0, self.k_deb, self.manning, self.cfl = Lpool, Lr, Wl, Q, Tw, heff, k_deb, manning, cfl
        nx, ny = self.nx, self.ny; self.i0 = int(round(Lpool / dx)) + 1; self.i1 = int(round((Lpool + Lr) / dx))      # route cells i0..i1 (1-based padded)
        self.zb = np.zeros((nx + 2, ny + 2)); self.zb[:self.i0] = pool_floor      # deep inert reservoir: NO fixed sill (the entrance control is each lane's own ice floor, so incision can lower the control)
        self.zd = np.zeros((nx + 2, ny + 2)); self.wall = np.zeros((nx + 2, ny + 2), bool); self.ice = np.zeros((nx + 2, ny + 2), bool)
        self.toe = np.zeros((nx + 2, ny + 2), bool); self.byp = np.zeros((nx + 2, ny + 2), bool)
        tl = slice(1, self.nw + 1); bl = slice(self.nw + 2, ny + 1)         # lane y-index ranges (toe low y, bypass high y), wall row nw+1
        for i in range(1, nx + 1):
            if i >= self.i0:      # routes and tail: toe lane floor z_toe, bypass floor z_byp, divider wall
                self.zb[i, tl] = z_toe; self.zb[i, bl] = z_byp; self.wall[i, self.nw + 1] = True
                if i > self.i1: dd = tail_drop + tail_slope * (i - self.i1 - 0.5) * dx; self.zb[i, tl] -= dd; self.zb[i, bl] -= dd      # free overfall (drop 1.5 m, no fixed downstream sill so incision can lower the control) then slope 0.5: critical control at the brink, supercritical tail -> zero-gradient outlet admissible (checked by outlet-length study)
                if i <= self.i1: self.toe[i, tl] = True; self.byp[i, bl] = True; self.ice[i, tl] = True; self.ice[i, bl] = True
        self.ramp_ratio = ramp_ratio; self.pool_floor = pool_floor
        if ramp_ratio:      # entrance ramps (run:rise = ramp_ratio) from the lane floors down to the deep reservoir; they are ICE (they incise with the lanes) so no fixed sill can form; abrupt steps were found to dissipate energy artificially (diagnose_v1.json)
            for i in range(1, self.i0):
                x = (i - 0.5) * dx
                for rows, zl, nm in ((tl, z_toe, "toe"), (bl, z_byp, "byp")):
                    z = max(pool_floor, zl - (Lpool - x) / ramp_ratio); self.zb[i, rows] = z
                    if z > pool_floor + 1e-9:
                        self.ice[i, rows] = True; (self.toe if nm == "toe" else self.byp)[i, rows] = True
                self.zb[i, self.nw + 1] = 0.5 * (self.zb[i, 1] + self.zb[i, ny])
        self.zb[:, 0] = self.zb[:, 1]; self.zb[:, ny + 1] = self.zb[:, ny]; self.zb[0] = self.zb[1]; self.zb[nx + 1] = self.zb[nx]
        self.porous = porous; self.closed_out = closed_out; self.wall_h = wall_h; self.t = 0.0
        self.h = np.zeros((nx + 2, ny + 2)); self.hu = np.zeros_like(self.h); self.hv = np.zeros_like(self.h); self.hT = np.zeros_like(self.h)
        self.solid = np.zeros((nx + 2, ny + 2), bool); self.pdep = np.zeros((nx + 2, ny + 2))     # solid-raised-bed debris height / porous debris height
        self.removed = False; self.export_log = []
        self.reset_state(0.15)
    def set_deposit(self, height, mode="solid"):
        """height: array (nx, ny) of deposit height [m] (0 = none); solid raised bed or porous drag layer."""
        H = np.zeros_like(self.zb); H[1:-1, 1:-1] = height; m = H > 0
        if mode == "solid": self.solid[m] = True; self.zd[m] = H[m]
        else: self.pdep[m] = H[m]
        self.deposit_volume = float(H.sum() * self.dx ** 2); self.deposit_footprint = float(m.sum() * self.dx ** 2)
    def zs(self):
        z = self.zb + np.where(self.solid, self.zd, 0.0); z = np.where(self.wall, self.wall_h, z); z[:, 0] = z[:, 1]; z[:, -1] = z[:, -2]; z[0] = z[1]; z[-1] = z[-2]; return z      # ghost layers copy the interior (otherwise deposits at the lateral boundary leak water)
    def reset_state(self, eta):
        zs = self.zs(); self.h = np.maximum(eta - zs, 0.0) * (~self.wall); self.hu[:] = 0; self.hv[:] = 0; self.hT = self.h * self.Tw
    def heff_field(self):
        he = np.full_like(self.zb, self.heff0); d = np.where(self.solid, self.zd, 0.0) + self.pdep
        cov = d > 0; he[cov] = 1.0 / (1.0 / self.heff0 + d[cov] / self.k_deb); return he
    def hydraulics(self, tmax=120.0, tol=2e-4, window=4.0, min_t=8.0, cool=True):
        he = self.heff_field(); zs = self.zs(); qacc = np.zeros_like(self.zb); fa = np.zeros(5); tt = 0.0; steps = 0; last = None; smax = 0.0
        pk, pphi, pD = (self.porous or (0.0, 0.4, 0.2)) if self.porous is not None else (0.0, 0.4, 0.2); pk = 1.0 if self.porous is not None else 0.0
        w = np.zeros(4)
        while tt < tmax:
            qacc[:] = 0; fa[:] = 0; self._V0 = float(self.h[1:-1, 1:-1].sum() * self.dx ** 2); self._E0 = float(RHO_W * CP * self.hT[1:-1, 1:-1].sum() * self.dx ** 2)
            t, s, sm = advance(self.h, self.hu, self.hv, self.hT, zs, self.pdep, he, self.ice, self.manning, self.Q / (self.ny * self.dx), self.Tw, self.dx, self.nx, self.ny, window, self.cfl, pk, pphi, pD, 1e-6, cool, qacc, fa, self.closed_out)
            self._dV = (float(self.h[1:-1, 1:-1].sum() * self.dx ** 2) - self._V0) / t; self._dE = (float(RHO_W * CP * self.hT[1:-1, 1:-1].sum() * self.dx ** 2) - self._E0) / t; self._win = t
            tt += t; steps += s; smax = max(smax, sm)
            cur = np.array([self.Qsec(self.toe_cols()), self.Qsec(self.byp_cols()), fa[1] / t, fa[3] / max(fa[2], 1e-300)])
            if last is not None and tt >= min_t and abs(cur[0] - last[0]) < tol * self.Q and abs(cur[1] - last[1]) < tol * self.Q and abs(cur[3] - last[3]) < 1e-4 and abs(self._dV) < 1e-3 * self.Q:
                break
            last = cur
        self.last_qacc = qacc / window; self.last_fa = fa / window; self.hyd_time = tt; self.hyd_steps = steps; self.smax = smax
        return dict(converged=tt < tmax, hyd_time=tt, steps=steps, dt_limit=self.cfl * self.dx / smax)
    def toe_cols(self): return slice(1, self.nw + 1)
    def byp_cols(self): return slice(self.nw + 2, self.ny + 1)
    def Qsec(self, ys, i=None):
        i = i or (self.i0 + (self.i1 - self.i0) // 2); return float((self.hu[i, ys]).sum() * self.dx)
    def diagnostics(self):
        d = {}; Qt, Qb = self.Qsec(self.toe_cols()), self.Qsec(self.byp_cols()); fa = self.last_fa
        d.update(Q_toe=Qt, Q_byp=Qb, f_toe=Qt / self.Q, Q_in=fa[0], Q_out=fa[1], water_resid=(fa[0] - fa[1] - self._dV) / self.Q, water_storage_rate=self._dV)
        wet = (self.h > 1e-3) & self.ice; d["wet_ice_area_toe"] = float((wet & self.toe).sum() * self.dx ** 2); d["wet_ice_area_byp"] = float((wet & self.byp).sum() * self.dx ** 2)
        P_in, P_out = fa[2], fa[3]; melt_rate = float((self.last_qacc * self.ice).sum() * self.dx ** 2)
        d.update(heat_in_W=P_in, heat_out_W=P_out, melt_heat_W=melt_rate, energy_resid_rel=(P_in - P_out - melt_rate - self._dE) / max(P_in, 1e-300), storage_rate_W=self._dE)
        d["melt_heat_toe_W"] = float((self.last_qacc * self.toe).sum() * self.dx ** 2); d["melt_heat_byp_W"] = float((self.last_qacc * self.byp).sum() * self.dx ** 2); return d
    def lower_bed(self, dtm):
        dz = self.last_qacc * dtm / (RHO_I * LF) * self.ice; self.zb = np.maximum(self.zb - dz, self.pool_floor); vol = float(dz.sum() * self.dx ** 2); return vol, float((dz * self.toe).sum() * self.dx ** 2), float((dz * self.byp).sum() * self.dx ** 2)
    def remove_debris(self):
        """prescribed complete removal; exported solid volume logged; debris assumed at T_m so exported enthalpy relative to T_m is zero (declared); water in the vacated cells keeps its own heat via the (continuing) flow solution."""
        v = float((np.where(self.solid, self.zd, 0.0) + self.pdep).sum() * self.dx ** 2); self.solid[:] = False; self.zd[:] = 0; self.pdep[:] = 0; self.removed = True; self.export_log.append(dict(t=self.t, volume_m3=v, enthalpy_J_rel_Tm=0.0)); return v
