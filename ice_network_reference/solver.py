"""Independent spatially resolved reference model: unsteady 1-D channel flow (full St Venant equations, convective terms included) with
transported water temperature, wall heat exchange and continuous melt-driven geometry change, coupled to well-mixed lakes.

INDEPENDENCE from the reduced model (ice_network): nothing is imported from it. Different numerics and approximations:
  * hydraulics: explicit finite volume, HLL Riemann flux with hydrostatic reconstruction (well balanced for variable bed AND width),
    semi-implicit Manning friction, CFL-limited time step; wetting/drying and sill overtopping emerge from the reconstruction;
    inertia, waves and channel water storage are RESOLVED (the reduced model assumes quasi-steady edges and neglects channel storage);
  * thermal: first-order upwind transport of temperature on the same fine grid, then exact exponential wall exchange per cell;
  * geometry: bed and width of every cell change at every time step (no macro-step operator splitting, no station averaging);
  * melt water is added to the flowing water in the cell where it is produced (the reduced model adds it at the downstream node).
SHARED, DECLARED CLOSURES (identical numbers on purpose, so a disagreement cannot come from a different closure): Manning roughness n,
entrance-loss coefficient k_e, critical depth at a free outfall, Dittus-Boelter heat-transfer coefficient (times nu_factor), equal heat flux
over the wetted perimeter with the wall share spread over the full bank height, temperate ice, constants below. These are modelling
assumptions of both codes; this reference does NOT test them.
Lakes are 0-D (shared assumption): prismatic or power-law basins, well mixed.
"""
import numpy as np
from numba import njit

G = 9.81; RHO_W = 1000.0; RHO_I = 917.0; CP_W = 4217.0; LF = 3.34e5
K_W = 0.561; NU_W = 1.79e-6; PR_W = 13.4
HMIN = 1e-7


@njit(cache=True)
def _hll(hl, ul, hr, ur):
    """HLL flux of the depth-based shallow-water system, (mass, momentum) per unit width."""
    if hl <= 0.0 and hr <= 0.0:
        return 0.0, 0.0
    cl = np.sqrt(G * hl); cr = np.sqrt(G * hr)
    if hl > 0.0 and hr > 0.0:
        cs = 0.5 * (cl + cr) + 0.25 * (ul - ur)              # Einfeldt/Roe-type intermediate wave speed estimate
        us = 0.5 * (ul + ur) + cl - cr
        sl = min(ul - cl, us - cs)
        sr = max(ur + cr, us + cs)
    elif hl > 0.0:
        sl = ul - cl; sr = ul + 2.0 * cl
    else:
        sl = ur - 2.0 * cr; sr = ur + cr
    fl0 = hl * ul; fl1 = hl * ul * ul + 0.5 * G * hl * hl
    fr0 = hr * ur; fr1 = hr * ur * ur + 0.5 * G * hr * hr
    if sl >= 0.0:
        return fl0, fl1
    if sr <= 0.0:
        return fr0, fr1
    d = sr - sl
    return ((sr * fl0 - sl * fr0 + sl * sr * (hr - hl)) / d,
            (sr * fl1 - sl * fr1 + sl * sr * (hr * ur - hl * ul)) / d)


@njit(cache=True)
def htc(q_abs, w, h, nu_factor):
    a = w * h; p = w + 2.0 * h; dh = 4.0 * a / p; v = q_abs / a
    re = max(v * dh / NU_W, 1.0)
    return nu_factor * 0.023 * re ** 0.8 * PR_W ** 0.3 * K_W / dh


@njit(cache=True)
def channel_step(h, Q, w, z, zs, T, dx, dt, n_man, ke, h_lake, t_lake, h_pool, nu_factor, melting, heat_diss, out):
    """Advance one channel by dt. Arrays are updated in place. out (len 8): [V_in_face_m3, V_out_face_m3, sensible_in_J, heat_to_ice_J, sensible_out_J,
    melt_water_m3, dissipation_J, ice_volume_melted_m3] increments over this step."""
    n = h.size
    u = np.zeros(n)
    for i in range(n):
        if h[i] > HMIN:
            u[i] = Q[i] / (w[i] * h[i])
    eta = z + h
    fm = np.zeros(n + 1); fq_r = np.zeros(n); fq_l = np.zeros(n)
    # face loop: face k between cell k-1 and k; k = 0 and n are boundaries
    for k in range(n + 1):
        if k == 0:
            # Subcritical reservoir inlet: Bernoulli with entrance loss, g h_f + (1 + ke) u_f^2/2 = g (H_lake - z_0), closed by the
            # Riemann invariant leaving the domain, u_f - 2 c_f = u_0 - 2 c_0 (characteristic from the first cell). Energy-conserving
            # apart from the explicit loss (a reservoir-at-rest Riemann problem would instead conserve u + 2c, not energy).
            hres = h_lake - z[0]
            wf = w[0]
            if hres > 1e-9:
                c0 = np.sqrt(G * max(h[0], 0.0)); inv = u[0] - 2.0 * c0
                lo, hi = 0.0, hres
                for _ in range(60):
                    hm = 0.5 * (lo + hi); uf = inv + 2.0 * np.sqrt(G * hm)
                    f = G * hm + (1.0 + ke) * uf * abs(uf) / 2.0 - G * hres
                    if f > 0.0:
                        hi = hm
                    else:
                        lo = hm
                hf = 0.5 * (lo + hi); uf = inv + 2.0 * np.sqrt(G * hf)
                fm[0] = wf * hf * uf
                fq_l[0] = wf * (hf * uf * uf + 0.5 * G * hf * hf) + 0.5 * G * (w[0] * h[0] * h[0] - wf * hf * hf) * 0.0
            else:
                fm[0] = 0.0
                fq_l[0] = 0.5 * G * w[0] * h[0] * h[0]                       # closed wall
        elif k == n:
            q_last = Q[n - 1]
            hc = (q_last * q_last / (G * w[n - 1] * w[n - 1])) ** (1.0 / 3.0)
            hg = max(h_pool - z[n - 1], hc)                              # free outfall: critical depth; submerged: tailwater
            ug = q_last / (w[n - 1] * hg) if hg > HMIN else 0.0
            zf = z[n - 1]; wf = w[n - 1]
            hl_s = max(eta[n - 1] - zf, 0.0)
            f0, f1 = _hll(hl_s, u[n - 1], hg, ug)
            fm[n] = wf * f0
            fq_r[n - 1] = wf * f1 + 0.5 * G * (w[n - 1] * h[n - 1] * h[n - 1] - wf * hl_s * hl_s)
        else:
            zf = max(z[k - 1], z[k]); wf = min(w[k - 1], w[k])
            hl_s = max(eta[k - 1] - zf, 0.0); hr_s = max(eta[k] - zf, 0.0)
            f0, f1 = _hll(hl_s, u[k - 1], hr_s, u[k])
            fm[k] = wf * f0
            fq_r[k - 1] = wf * f1 + 0.5 * G * (w[k - 1] * h[k - 1] * h[k - 1] - wf * hl_s * hl_s)
            fq_l[k] = wf * f1 + 0.5 * G * (w[k] * h[k] * h[k] - wf * hr_s * hr_s)
    # conservative update (areas and discharges)
    a_old = w * h
    a_new = a_old.copy(); q_new = Q.copy()
    for i in range(n):
        a_new[i] = a_old[i] - dt / dx * (fm[i + 1] - fm[i])
        q_new[i] = Q[i] - dt / dx * (fq_r[i] - fq_l[i])
        if a_new[i] < 0.0:
            a_new[i] = 0.0
    # semi-implicit Manning friction, dissipation, thermal transport and wall exchange
    rc = RHO_W * CP_W
    t_new = T.copy()
    heat_step = 0.0; diss_step = 0.0; melt_vol = 0.0
    for i in range(n):
        hn = a_new[i] / w[i]
        if hn > HMIN:
            p = w[i] + 2.0 * hn; r = a_new[i] / p
            fric = G * dt * n_man * n_man * abs(q_new[i]) / (a_new[i] * r ** (4.0 / 3.0))
            qf = q_new[i] / (1.0 + fric)
            sf = n_man * n_man * qf * abs(qf) / (a_new[i] * a_new[i] * r ** (4.0 / 3.0))
            diss = RHO_W * G * abs(qf) * abs(sf) * dx              # W: work of the friction force in the cell
            q_new[i] = qf
        else:
            q_new[i] = 0.0; diss = 0.0
        # upwind advection of temperature (conservative in A*T)
        fin = fm[i]; fout = fm[i + 1]
        t_in = (t_lake if fin > 0.0 else T[i]) if i == 0 else (T[i - 1] if fin > 0.0 else T[i])
        t_out = (T[i] if fout > 0.0 else T[i + 1]) if i < n - 1 else T[i]
        if hn > HMIN:
            at = a_old[i] * T[i] - dt / dx * (fout * t_out - fin * t_in)
            t_adv = at / a_new[i] if a_new[i] > 0 else 0.0
        else:
            t_adv = T[i]
        # entrance-loss heating goes into the first cell
        phi_cell = diss
        if i == 0:
            phi_cell += RHO_W * abs(fm[0]) * ke * u[0] * u[0] / 2.0          # entrance-loss dissipation (W), heated into the first cell
        if not heat_diss:
            phi_cell = 0.0
        if hn > HMIN and melting:
            qa = abs(q_new[i])
            if qa > 1e-14:
                h_t = htc(qa, w[i], hn, nu_factor)
                p = w[i] + 2.0 * hn
                lam = h_t * p / (rc * a_new[i])
                s = phi_cell / (rc * a_new[i] * dx)
                e = np.exp(-lam * dt)
                t_after = t_adv * e + (s / lam) * (-np.expm1(-lam * dt)) if lam > 0 else t_adv + s * dt
                q_ice_j = rc * a_new[i] * dx * (t_adv - t_after) + phi_cell * dt      # J delivered to ice in this step
                t_new[i] = t_after
                dv_ice = q_ice_j / (RHO_I * LF)
                wall = 2.0 * hn / p
                dz = dv_ice * (1.0 - wall) / (w[i] * dx)
                dw = dv_ice * wall / ((zs[i] - z[i] + dz) * dx)
                z[i] -= dz; w[i] += dw
                # water: melt water added where produced (kg conserved: rho_i * dv_ice = rho_w * dv_w)
                a_pre = a_new[i]
                a_new[i] += RHO_I / RHO_W * dv_ice / dx
                t_new[i] = t_after * a_pre / a_new[i]                   # melt water enters at the melting point: dilutes, adds no sensible heat
                heat_step += q_ice_j; melt_vol += RHO_I / RHO_W * dv_ice
                diss_step += phi_cell * dt
                continue
        t_new[i] = t_adv + (phi_cell * dt / (rc * a_new[i] * dx) if hn > HMIN and a_new[i] > 0 else 0.0)
        diss_step += phi_cell * dt
    out[0] = fm[0] * dt; out[1] = fm[n] * dt
    out[2] = rc * fm[0] * dt * (t_lake if fm[0] > 0 else T[0])            # sensible heat entering from the lake (negative for back flow)
    out[3] = heat_step; out[4] = rc * fm[n] * dt * T[n - 1]                  # sensible heat leaving through the exit face
    out[5] = melt_vol; out[6] = diss_step; out[7] = heat_step / (RHO_I * LF)
    for i in range(n):
        h[i] = a_new[i] / w[i]; Q[i] = q_new[i]; T[i] = t_new[i]


def lake_stage(v, z_bottom, area, h_ref=1.0, q=0.0):
    return z_bottom + (((q + 1.0) * h_ref ** q * v / area) ** (1.0 / (q + 1.0)) if v > 0 else 0.0)


class RefNetwork:
    """One or more 0-D lakes and 1-D channels. lakes: dict(z_bottom, area, [q,h_ref]); channels: dict(lake, z (cells), w, zs, length, n_man, ke, pool)."""

    def __init__(self, lakes, channels, melting=True, heat_diss=True, nu_factor=1.0):
        self.lakes = lakes; self.melting = melting; self.heat_diss = heat_diss; self.nu = nu_factor
        self.ch = []
        for c in channels:
            n = len(c["z"]); d = dict(c)
            d["z"] = np.array(c["z"], float); d["w"] = np.array(c["w"], float); d["zs"] = np.array(c["zs"], float)
            d["dx"] = c["length"] / n; d["h"] = np.zeros(n); d["Q"] = np.zeros(n); d["T"] = np.zeros(n)
            d["z0"] = d["z"].copy(); d["w0"] = d["w"].copy()
            self.ch.append(d)

    def run(self, v0, t0_lake, t_end, inflow, inflow_temp, cfl=0.4, dt_max=5.0, record_every=60.0, h_pool_default=-1e9):
        v = np.array(v0, float); sv = v * np.array(t0_lake, float)
        t = 0.0; out = np.zeros(8); rec = []
        tot = dict(input=0.0, delivered=np.zeros(len(self.ch)), heat=0.0, melt_w=0.0, diss=0.0, sens_in_ext=0.0, sens_out=0.0, sens_from_lake=np.zeros(len(self.ch)))
        rc = RHO_W * CP_W
        next_rec = 0.0
        nsteps = 0
        while t < t_end - 1e-12:
            st = np.array([lake_stage(v[i], L["z_bottom"], L["area"], L.get("h_ref", 1.0), L.get("q", 0.0)) for i, L in enumerate(self.lakes)])
            tl = np.where(v > 1e-9, sv / np.maximum(v, 1e-9), 0.0)
            if t >= next_rec - 1e-12:
                rec.append(self._record(t, v, st, tl, tot, sv)); next_rec += record_every
            dt = dt_max
            for c in self.ch:
                hh = c["h"]; wet = hh > 1e-6
                u = np.abs(c["Q"] / (c["w"] * np.maximum(hh, 1e-6)))
                cmax = np.max(np.where(wet, u + np.sqrt(G * hh), 0.0))
                # a dry channel next to a reservoir still needs a step bound: use the inlet wave speed
                hin = max(st[c["lake"]] - c["z"][0], 0.0)
                cmax = max(cmax, np.sqrt(G * hin))
                if cmax > 0:
                    dt = min(dt, cfl * c["dx"] / cmax)
            dt = min(dt, t_end - t, next_rec - t if next_rec > t else dt)
            dv = np.zeros(len(v)); ds = np.zeros(len(v))
            I = inflow(t); TI = inflow_temp(t)
            dv += I * dt; ds += I * TI * dt
            tot["input"] += I.sum() * dt; tot["sens_in_ext"] += rc * (I * TI).sum() * dt
            for k, c in enumerate(self.ch):
                li = c["lake"]
                channel_step(c["h"], c["Q"], c["w"], c["z"], c["zs"], c["T"], c["dx"], dt, c["n_man"], c["ke"], st[li], tl[li],
                             c.get("pool", h_pool_default), self.nu, self.melting, self.heat_diss, out)
                dv[li] -= out[0]; ds[li] -= out[2] / rc
                tot["delivered"][k] += out[1]; tot["heat"] += out[3]; tot["sens_out"] += out[4]; tot["melt_w"] += out[5]; tot["diss"] += out[6]
            v = v + dv; sv = sv + ds
            t += dt; nsteps += 1
        st = np.array([lake_stage(v[i], L["z_bottom"], L["area"], L.get("h_ref", 1.0), L.get("q", 0.0)) for i, L in enumerate(self.lakes)])
        tl = np.where(v > 1e-9, sv / np.maximum(v, 1e-9), 0.0)
        rec.append(self._record(t, v, st, tl, tot, sv))
        self.nsteps = nsteps
        return dict(rec=rec, totals=tot, V=v, S=sv, nsteps=nsteps)

    def _record(self, t, v, st, tl, tot, sv_):
        return dict(t=t, V=v.copy(), stage=st.copy(), T=tl.copy(), Q_in=np.array([c["Q"][0] for c in self.ch]), Q_out=np.array([c["Q"][-1] for c in self.ch]),
                    delivered=tot["delivered"].copy(), sill=np.array([c["z"][0] for c in self.ch]), heat=tot["heat"], input=tot["input"],
                    melt_w=tot["melt_w"], diss=tot["diss"], sens_out=tot["sens_out"], sens_in_ext=tot["sens_in_ext"],
                    channel_volume=np.array([float(np.sum(c["w"] * c["h"]) * c["dx"]) for c in self.ch]),
                    chan_sens=sum(RHO_W * CP_W * float(np.sum(c["w"] * c["h"] * c["T"]) * c["dx"]) for c in self.ch),
                    lake_sens=RHO_W * CP_W * float(np.sum(sv_)),
                    ice_geom=sum(float(np.sum(((c["zs"] - c["z"]) * c["w"] - (c["zs"] - c["z0"]) * c["w0"])) * c["dx"]) for c in self.ch))
