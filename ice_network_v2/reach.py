"""Prescribed-discharge reach step (v2 baseline). DERIVED FROM ice_network_reference/solver.py::channel_step (same HLL hydrostatic-reconstruction hydraulics, semi-implicit Manning friction,
first-order upwind temperature, exact exponential wall exchange with Dittus-Boelter h_t, melt water added in the producing cell, one q_ice used for BOTH water cooling and ice melt).
Changes: (1) inlet face: discharge q_in is FORCING (subcritical characteristic closure for depth, no entrance loss); (2) inlet temperature is an external flux-weighted bulk temperature;
(3) optional fixed-width mode (evolve_width=False): all interface energy lowers the bed, width is unchanged and NOT an output. The frozen reference file is not modified."""
import numpy as np
from numba import njit
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from ice_network_reference.solver import _hll, htc, G, RHO_W, RHO_I, CP_W, LF, HMIN, NU_W


@njit(cache=True)
def channel_step_q(h, Q, w, z, zs, T, dx, dt, n_man, ke, q_in, t_in_ext, h_pool, nu_factor, melting, heat_diss, evolve_width, out):
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
            # PRESCRIBED-DISCHARGE inlet (v2): Q_in is forcing. Depth from the characteristic leaving the domain, u_f - 2 c_f = u_0 - 2 c_0, with w h_f u_f = Q_in.
            wf = w[0]
            c0 = np.sqrt(G * max(h[0], 0.0)); inv = u[0] - 2.0 * c0
            lo, hi = 1e-9, 50.0
            for _ in range(80):
                hm = 0.5 * (lo + hi); uf = q_in / (wf * hm); f = uf - 2.0 * np.sqrt(G * hm) - inv      # decreasing in hm
                if f > 0.0:
                    lo = hm
                else:
                    hi = hm
            hf = 0.5 * (lo + hi); uf = q_in / (wf * hf)
            fm[0] = q_in
            fq_l[0] = wf * (hf * uf * uf + 0.5 * G * hf * hf)
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
        t_in = (t_in_ext if fin > 0.0 else T[i]) if i == 0 else (T[i - 1] if fin > 0.0 else T[i])
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
                if evolve_width:
                    dz = dv_ice * (1.0 - wall) / (w[i] * dx)
                    dw = dv_ice * wall / ((zs[i] - z[i] + dz) * dx)
                else:
                    dz = dv_ice / (w[i] * dx); dw = 0.0     # fixed-width mode: the whole interface energy lowers the bed (declared approximation; walls not eroded)
                z[i] -= dz
                if evolve_width: w[i] += dw
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
    out[2] = rc * fm[0] * dt * (t_in_ext if fm[0] > 0 else T[0])            # sensible heat entering from the lake (negative for back flow)
    out[3] = heat_step; out[4] = rc * fm[n] * dt * T[n - 1]                  # sensible heat leaving through the exit face
    out[5] = melt_vol; out[6] = diss_step; out[7] = heat_step / (RHO_I * LF)
    for i in range(n):
        h[i] = a_new[i] / w[i]; Q[i] = q_new[i]; T[i] = t_new[i]


