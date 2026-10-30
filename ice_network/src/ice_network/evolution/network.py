"""Coupled network: node storage, quasi-steady edge hydraulics, thermal transport, melt-driven geometry update.

State vector (all extensive quantities, so the discrete budgets are linear invariants of the ODE):
  V[NL]            lake volumes (m3)                          S[NL]  = V * T  (m3 K), T relative to the melting point
  C[NO]            cumulative water delivered to each outlet  Eout[NO] cumulative sensible heat leaving (J)
  Win, Ein         cumulative external water (m3) and sensible heat (J)
  Qm               cumulative meltwater volume added (m3)     Phi  cumulative dissipation rho g sum |Q| (H_a - H_b) (J)
  Hseg[E*(n-1)]    cumulative heat delivered to the ice per segment (J), physical station order
  Wseg[E*(n-1)]    the part of it delivered to the wetted walls (J)
Geometry (z, w) is frozen inside a macro-step and updated afterwards (first-order operator split; refine dt_geom).
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

from ..constants import CP_W, G, LF, RHO_I, RHO_W
from ..geometry.storage import area_from_depth, depth_from_volume, volume_from_depth
from ..hydraulics.channel import edge_flow
from ..thermal.channel_heat import reach_heat


class Network:
    def __init__(self, lakes, outlets, edges, melting=True, nu_factor=1.0, wall_melt=True, heat_dissipation=True):
        """lakes: list of dict(z_bottom, a_ref, h_ref, q); outlets: list of dict(h_fixed);
        edges: list of dict(a, b, z, w, zs, n_man, ke) with a, b node ids: 0..NL-1 lakes, NL..NL+NO-1 outlets.
        z, w, zs are station arrays (physical orientation a -> b)."""
        self.lakes, self.outlets = lakes, outlets
        self.NL, self.NO, self.NE = len(lakes), len(outlets), len(edges)
        self.edges = [dict(e) for e in edges]
        for e in self.edges:
            for k in ("z", "w", "zs"):
                e[k] = np.array(e[k], float)
            e["n"] = e["z"].size
            e["dx"] = float(e["length"]) / (e["n"] - 1)
        self.melting, self.nu_factor, self.wall_melt = melting, nu_factor, wall_melt
        self.heat_dissipation = heat_dissipation   # False: dissipation is excluded from the thermal problem AND from the energy budget (isolates advection/mixing/melting)
        self.nseg = [e["n"] - 1 for e in self.edges]
        self.seg_off = np.concatenate([[0], np.cumsum(self.nseg)])
        self.NS = int(self.seg_off[-1])
        for e in self.edges:      # an outlet is a free-outfall sink: its stage must lie below the exit bed of its edges
            if e["b"] >= self.NL:
                assert self.outlets[e["b"] - self.NL]["h_fixed"] <= e["z"][-1] + 1e-9, "outlet stage above exit bed"
        self.reset_state()

    # ----------------------------------------------------------------------------- state layout
    def reset_state(self, v0=None, t0=None):
        NL, NO, NS = self.NL, self.NO, self.NS
        self.i_V = slice(0, NL)
        self.i_S = slice(NL, 2 * NL)
        self.i_C = slice(2 * NL, 2 * NL + NO)
        self.i_E = slice(2 * NL + NO, 2 * NL + 2 * NO)
        base = 2 * NL + 2 * NO
        self.i_Win, self.i_Ein, self.i_Qm, self.i_Phi, self.i_PhiC = base, base + 1, base + 2, base + 3, base + 4
        self.i_H = slice(base + 5, base + 5 + NS)
        self.i_W = slice(base + 5 + NS, base + 5 + 2 * NS)
        self.n_state = base + 5 + 2 * NS
        y = np.zeros(self.n_state)
        if v0 is not None:
            y[self.i_V] = v0
        if t0 is not None:
            y[self.i_S] = np.asarray(v0) * np.asarray(t0)
        return y

    # ----------------------------------------------------------------------------- levels
    def lake_stage(self, v):
        return np.array([L["z_bottom"] + depth_from_volume(v[i], L["a_ref"], L["h_ref"], L["q"]) for i, L in enumerate(self.lakes)])

    def stages(self, v):
        return np.concatenate([self.lake_stage(v), [o["h_fixed"] for o in self.outlets]])

    # ----------------------------------------------------------------------------- right-hand side
    def edge_state(self, e, h_a, h_b, t_a, t_b=0.0):
        """Hydraulics + thermal for one edge. Returns dict with flow-direction quantities."""
        n = e["n"]
        h = np.zeros(n); hf = np.zeros(n - 1); ch = np.zeros(n, np.int64)
        q, d, e_entry, e_exit = edge_flow(h_a, h_b, e["z"], e["w"], e["n_man"], e["dx"], e["ke"], int(e.get("msub", 16)), h, hf, ch)
        out = dict(Q=q, dir=d, h=h, hf=hf, choked=ch)
        if q == 0.0:
            out.update(Qabs=0.0, phi_in=0.0, phi_exit=0.0, q_ice=np.zeros(n - 1), wall_frac=np.zeros(n - 1), t_out=t_a, q_melt=0.0, phi_seg=np.zeros(n - 1))
            return out
        qa = abs(q)
        t_a = t_a if d > 0 else t_b          # inflow temperature is that of the SOURCE node of the flow
        wf = e["w"] if d > 0 else e["w"][::-1]
        v0 = qa / (wf[0] * h[0])
        phi_in = RHO_W * qa * e["ke"] * v0 * v0 / 2.0
        phi_seg = RHO_W * G * qa * hf
        drop = max(h_a, h_b) - min(h_a, h_b)
        phi_tot = RHO_W * G * qa * drop
        phi_exit = phi_tot - phi_in - phi_seg.sum()
        if not self.heat_dissipation:
            phi_in, phi_seg, phi_exit, phi_tot = 0.0, np.zeros_like(phi_seg), 0.0, 0.0
        out.update(Qabs=qa, phi_in=phi_in, phi_seg=phi_seg, phi_exit=phi_exit, phi_tot=phi_tot)
        t_arr = np.zeros(n); q_ice = np.zeros(n - 1)
        if self.melting:
            t_out = reach_heat(qa, t_a, 0.0, h, wf, e["dx"], phi_seg, phi_in, self.nu_factor, t_arr, q_ice)
            wm = 0.5 * (wf[:-1] + wf[1:]); hm = 0.5 * (h[:-1] + h[1:])
            wall_frac = (2.0 * hm / (wm + 2.0 * hm)) if self.wall_melt else np.zeros(n - 1)
            out["q_melt"] = q_ice.sum() / (RHO_W * LF)
        else:
            # no-melting control: dissipation is not allowed to melt ice; all heat stays in the water
            t_out = t_a + (phi_in + phi_seg.sum()) / (RHO_W * CP_W * qa)
            wall_frac = np.zeros(n - 1)
            out["q_melt"] = 0.0
            q_ice[:] = 0.0
        out.update(q_ice=q_ice, wall_frac=wall_frac, t_out=t_out)
        return out

    def rhs(self, t, y, inflow, inflow_temp):
        NL, NO = self.NL, self.NO
        v = np.maximum(y[self.i_V], 0.0)
        s = y[self.i_S]
        temps = np.where(v > 1e-9, s / np.maximum(v, 1e-9), 0.0)
        stage = self.stages(v)
        dy = np.zeros_like(y)
        dV, dS = dy[self.i_V], dy[self.i_S]
        dC, dE = dy[self.i_C], dy[self.i_E]
        I = inflow(t); TI = inflow_temp(t)
        dV += I; dS += I * TI
        dy[self.i_Win] = I.sum(); dy[self.i_Ein] = RHO_W * CP_W * (I * TI).sum()
        phi_total = 0.0; qm_total = 0.0; phi_chan = 0.0
        for ie, e in enumerate(self.edges):
            a, b = e["a"], e["b"]
            st = self.edge_state(e, stage[a], stage[b], temps[a] if a < NL else 0.0, temps[b] if b < NL else 0.0)
            # source node of the flow
            if st["Q"] == 0.0:
                continue
            src, dst = (a, b) if st["dir"] > 0 else (b, a)
            assert src < NL, "flow originates at an outlet (stage above exit bed)"
            qa = st["Qabs"]
            dV[src] -= qa; dS[src] -= qa * temps[src]
            gain = qa + st["q_melt"]
            heat_s = qa * st["t_out"] + st["phi_exit"] / (RHO_W * CP_W)
            if dst < NL:
                dV[dst] += gain; dS[dst] += heat_s
            else:
                dC[dst - NL] += gain; dE[dst - NL] += RHO_W * CP_W * qa * st["t_out"] + st["phi_exit"]
            phi_total += st["phi_tot"]; qm_total += st["q_melt"]; phi_chan += st["phi_in"] + st["phi_seg"].sum()
            ns = e["n"] - 1
            o = self.seg_off[ie]
            qi, wf_ = st["q_ice"], st["wall_frac"]
            if st["dir"] < 0:
                qi, wf_ = qi[::-1], wf_[::-1]      # flow-direction order -> physical station order
            dy[self.i_H][o:o + ns] += qi
            dy[self.i_W][o:o + ns] += qi * wf_
        dy[self.i_Qm] = qm_total; dy[self.i_Phi] = phi_total; dy[self.i_PhiC] = phi_chan
        return dy

    # ----------------------------------------------------------------------------- geometry update
    def update_geometry(self, d_hseg, d_wseg):
        """Apply melt to bed and walls from the heat delivered to each segment over the last macro-step.
        Returns (ice volume from energy, ice volume change implied by the new geometry) for the consistency audit."""
        v_energy = v_geom = 0.0
        for ie, e in enumerate(self.edges):
            ns, o = e["n"] - 1, self.seg_off[ie]
            dh, dw = d_hseg[o:o + ns], d_wseg[o:o + ns]
            dv_ice = dh / (RHO_I * LF)
            dv_wall = dw / (RHO_I * LF)
            dv_bed = dv_ice - dv_wall
            wm = 0.5 * (e["w"][:-1] + e["w"][1:])
            dm = 0.5 * ((e["zs"] - e["z"])[:-1] + (e["zs"] - e["z"])[1:])
            dz_seg = dv_bed / (wm * e["dx"])
            dw_seg = dv_wall / ((dm + dz_seg) * e["dx"])     # walls retreat at the depth AFTER bed lowering, so the area change is exactly w dz + D' dw
            dz = np.zeros(e["n"]); dwk = np.zeros(e["n"])
            dz[:-1] += 0.5 * dz_seg; dz[1:] += 0.5 * dz_seg; dz[0] += 0.5 * dz_seg[0]; dz[-1] += 0.5 * dz_seg[-1]
            dwk[:-1] += 0.5 * dw_seg; dwk[1:] += 0.5 * dw_seg; dwk[0] += 0.5 * dw_seg[0]; dwk[-1] += 0.5 * dw_seg[-1]
            a_old = (e["zs"] - e["z"]) * e["w"]
            e["z"] = e["z"] - dz
            e["w"] = e["w"] + dwk
            a_new = (e["zs"] - e["z"]) * e["w"]
            v_geom += np.sum(0.5 * (a_new - a_old)[:-1] * e["dx"] + 0.5 * (a_new - a_old)[1:] * e["dx"])
            v_energy += dv_ice.sum()
        return v_energy, v_geom

    # ----------------------------------------------------------------------------- driver
    def run(self, v0, t0, t_end, dt_geom, inflow, inflow_temp, rtol=1e-9, atol=1e-10, evolve=True, record=None):
        y = self.reset_state(v0, t0)
        self._area0 = [((e["zs"] - e["z"]) * e["w"]).copy() for e in self.edges]      # initial channel section areas (for the geometry-derived ice inventory)
        v_init = np.array(v0, float).copy(); s_init = (np.array(v0) * np.array(t0)).copy()
        ts = [0.0]; hist = [self._snapshot(0.0, y, inflow, inflow_temp)]
        ve_tot = vg_tot = 0.0
        t = 0.0
        while t < t_end - 1e-9:
            t1 = min(t + dt_geom, t_end)
            h_before = y[self.i_H].copy(); w_before = y[self.i_W].copy()
            sol = solve_ivp(lambda tt, yy: self.rhs(tt, yy, inflow, inflow_temp), (t, t1), y, method="LSODA", rtol=rtol, atol=atol)
            if not sol.success:
                raise RuntimeError(f"ODE failure at t={t}: {sol.message}")
            y = sol.y[:, -1].copy()
            if evolve and self.melting:
                ve, vg = self.update_geometry(y[self.i_H] - h_before, y[self.i_W] - w_before)
                ve_tot += ve; vg_tot += vg
                for ie, e in enumerate(self.edges):
                    src = e["a"]
                    if src < self.NL:
                        assert e["z"][0] > self.lakes[src]["z_bottom"] + 1e-6, "sill incised below lake bottom: lake would drain completely (beyond model scope)"
            t = t1
            ts.append(t); hist.append(self._snapshot(t, y, inflow, inflow_temp))
        self.final_state = y
        budgets = self.budgets(y, v_init, s_init)
        vg_t = np.array([h["ice_vol_geometry"] for h in hist]); ve_t = np.array([h["ice_vol_energy"] for h in hist])
        d_kg = RHO_I * (vg_t - ve_t)
        budgets.update(geometry_vs_latent_kg_abs_max=float(np.max(np.abs(d_kg))),
                       geometry_vs_latent_rel_max=float(np.max(np.abs(vg_t - ve_t) / np.maximum(np.abs(ve_t), 1e-30)[0:] ) if ve_t.max() > 0 else 0.0),
                       geometry_vs_latent_rel_of_final=float(np.max(np.abs(vg_t - ve_t)) / max(ve_t[-1], 1e-30)),
                       ice_vol_geometry_series=vg_t.tolist(), ice_vol_energy_series=ve_t.tolist())
        budgets.update(ice_volume_energy=ve_tot, ice_volume_geometry=vg_tot, ice_mass_from_energy_update_kg=RHO_I * ve_tot,
                       ice_mass_from_geometry_kg=RHO_I * vg_tot, geometry_energy_volume_ratio=(vg_tot / ve_tot) if ve_tot else float("nan"))
        return dict(t=np.array(ts), hist=hist, budgets=budgets, y=y)

    # ----------------------------------------------------------------------------- diagnostics
    def _snapshot(self, t, y, inflow, inflow_temp):
        v = np.maximum(y[self.i_V], 0.0)
        stage = self.stages(v)
        temps = np.where(v > 1e-9, y[self.i_S] / np.maximum(v, 1e-9), 0.0)
        q_edges = []
        for e in self.edges:
            st = self.edge_state(e, stage[e["a"]], stage[e["b"]], temps[e["a"]] if e["a"] < self.NL else 0.0, temps[e["b"]] if e["b"] < self.NL else 0.0)
            q_edges.append(st["Q"])
        # ice inventory two ways at every snapshot: from the channel GEOMETRY (section-area change x dx, trapezoid) and from the
        # LATENT HEAT delivered (cumulative heat / (rho_i L_f)); their difference is the geometry-reconstruction error
        vg = 0.0
        for e, a0 in zip(self.edges, self._area0):
            da = (e["zs"] - e["z"]) * e["w"] - a0
            vg += float(np.sum(0.5 * (da[:-1] + da[1:])) * e["dx"])
        ve = float(y[self.i_H].sum() / (RHO_I * LF))
        return dict(t=t, V=v.copy(), stage=stage.copy(), T=temps.copy(), Q=np.array(q_edges), C=y[self.i_C].copy(),
                    Win=float(y[self.i_Win]), Qm=float(y[self.i_Qm]), heat_to_ice=float(y[self.i_H].sum()), Phi=float(y[self.i_Phi]),
                    ice_vol_geometry=vg, ice_vol_energy=ve,
                    sill=np.array([e["z"][0] if e["dir0"] > 0 else e["z"][-1] for e in self._edges_dir()]),
                    z=[e["z"].copy() for e in self.edges], w=[e["w"].copy() for e in self.edges])

    def _edges_dir(self):
        for e in self.edges:
            e["dir0"] = 1
        return self.edges

    def budgets(self, y, v_init, s_init):
        rc = RHO_W * CP_W
        v = y[self.i_V]; s = y[self.i_S]
        water_res = v.sum() + y[self.i_C].sum() - v_init.sum() - y[self.i_Win] - y[self.i_Qm]
        water_scale = max(abs(y[self.i_Win]) + abs(y[self.i_Qm]) + v_init.sum() + v.sum() + y[self.i_C].sum(), 1e-30)
        lhs = rc * s.sum() + y[self.i_E].sum() + y[self.i_H].sum()
        rhs = rc * s_init.sum() + y[self.i_Ein] + y[self.i_Phi]
        energy_res = lhs - rhs
        energy_scale = max(abs(rc * s_init.sum()) + abs(y[self.i_Ein]) + abs(y[self.i_Phi]) + abs(rc * s.sum()) + abs(y[self.i_E].sum()) + abs(y[self.i_H].sum()), 1e-30)
        m_in = RHO_W * y[self.i_Win]; m_out = RHO_W * y[self.i_C].sum(); m_stored_change = RHO_W * (v.sum() - v_init.sum())
        m_melt_water = RHO_W * y[self.i_Qm]; m_melt_energy = y[self.i_H].sum() / LF
        mass_res = m_stored_change + m_out - m_in - m_melt_water
        mass_scale = max(abs(m_in) + abs(m_melt_water) + abs(m_out) + abs(m_stored_change), 1e-30)
        return dict(mass_residual_kg=mass_res, mass_scale_kg=mass_scale, mass_residual_rel=mass_res / mass_scale,
                    mass_in_kg=m_in, mass_out_kg=m_out, mass_stored_change_kg=m_stored_change,
                    ice_mass_from_meltwater_kg=m_melt_water, ice_mass_from_heat_kg=m_melt_energy,
                    water_residual_m3=water_res, water_scale_m3=water_scale, water_residual_rel=water_res / water_scale,
                    energy_residual_J=energy_res, energy_scale_J=energy_scale, energy_residual_rel=energy_res / energy_scale,
                    delivered_m3=y[self.i_C].copy(), input_m3=y[self.i_Win], melt_water_m3=y[self.i_Qm], dissipation_J=y[self.i_Phi],
                    dissipation_channel_J=y[self.i_PhiC], dissipation_exit_J=y[self.i_Phi] - y[self.i_PhiC], sensible_out_J=y[self.i_E].sum(), sensible_in_J=y[self.i_Ein],
                    sensible_storage_change_J=rc * (s.sum() - s_init.sum()),
                    heat_to_ice_J=y[self.i_H].sum(), ice_melted_m3=y[self.i_H].sum() / (RHO_I * LF))
