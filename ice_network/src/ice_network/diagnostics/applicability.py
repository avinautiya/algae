"""Inexpensive automatic applicability diagnostics for a finished run. They EXPOSE where the stated assumptions are strained; they never
change the solution and there is no universal pass/fail threshold (ratios are reported with their meaning).

Quantities, per snapshot and active edge, using the geometry stored at that snapshot:
  tau_transit   = sum_k dx / v_k                    water transit time along the reach (s)
  tau_wave      = L / sqrt(g h_mean)                gravity-wave crossing time, the time for the channel to adjust hydraulically (s)
  tau_stage     = (H_a - z_0) / |dH_a/dt|           time over which the feeding node's head over the sill changes (finite differences of snapshots)
  r_qs          = tau_wave / tau_stage              quasi-steady hydraulics is a better description the smaller this is; meaning only, no threshold
  channel_over_node = V_channel / V_node            channel water volume relative to the feeding node (neglected storage)
  dV_channel_over_throughput                        change of channel volume between snapshots relative to the volume that flowed through
  melt_over_flow = Q_melt / Q                       melt-generated discharge relative to through-flow (constant-Q-along-edge assumption)
  froude        = v / sqrt(g h) at stations; n_choked = stations set to critical depth by the control-point rule (supercritical flow and
                  hydraulic jumps are NOT represented)
  Re, Pr, D_h/L, and the correlation assumptions. A Reynolds number inside the nominal range of a pipe correlation does not establish
  validity for a melting open channel (free surface, ice roughness/scallops, evolving cross-section are not in the correlation).
"""
import numpy as np

from ..constants import G, NU_W, PR_W

CLOSURE_ASSUMPTIONS = (
    "Dittus-Boelter Nu = 0.023 Re^0.8 Pr^0.3 is a smooth-pipe, fully developed, turbulent, constant-property correlation (nominal Re about "
    "1e4 to 1.2e5, 0.7 < Pr < 160). It is applied here with D_h = 4A/P to an open rectangular channel with a free surface and a melting, "
    "possibly rough, evolving wall. Entrance effects, free-surface effects and ice roughness are not represented; its use is an assumption, "
    "not a validated closure.")


def diagnose(net, result):
    hist = result["hist"]; t = np.array(result["t"])
    rows = []
    for ie, e in enumerate(net.edges):
        for it, snap in enumerate(hist):
            ed = dict(e); ed["z"] = snap["z"][ie]; ed["w"] = snap["w"][ie]
            a, b = e["a"], e["b"]
            ta = snap["T"][a] if a < net.NL else 0.0; tb = snap["T"][b] if b < net.NL else 0.0
            st = net.edge_state(ed, snap["stage"][a], snap["stage"][b], ta, tb)
            if st["Q"] == 0.0:
                continue
            q = st["Qabs"]; d = st["dir"]
            wf = ed["w"] if d > 0 else ed["w"][::-1]
            h = st["h"]; dx = e["dx"]
            v = q / (wf * h)
            src = a if d > 0 else b
            fr = v / np.sqrt(G * h)
            L = dx * (e["n"] - 1)
            tau_transit = float(np.sum(dx / (0.5 * (v[:-1] + v[1:]))))
            tau_wave = float(L / np.sqrt(G * h.mean()))
            vol_ch = float(np.sum(0.5 * (wf[:-1] * h[:-1] + wf[1:] * h[1:]) * dx))
            head = snap["stage"][src] - (ed["z"][0] if d > 0 else ed["z"][-1])
            if it > 0:
                dt = t[it] - t[it - 1]
                dH = abs(snap["stage"][src] - hist[it - 1]["stage"][src]) / dt
                tau_stage = head / dH if dH > 0 else np.inf
            else:
                tau_stage = np.nan
            vnode = snap["V"][src] if src < net.NL else np.nan
            D_h = 4 * (wf * h) / (wf + 2 * h)
            re = v * D_h / NU_W
            rows.append(dict(edge=ie, t=float(t[it]), Q=float(q), tau_transit_s=tau_transit, tau_wave_s=tau_wave, tau_stage_s=float(tau_stage),
                             r_qs=float(tau_wave / tau_stage) if np.isfinite(tau_stage) and tau_stage > 0 else float("nan"),
                             channel_volume_m3=vol_ch, channel_over_node=float(vol_ch / vnode) if vnode > 0 else float("nan"),
                             melt_over_flow=float(st["q_melt"] / q), froude_max=float(fr.max()), n_stations_critical_control=int(st["choked"].sum()),
                             n_stations_froude_ge_1=int((fr >= 1.0).sum()), Re_min=float(re.min()), Re_max=float(re.max()), Pr=PR_W,
                             Dh_over_L_max=float(D_h.max() / L)))
    out = dict(rows=rows, closure_assumptions=CLOSURE_ASSUMPTIONS)
    if rows:
        def agg(key, f):
            vals = [r[key] for r in rows if np.isfinite(r[key])]
            return float(f(vals)) if vals else float("nan")
        out["summary"] = dict(
            n_active_edge_snapshots=len(rows),
            r_qs_max=agg("r_qs", max), tau_wave_s_max=agg("tau_wave_s", max), tau_transit_s_max=agg("tau_transit_s", max),
            channel_over_node_max=agg("channel_over_node", max), melt_over_flow_max=agg("melt_over_flow", max),
            froude_max=agg("froude_max", max), snapshots_with_critical_control=int(sum(r["n_stations_critical_control"] > 0 for r in rows)),
            snapshots_with_froude_ge_1=int(sum(r["n_stations_froude_ge_1"] > 0 for r in rows)), Re_min=agg("Re_min", min), Re_max=agg("Re_max", max),
            Re_outside_nominal_1e4_1p2e5=bool(agg("Re_min", min) < 1e4 or agg("Re_max", max) > 1.2e5))
    return out
