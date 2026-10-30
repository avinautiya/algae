"""Reduced network versus independent resolved reference on the frozen competing-outlets scenarios (competing_outlets_v1.json)."""
import json, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, "..")
sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "ice_network", "src"))
from ice_network.evolution.network import Network
from ice_network_reference.solver import RefNetwork

CFG = json.load(open(os.path.join(HERE, "competing_outlets_v1.json")))
T_END = CFG["duration_h"] * 3600.0
BANK = CFG["bank_height_m"]
G_COMMON = dict(n_man=CFG["physics_shared_by_both_models"]["manning_n"], ke=CFG["physics_shared_by_both_models"]["entrance_loss_ke"])


def bed(ch, x, sill):
    c = ch["crest_length_m"]
    return sill - ch["slope"] * np.maximum(np.asarray(x) - c, 0.0)


def sills(scn):
    s1 = CFG["channels"][0]["sill_m"]
    return [s1, s1 + scn["sill_offset_m"]]


def run_reduced(scn, n_st, dt_geom, msub=16, rtol=1e-7, nu_factor=1.0, melting=True, heat_dissipation=True, sample_dt=60.0):
    lk = CFG["lake"]; sl = sills(scn)
    edges = []
    for k, ch in enumerate(CFG["channels"]):
        x = np.linspace(0.0, ch["length_m"], n_st); z = bed(ch, x, sl[k])
        edges.append(dict(a=0, b=1 + k, z=z, w=np.full(n_st, ch["width_m"]), zs=z + BANK, length=ch["length_m"], msub=msub, **G_COMMON))
    outlets = [dict(h_fixed=float(e["z"][-1]) - 5.0) for e in edges]
    net = Network([dict(z_bottom=lk["bottom_m"], a_ref=lk["area_m2"], h_ref=1.0, q=0.0)], outlets, edges, melting=melting, nu_factor=nu_factor, heat_dissipation=heat_dissipation)
    v0 = lk["area_m2"] * (sl[0] - lk["initial_stage_below_sill1_m"] - lk["bottom_m"])
    I, T = scn["inflow_m3_s"], scn["inflow_T_C"]
    # COST run: no intermediate sampling (sampling is a diagnostic overhead, not part of the method); count discharge-law evaluations
    import ice_network.evolution.network as NW
    orig = NW.edge_flow; calls = [0]
    def counted(*a):
        calls[0] += 1; return orig(*a)
    NW.edge_flow = counted
    t0 = time.perf_counter()
    net.run([v0], [T], T_END, dt_geom, lambda t: np.array([I]), lambda t: np.array([T]), rtol=rtol, atol=1e-9)
    wall = time.perf_counter() - t0; nfev = calls[0]; NW.edge_flow = orig
    # OUTPUT run: identical settings with dense output samples every sample_dt (does not change the solution)
    net2 = Network([dict(z_bottom=lk["bottom_m"], a_ref=lk["area_m2"], h_ref=1.0, q=0.0)], outlets, [dict(e) for e in [dict(a=0, b=1 + k, z=bed(ch, np.linspace(0.0, ch["length_m"], n_st), sl[k]), w=np.full(n_st, ch["width_m"]), zs=bed(ch, np.linspace(0.0, ch["length_m"], n_st), sl[k]) + BANK, length=ch["length_m"], msub=msub, **G_COMMON) for k, ch in enumerate(CFG["channels"])]], melting=melting, nu_factor=nu_factor, heat_dissipation=heat_dissipation)
    r = net2.run([v0], [T], T_END, dt_geom, lambda t: np.array([I]), lambda t: np.array([T]), rtol=rtol, atol=1e-9, sample_dt=sample_dt)
    net = net2
    t = np.array(r["t"]); h = r["hist"]
    return dict(t=t, Q1=np.array([abs(s["Q"][0]) for s in h]), Q2=np.array([abs(s["Q"][1]) for s in h]), stage=np.array([s["stage"][0] for s in h]),
                sill1=np.array([e0 for e0 in [s["z"][0][0] for s in h]]), sill2=np.array([s["z"][1][0] for s in h]),
                delivered=r["budgets"]["delivered_m3"], wall=wall, budgets={k: (float(v) if np.ndim(v) == 0 else None) for k, v in r["budgets"].items()},
                heat=float(r["budgets"]["heat_to_ice_J"]), net=net, result=r, nfev=nfev)


def run_reference(scn, ncell, nu_factor=1.0, melting=True, heat_dissipation=True):
    lk = CFG["lake"]; sl = sills(scn)
    chans = []
    for k, ch in enumerate(CFG["channels"]):
        dx = ch["length_m"] / ncell; x = (np.arange(ncell) + 0.5) * dx; z = bed(ch, x, sl[k])
        chans.append(dict(lake=0, z=z, w=np.full(ncell, ch["width_m"]), zs=z + BANK, length=ch["length_m"], pool=float(bed(ch, ch["length_m"], sl[k])) - 5.0, **G_COMMON))
    net = RefNetwork([dict(z_bottom=lk["bottom_m"], area=lk["area_m2"])], chans, melting=melting, heat_diss=heat_dissipation, nu_factor=nu_factor)
    v0 = lk["area_m2"] * (sl[0] - lk["initial_stage_below_sill1_m"] - lk["bottom_m"])
    I, T = scn["inflow_m3_s"], scn["inflow_T_C"]
    t0 = time.perf_counter()
    r = net.run([v0], [T], T_END, lambda t: np.array([I]), lambda t: np.array([T]), record_every=60.0)
    wall = time.perf_counter() - t0
    rec = r["rec"]
    return dict(t=np.array([x["t"] for x in rec]), Q1=np.array([x["Q_in"][0] for x in rec]), Q2=np.array([x["Q_in"][1] for x in rec]),
                stage=np.array([x["stage"][0] for x in rec]), sill1=np.array([x["sill"][0] for x in rec]), sill2=np.array([x["sill"][1] for x in rec]),
                delivered=rec[-1]["delivered"], wall=wall, nsteps=r["nsteps"], heat=rec[-1]["heat"], diss=rec[-1]["diss"])


def crossing(t, y, thr, rising=True):
    """first time y crosses thr upward (linear interpolation); None if never."""
    idx = np.flatnonzero((y[:-1] <= thr) & (y[1:] > thr))
    if idx.size == 0:
        return None
    i = idx[0]; return float(t[i] + (thr - y[i]) / (y[i + 1] - y[i]) * (t[i + 1] - t[i]))


def last_above(t, y, thr):
    if y[-1] > thr: return None
    idx = np.flatnonzero(y > thr)
    if idx.size == 0: return None
    i = idx[-1]; return float(t[i] + (y[i] - thr) / (y[i] - y[i + 1]) * (t[i + 1] - t[i]))


def metrics(m, scn):
    I = scn["inflow_m3_s"]; t = m["t"]; thr = CFG["metrics"]["activity_threshold_fraction_of_inflow"] * I
    return dict(t_act1_s=crossing(t, m["Q1"], thr), t_act2_s=crossing(t, m["Q2"], thr), t_end_act2_s=last_above(t, m["Q2"], thr),
                t_dom2_s=crossing(t, m["Q2"] - m["Q1"], 0.0),
                active_at_end=[bool(m["Q1"][-1] > thr), bool(m["Q2"][-1] > thr)], peak_Q2_over_I=float(m["Q2"].max() / I),
                sill_lowering1_m=float(m["sill1"][0] - m["sill1"][-1]), sill_lowering2_m=float(m["sill2"][0] - m["sill2"][-1]),
                delivered1_m3=float(m["delivered"][0]), delivered2_m3=float(m["delivered"][1]), heat_to_ice_J=float(m["heat"]),
                stage_end_m=float(m["stage"][-1]), stage_max_m=float(m["stage"].max()))


def series_error(red, ref, I):
    """errors of the reduced discharge series against the reference on the reference's 60 s grid, normalised by the inflow."""
    t = ref["t"]
    out = {}
    for key in ("Q1", "Q2"):
        q = np.interp(t, red["t"], red[key]); d = q - ref[key]
        out[key] = dict(rmse_over_I=float(np.sqrt(np.mean(d ** 2)) / I), max_abs_over_I=float(np.abs(d).max() / I), volume_error_m3=float(np.trapezoid(d, t)))
    out["stage_rmse_m"] = float(np.sqrt(np.mean((np.interp(t, red["t"], red["stage"]) - ref["stage"]) ** 2)))
    return out


def main(out_path, scenario_ids=None, ref_cells=(100, 200, 400)):
    results = {}
    for scn in CFG["scenarios"]:
        if scenario_ids and scn["id"] not in scenario_ids: continue
        print("scenario", scn["id"], flush=True)
        refs = {n: run_reference(scn, n) for n in ref_cells}
        print("  reference done", {n: round(r["wall"], 1) for n, r in refs.items()}, flush=True)
        ref_m = {n: metrics(r, scn) for n, r in refs.items()}
        reds = {}
        for rc in CFG["reduced_model_configurations"]:
            m = run_reduced(scn, rc["stations"], rc["dt_geom_s"])
            reds[rc["id"]] = dict(metrics=metrics(m, scn), discharge_law_evaluations=m["nfev"], series_error_vs_finest_reference=series_error(m, refs[ref_cells[-1]], scn["inflow_m3_s"]), wall_s=m["wall"],
                                  conservation={k: m["budgets"][k] for k in ("water_residual_rel", "mass_residual_rel", "energy_residual_rel", "geometry_vs_latent_rel_max")})
            print("  reduced", rc["id"], round(m["wall"], 2), "s", flush=True)
        results[scn["id"]] = dict(scenario=scn, reference_metrics=ref_m, reference_wall_s={n: r["wall"] for n, r in refs.items()}, reference_steps={n: r["nsteps"] for n, r in refs.items()},
                                  reduced=reds)
        json.dump(results, open(out_path, "w"), indent=1, default=float)
    return results


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "records", "eval_results.json")
    ids = sys.argv[2].split(",") if len(sys.argv) > 2 else None
    cells = tuple(int(x) for x in sys.argv[3].split(",")) if len(sys.argv) > 3 else (100, 200, 400)
    main(out, ids, cells)
