"""Development and test configurations of the forcing-sequence study and their conversion to campaign cases (ice_network_eval/campaign.py runner).
Every configuration has episodes A and B of equal duration tau (default 6 h), a common 6 h continuation at the discharge-weighted mean forcing, the same initial state (lake stage 0.10 m above sill 1, lake temperature = discharge-weighted mean temperature),
three forcings (AB, BA, AVG) and the treatments thermal (production network), fixed geometry, dissipation-only incision."""
import json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "ice_network_capture")); sys.path.insert(0, os.path.join(HERE, "..", "ice_network_eval"))
import forcing, rtm
P1 = ((0.03, 2.0), (0.09, 0.2)); P2 = ((0.05, 1.5), (0.05, 0.1)); P3 = ((0.02, 1.0), (0.10, 1.0))
SYM = dict(w2=0.4, L2=100.0, S2=0.005)
def cfg(geom, pair, tau_h=6.0, tail_h=6.0, **kw): return dict(geom=geom, pair=pair, tau_h=tau_h, tail_h=tail_h, **kw)
DEV = {
    "S01": ("geom A, gap 0.05, lake 2000 m2, warm-low vs cold-high, 6 h episodes", cfg(dict(offset=0.05), P1)),
    "S02": ("S01, small lake 300 m2", cfg(dict(offset=0.05, area=300.0), P1)),
    "S03": ("S01, large lake 6000 m2", cfg(dict(offset=0.05, area=6000.0), P1)),
    "S04": ("geom A, gap 0.18 (near capture boundary), lake 300", cfg(dict(offset=0.18, area=300.0), P1)),
    "S05": ("geom B (S2 0.01), gap 0.115 (near boundary), lake 300", cfg(dict(offset=0.115, S2=0.01, area=300.0), P1)),
    "S06": ("symmetric routes, gap 0.02", cfg(dict(SYM, offset=0.02), P1)),
    "S07": ("S01 with temperature contrast only (Q equal)", cfg(dict(offset=0.05), P2)),
    "S08": ("S01 with discharge contrast only (T equal)", cfg(dict(offset=0.05), P3)),
    "S09": ("S01 with 2 h episodes", cfg(dict(offset=0.05), P1, tau_h=2.0)),
    "S10": ("geom A gap 0.10, 12 h episodes", cfg(dict(offset=0.10), P1, tau_h=12.0)),
}
TEST = {   # frozen before any test run; configurations chosen from the RTM scan (scan_rtm.py) by predicted-effect strata; no simulator was run for the selection
    "X01": ("scan config 2: predicted share effect -0.0019", cfg({"w2": 0.8, "offset": 0.12, "area": 1859.4476040701634}, ((0.03085806418120406, 2.7593023466698443), (0.03085806418120406, 0.4395907665106825)), tau_h=3.0)),
    "X02": ("scan config 13: predicted share effect +0.0000", cfg({"w2": 0.8, "offset": 0.12, "area": 2109.713666937508}, ((0.10144654627240196, 2.0027144088342523), (0.10144654627240196, 2.0027144088342523)), tau_h=6.0)),
    "X03": ("scan config 23: predicted share effect -0.0001", cfg({"w2": 0.8, "offset": 0.05, "area": 2340.8978701912793}, ((0.04028063317107432, 1.208842411456948), (0.10558812573934266, 0.5114721226072387)), tau_h=3.0)),
    "X04": ("scan config 42: predicted share effect +0.0240", cfg({"offset": 0.05, "area": 2536.242032627394}, ((0.043294682635434054, 2.891316597848004), (0.07138528046417403, 0.16238176086501102)), tau_h=6.0)),
    "X05": ("scan config 77: predicted share effect -0.0223", cfg({"w2": 0.8, "offset": 0.16, "area": 1298.3559331781528}, ((0.03981758258419583, 2.7279354470287314), (0.03981758258419583, 0.15085519825946422)), tau_h=6.0)),
    "X06": ("scan config 74: predicted share effect -0.0156", cfg({"w2": 0.8, "offset": 0.16, "area": 1636.5062661692286}, ((0.04283682516212432, 1.2209866115801564), (0.06654471417996842, 0.589436676331678)), tau_h=6.0)),
    "X07": ("scan config 51: predicted share effect -0.0180", cfg({"w2": 0.8, "offset": 0.12, "area": 3683.377203914502}, ((0.02470681570213438, 1.9547128454212397), (0.11199699481579614, 0.27748237566592293)), tau_h=3.0)),
    "X08": ("scan config 43: predicted share effect +0.0093", cfg({"S2": 0.01, "offset": 0.05, "area": 4468.848189901482}, ((0.020251525040711153, 2.820531790048876), (0.06543425447272558, 0.2531366281238666)), tau_h=12.0)),
    "X09": ("scan config 64: predicted share effect -0.0471", cfg({"L2": 120.0, "offset": 0.2, "area": 1125.4524383456217}, ((0.020634312269279746, 2.123019044796103), (0.11664444468420922, 0.11895698634794147)), tau_h=12.0)),
    "X10": ("scan config 96: predicted share effect +0.0451", cfg({"offset": 0.05, "area": 6042.329573281198}, ((0.0658895032695474, 1.944750319637875), (0.0658895032695474, 0.5872767891467314)), tau_h=12.0)),
}
TEST_REF = ('X02', 'X04', 'X07', 'X08', 'X09', 'X10')   # informative cases for the independent resolved-model comparison: X02 (scan 13, low), X07 (scan 51, near tie), X08 (scan 43, near tie, large heat effect), X09 (scan 64, 12 h), X10 (scan 96, 12 h); X04 medium is added for ref400
TEST_REF400 = ('X07', 'X04', 'X10')
TEST_CLOSURE = ('X04', 'X08', 'X09')
TEST_TOL = ('X04', 'X09')
TREATMENTS = {"thermal": {}, "fixed": dict(melting=False), "dissipation": dict(melt_mode="dissipation_all")}


def physical(c, ramp_s=900.0):
    A, B = c["pair"]; seq, mean = forcing.sequences(A, B, c["tau_h"] * 3600, c["tau_h"] * 3600, c["tail_h"] * 3600, ramp_s)
    return seq, mean


def spec(cid, c, forcing_name, treat, kind, stations=17, ncell=200, overrides_extra=None):
    seq, mean = physical(c); q = rtm.full(dict(c["geom"], deficit=-0.10, duration_h=mean["duration"] / 3600.0))
    patch = dict(lake=dict(area_m2=q["area"], initial_stage_below_sill1_m=q["deficit"]), duration_h=q["duration_h"],
                 channels=[dict(width_m=q["w1"], length_m=q["L1"], slope=q["S1"], crest_length_m=q["crest"]), dict(width_m=q["w2"], length_m=q["L2"], slope=q["S2"], crest_length_m=q["crest"])])
    scn = dict(id=f"{cid}{forcing_name}", sill_offset_m=q["offset"], inflow_m3_s=mean["Q_bar"], inflow_T_C=mean["T_bar"], lake_T0=mean["T_bar"], segments=seq[forcing_name], ramp_s=900.0)
    ov = dict(n_man=q["n_man"], nu_factor=q["nu_factor"]); t = TREATMENTS[treat]
    if kind == "reduced": ov.update(t)
    else:
        if treat == "dissipation": return None
        ov.update({k: v for k, v in t.items() if k == "melting"})
    if overrides_extra: ov.update(overrides_extra)
    base = dict(kind=kind, scn=scn["id"], scn_override=scn, cfg_patch=patch, overrides=ov, params=dict(q, **{"cid": cid, "forcing": forcing_name, "mean": mean}), version="sequence_v1")
    if kind == "reduced": return dict(base, stations=stations, dt_geom=450 if stations == 17 else 225, label=f"seq_{cid}_{forcing_name}_{treat}_r{stations}")
    return dict(base, ncell=ncell, label=f"seq_{cid}_{forcing_name}_{treat}_ref{ncell}")


def cases(table, treatments=("thermal", "fixed", "dissipation"), ref_ids=(), ref_cells=(200,), r33_ids=(), ref400_ids=()):
    out = []
    for cid, (_, c) in table.items():
        for fn in ("AB", "BA", "AVG"):
            for tr in treatments: out.append(spec(cid, c, fn, tr, "reduced"))
            if cid in r33_ids: out.append(spec(cid, c, fn, "thermal", "reduced", stations=33))
            if cid in ref_ids:
                for n in ref_cells: out.append(spec(cid, c, fn, "thermal", "reference", ncell=n))
            if cid in ref400_ids and fn != "AVG": out.append(spec(cid, c, fn, "thermal", "reference", ncell=400))
    return [c for c in out if c]


if __name__ == "__main__":
    if sys.argv[1] == "test":
        out = cases(TEST, treatments=("thermal", "fixed"), ref_ids=TEST_REF, ref400_ids=TEST_REF400, r33_ids=())
        for cid in ("X07", "X04", "X10"):    # r33 for the effect's discretisation uncertainty (AB and BA only)
            for fn in ("AB", "BA"): out.append(spec(cid, TEST[cid][1], fn, "thermal", "reduced", stations=33))
        for cid in TEST_CLOSURE:             # closure multipliers on the thermal network (matched in both orders and the control)
            for tag, ov in (("nu0.5", dict(nu_factor=0.5)), ("nu2", dict(nu_factor=2.0))):
                for fn in ("AB", "BA", "AVG"):
                    o = spec(cid, TEST[cid][1], fn, "thermal", "reduced", overrides_extra=ov); o["label"] = o["label"].replace("_thermal_", f"_thermal{tag}_"); out.append(o)
        for cid in TEST_TOL:                 # corrected energy tolerance (atol_energy 1e-3 J) against the legacy setting
            for fn in ("AB", "BA"):
                o = spec(cid, TEST[cid][1], fn, "thermal", "reduced", overrides_extra=dict(atol_energy=1e-3)); o["label"] = o["label"].replace("_thermal_", "_thermalTOL_"); out.append(o)
        json.dump(out, open(os.path.join(HERE, "cases_test.json"), "w"), indent=1); print(len(out), "test cases")
    if sys.argv[1] == "dev":
        c = cases(DEV); json.dump(c, open(os.path.join(HERE, "cases_dev.json"), "w"), indent=1); print(len(c), "dev cases")
