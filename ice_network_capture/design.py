"""Development and test configurations for the two-outlet capture study and their conversion to campaign cases (ice_network_eval/campaign.py runner: per-run dirs, progress, limits).
Configurations are dicts over rtm.DEFAULT (12 h horizon). Treatments: thermal (production network), dissipation (T2: all channel dissipation melts, SaDS-type surrogate), fixed (no melting)."""
import json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "ice_network_eval"))
import rtm
H12 = dict(duration_h=12.0)
SYM = dict(w2=0.4, L2=100.0, S2=0.005)       # route 2 identical to route 1
DEV = {   # id: (description, overrides)
    "D01": ("symmetric routes, sill offset 0.02 m, 1 C, steady", dict(SYM, offset=0.02)),
    "D02": ("symmetric, offset 0.05 m, 2 C, steady", dict(SYM, offset=0.05, T=2.0)),
    "D03": ("symmetric, offset 0.05 m, 1 C, pulsed 50 %, 3 h", dict(SYM, offset=0.05, amp=0.5, period_h=3.0)),
    "D04": ("S1-like: steep short route 2, offset 0.05, 1 C", dict(offset=0.05)),
    "D05": ("S2-like: offset 0.02, 1 C", dict(offset=0.02)),
    "D06": ("S3-like: offset 0.10 m, inflow 0.08", dict(offset=0.10, I=0.08)),
    "D07": ("D04 at 2 C", dict(offset=0.05, T=2.0)),
    "D08": ("D04 at 0.5 C", dict(offset=0.05, T=0.5)),
    "D09": ("route 2 less steep (0.01), offset 0.02", dict(offset=0.02, S2=0.01)),
    "D10": ("route 2 long (120 m), offset 0.02", dict(offset=0.02, L2=120.0)),
    "D11": ("route 2 wide (0.8 m), offset 0.05", dict(offset=0.05, w2=0.8)),
    "D12": ("route 1 steep-short, route 2 flat-long, offset 0.05", dict(offset=0.05, S1=0.02, L1=60.0, S2=0.005, L2=100.0)),
    "D13": ("pulsed 60 %, 2 h, offset 0.02", dict(offset=0.02, amp=0.6, period_h=2.0)),
    "D14": ("pulsed 40 %, 6 h, offset 0.05", dict(offset=0.05, amp=0.4, period_h=6.0)),
    "D15": ("low inflow 0.02, 2 C, offset 0.02", dict(offset=0.02, I=0.02, T=2.0)),
    "D16": ("route 2 the lower sill (-0.03), steep route 2", dict(offset=-0.03)),
}
TEST = {   # frozen BEFORE any test run (see TEST_FREEZE.md): both sides of the predicted boundaries
    "T01": ("geom A, small lake (300 m2), 3 C, gap 0.15 m (below QS threshold 0.166)", dict(area=300.0, T=3.0, offset=0.15)),
    "T02": ("geom A, small lake, 3 C, gap 0.18 m (between QS threshold and horizon boundary)", dict(area=300.0, T=3.0, offset=0.18)),
    "T03": ("geom A, small lake, 3 C, gap 0.21 m (above both)", dict(area=300.0, T=3.0, offset=0.21)),
    "T04": ("geom B (route-2 slope 0.01), small lake, 3 C, gap 0.085 (below QS threshold 0.0965)", dict(S2=0.01, area=300.0, T=3.0, offset=0.085)),
    "T05": ("geom B, small lake, 3 C, gap 0.115 (between)", dict(S2=0.01, area=300.0, T=3.0, offset=0.115)),
    "T06": ("geom B, small lake, 3 C, gap 0.15 (above both)", dict(S2=0.01, area=300.0, T=3.0, offset=0.15)),
    "T07": ("geom A, 2 C, gap 0.05 (switching expected)", dict(T=2.0, offset=0.05)),
    "T08": ("geom A, 2 C, gap 0.065", dict(T=2.0, offset=0.065)),
    "T09": ("geom A, 2 C, gap 0.085", dict(T=2.0, offset=0.085)),
    "T10": ("geom A, 2 C, gap 0.10", dict(T=2.0, offset=0.10)),
    "T11": ("geom A, 2 C, gap 0.05, pulsed 50 % / 3 h", dict(T=2.0, offset=0.05, amp=0.5, period_h=3.0)),
    "T12": ("geom A, 2 C, gap 0.085, pulsed 50 % / 3 h", dict(T=2.0, offset=0.085, amp=0.5, period_h=3.0)),
    "T13": ("symmetric routes, 3 C, gap 0.02", dict(SYM, T=3.0, offset=0.02)),
    "T14": ("symmetric steep-short routes (S 0.05, L 40 m), 3 C, gap 0.02", dict(w1=0.4, L1=40.0, S1=0.05, w2=0.4, L2=40.0, S2=0.05, T=3.0, offset=0.02)),
    "T15": ("route 1 steep-short, route 2 flat-long, 2 C, gap 0.05 (no overtaking expected)", dict(S1=0.02, L1=60.0, S2=0.005, L2=100.0, T=2.0, offset=0.05)),
    "T16": ("geom A, 0.3 C, gap 0.05 (cold: incision too slow to overtake?)", dict(T=0.3, offset=0.05)),
}
TEST_REF = ("T01", "T02", "T04", "T05", "T07", "T13")        # informative cases for the independent resolved-model comparison
TREATMENTS = {"thermal": {}, "dissipation": dict(melt_mode="dissipation_all"), "fixed": dict(melting=False)}


def spec(cid, p, treat, kind, stations=17, ncell=200):
    q = rtm.full(dict(H12, **p)); sid = cid
    patch = dict(lake=dict(area_m2=q["area"], initial_stage_below_sill1_m=q["deficit"]), duration_h=q["duration_h"],
                 channels=[dict(width_m=q["w1"], length_m=q["L1"], slope=q["S1"], crest_length_m=q["crest"]), dict(width_m=q["w2"], length_m=q["L2"], slope=q["S2"], crest_length_m=q["crest"])])
    scn = dict(id=sid, sill_offset_m=q["offset"], inflow_m3_s=q["I"], inflow_T_C=q["T"], pulse_amp=q["amp"], pulse_period_s=q["period_h"] * 3600.0)
    ov = dict(n_man=q["n_man"], nu_factor=q["nu_factor"])
    t = TREATMENTS[treat]
    if kind == "reduced": ov.update(t)
    else:
        if treat == "dissipation": return None          # the resolved reference has no dissipation-only closure
        ov.update({k: v for k, v in t.items() if k == "melting"})
    base = dict(kind=kind, scn=sid, scn_override=scn, cfg_patch=patch, overrides=ov, params=q, version="capture_v1")
    if kind == "reduced": return dict(base, stations=stations, dt_geom=450 if stations == 17 else 225, label=f"cap_{cid}_{treat}_r{stations}")
    return dict(base, ncell=ncell, label=f"cap_{cid}_{treat}_ref{ncell}")


def cases(table, treatments=("thermal", "dissipation", "fixed"), ref_ids=(), ref_cells=(200,), r33_ids=()):
    out = []
    for cid, (_, p) in table.items():
        for tr in treatments: out.append(spec(cid, p, tr, "reduced"))
        if cid in r33_ids: out.append(spec(cid, p, "thermal", "reduced", stations=33))
        if cid in ref_ids:
            for n in ref_cells:
                for tr in ("thermal", "fixed"): out.append(spec(cid, p, tr, "reference", ncell=n))
    return [c for c in out if c]


if __name__ == "__main__":
    which = sys.argv[1]
    if which == "dev":
        c = cases(DEV); json.dump(c, open(os.path.join(HERE, "cases_dev.json"), "w"), indent=1); print(len(c), "dev cases")
    if which == "test":
        c = cases(TEST, ref_ids=TEST_REF, ref_cells=(200, 400), r33_ids=TEST_REF); json.dump(c, open(os.path.join(HERE, "cases_test.json"), "w"), indent=1); print(len(c), "test cases")
