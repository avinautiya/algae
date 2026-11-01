"""Cost benchmark (M3): short strip runs at coarse (dx=0.05) and medium (dx=0.025) resolution; `python -m ice_cliff_survival.reactivation.benchmark OUT.json`."""
import json, sys, time, resource, subprocess
def one(dx, days, dtm):
    from . import cases
    t = time.time(); r = cases.run("strip", dx=dx, dtm=dtm, days=days, removal_day=days / 2)
    rec = r["records"]; wall = time.time() - t; spin = rec[0]["wall_s"]
    return dict(dx=dx, cells=r["cells"], days=days, dtm_h=dtm / 3600, wall_total_s=wall, wall_spinup_s=spin, wall_per_sim_hour_s=(wall - spin) / (days * 24), hyd_steps=r["hyd_steps"], dt_limit_min_s=min(x["dt_limit"] for x in rec), dt_limit_max_s=max(x["dt_limit"] for x in rec), maxrss_MB=r["maxrss_MB"],
                hyd_time_per_update_s=sum(x["hyd_time"] for x in rec[1:]) / (len(rec) - 1), steps_per_update=r["hyd_steps"] / max(1, len(rec) - 1))
if __name__ == "__main__":
    out = []
    for dx in (0.05, 0.025):
        out.append(one(dx, 3.0, 6 * 3600.0)); print(out[-1], flush=True)
    json.dump(out, open(sys.argv[1], "w"), indent=1)
