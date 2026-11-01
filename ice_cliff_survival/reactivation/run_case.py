"""python -m ice_cliff_survival.reactivation.run_case OUTDIR NAME layout dx dtm_h days [removal_day|none]  -- one run, atomic JSON result, skips if the result exists (restart-safe at run granularity)."""
import json, os, sys, time
from . import cases
def main(a):
    out, name, layout, dx, dtm_h, days, rem = a[0], a[1], a[2], float(a[3]), float(a[4]), float(a[5]), (None if a[6] == "none" else float(a[6]))
    path = os.path.join(out, name + ".json")
    if os.path.exists(path): print("exists, skipped", path); return 0
    os.makedirs(out, exist_ok=True)
    try:
        r = cases.run(layout, dx=dx, dtm=dtm_h * 3600.0, days=days, removal_day=rem); r["status"] = "complete"; r["name"] = name
    except Exception as e:
        import traceback; r = dict(status="FAILED", name=name, error=repr(e), trace=traceback.format_exc(), layout=layout, dx=dx)
    json.dump(r, open(path + ".tmp", "w"), default=float); os.replace(path + ".tmp", path); print(name, r["status"]); return 0
if __name__ == "__main__": sys.exit(main(sys.argv[1:]))
