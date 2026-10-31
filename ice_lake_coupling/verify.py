"""Small timestep-refinement campaign, not a scientific performance campaign."""
import argparse
import copy
import hashlib
import json
import math
import platform
import time
from pathlib import Path

from .simulation import case_digest, run_case


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True, help="New JSON file")
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Output exists; choose a new filename")
    root = Path(__file__).resolve().parent
    case = json.loads((root/"cases/synthetic_lake.json").read_text())
    rows = []
    for dt in (60.0, 30.0, 15.0, 7.5):
        run = copy.deepcopy(case)
        run["max_dt_s"] = dt
        start = time.perf_counter()
        _, summary = run_case(run, True)
        rows.append(dict(dt_s=dt, wall_s=time.perf_counter()-start,
            final_stage_m=summary["final_stage_m"], sill_lowering_m=summary["sill_lowering_m"],
            delivered_m3=summary["totals"]["delivered_m3"], exported_J=summary["totals"]["exported_J"],
            max_water_residual_relative=summary["max_water_residual_relative"],
            max_energy_residual_relative=summary["max_energy_residual_relative"]))
    orders = {}
    for key in ("final_stage_m", "sill_lowering_m", "delivered_m3", "exported_J"):
        diffs = [abs(a[key]-b[key]) for a,b in zip(rows, rows[1:])]
        orders[key] = dict(successive_absolute_differences=diffs,
            observed_orders=[math.log(a/b,2) if a>0 and b>0 else None for a,b in zip(diffs,diffs[1:])])
    report = dict(scope="Numerical verification on one synthetic case, not environmental validation",
        python=platform.python_version(), case_sha256=case_digest(case), levels=rows, convergence=orders,
        source_sha256={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest()
                       for p in sorted(root.rglob("*.py"))})
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open("x") as stream:
        json.dump(report,stream,indent=2,allow_nan=False); stream.write("\n")
    print(json.dumps({"output":str(args.output),"convergence":orders},indent=2))


if __name__ == "__main__":
    main()
