"""Run a versioned synthetic comparison without optional dependencies."""
import argparse
import copy
import json
import platform
import sys
from pathlib import Path

from .kernel import DomainError
from .simulation import case_digest, run_case, validate_case, write_run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True, help="New directory; existing results are never overwritten")
    parser.add_argument("--dt", type=float, help="Override maximum timestep in seconds, recorded in manifest")
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Output path exists; choose a new directory to preserve results")
    try:
        case = json.loads(args.case.read_text())
        original_hash = case_digest(case)
        if args.dt is not None:
            case = copy.deepcopy(case)
            case["max_dt_s"] = args.dt
        validate_case(case)
        args.output.mkdir(parents=True)
        manifest = dict(original_case_sha256=original_hash, executed_case_sha256=case_digest(case),
                        python=platform.python_version(), argv=sys.argv, executed_case=case)
        (args.output/"manifest.json").write_text(json.dumps(manifest, indent=2, allow_nan=False)+"\n")
        summaries = {}
        for name, melting in (("fixed_control", False), ("thermal_control", True)):
            rows, summary = run_case(case, melting)
            write_run(args.output/name, rows, summary)
            summaries[name] = summary
        fixed, thermal = summaries["fixed_control"], summaries["thermal_control"]
        comparison = dict(scope="Synthetic illustration only; no field skill, novelty or hazard validation",
            thermal_minus_fixed_final_stage_m=thermal["final_stage_m"]-fixed["final_stage_m"],
            thermal_minus_fixed_delivered_m3=thermal["totals"]["delivered_m3"]-fixed["totals"]["delivered_m3"],
            thermal_minus_fixed_heat_export_J=thermal["totals"]["exported_J"]-fixed["totals"]["exported_J"],
            thermal_sill_lowering_m=thermal["sill_lowering_m"],
            max_relative_water_residual=max(s["max_water_residual_relative"] for s in summaries.values()),
            max_relative_energy_residual=max(s["max_energy_residual_relative"] for s in summaries.values()))
        (args.output/"comparison.json").write_text(json.dumps(comparison, indent=2, allow_nan=False)+"\n")
        print(json.dumps(comparison, indent=2))
    except (DomainError, ValueError, KeyError, TypeError, OSError) as exc:
        if args.output.is_dir() and (args.output/"manifest.json").exists():
            (args.output/"FAILURE.json").write_text(json.dumps({"error": str(exc), "type": type(exc).__name__}, indent=2)+"\n")
        parser.exit(2, f"Run failed: {exc}\n")


if __name__ == "__main__":
    main()
