"""Auditable dimensionless linear-forecast CLI; output paths are preserved."""

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import numpy as np

from .model import BaseConfig, SealModel, build_base


def _write_new(path, data):
    p=Path(path)
    p.parent.mkdir(parents=True,exist_ok=True)
    with p.open("x") as f: json.dump(data,f,indent=2,allow_nan=False); f.write("\n")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    subs=parser.add_subparsers(dest="command",required=True)
    demo=subs.add_parser("demo-case"); demo.add_argument("--output",required=True)
    predict=subs.add_parser("predict"); predict.add_argument("--case",required=True)
    predict.add_argument("--output",required=True)
    predict.add_argument("--engine",choices=["linear","nonlinear","matched"],default="linear")
    args=parser.parse_args()
    if args.command=="demo-case":
        case={"purpose":"Manufactured dimensionless verification example; not a glacier forecast",
              "base":asdict(BaseConfig(points=65)),"lake_capacity":2.0,
              "initial_head":0.0,"initial_bed":None,"sample_step":0.1,
              "relative_limit":0.1,"segments":[[0.5,0.001],[0.5,0.0]]}
        _write_new(args.output,case)
        print(json.dumps({"case":str(Path(args.output).resolve()),"scope":case["purpose"]}))
        return
    with Path(args.case).open() as f: case=json.load(f)
    output=Path(args.output)
    if output.exists(): raise FileExistsError("Forecast output exists; choose a new path")
    base=build_base(BaseConfig(**case["base"]))
    model=SealModel(base,case.get("lake_capacity",2.0))
    bed=case.get("initial_bed")
    initial=model.initial_state(case.get("initial_head",0.),bed)
    if args.engine=="linear":
        result=model.forecast(case["segments"],initial,case.get("sample_step",0.1),
                              case.get("relative_limit",0.1))
    elif args.engine=="nonlinear":
        from .nonlinear import NonlinearSealModel
        nonlinear=NonlinearSealModel(base,case.get("lake_capacity",2.))
        result=nonlinear.forecast(case["segments"],initial,case.get("maximum_step",.02),
                                  case.get("relative_limit",.1))
    else:
        from .matched import MatchedSealModel
        matched=MatchedSealModel(base,case.get("lake_capacity",2.))
        start=matched.initial_state(case.get("initial_outer_relative_head",case.get("initial_head",0.)),
                                    case.get("initial_shape"))
        result=matched.forecast(case["segments"],start,case.get("sample_step",.1),
                                case.get("relative_limit",.1))
    result["input_case"]=case
    result["forecast_origin_rule"]="All segment inputs are supplied scenarios; no future measured output is used"
    _write_new(output,result)
    print(json.dumps({"output":str(output.resolve()),
                      "maximum_total_head_ledger_residual":result.get("maximum_total_head_ledger_residual"),
                      "maximum_matched_limit_storage_ledger_residual":result.get("maximum_matched_limit_storage_ledger_residual"),
                      "maximum_surface_stage_storage_residual":result.get("maximum_surface_stage_storage_residual"),
                      "all_samples_within_declared_linear_screen":result.get("all_samples_within_declared_linear_screen"),
                      "environmental_validation":False}))


if __name__=="__main__": main()
