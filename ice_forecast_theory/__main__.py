"""Run offline forcing analysis or environmental screening; never launches CFD."""
import argparse
import hashlib
import json
from pathlib import Path

from .forcing import optimal_forcing
from .regime import regime_report


def save(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("x") as stream:
        json.dump(data,stream,indent=2,allow_nan=False); stream.write("\n")


def demos(output):
    output.mkdir(parents=True,exist_ok=False)
    forcing=dict(
        intervals=[dict(duration_s=1.0,A_per_s=[[-1.0]],B_per_s=[[1.0]]) for _ in range(2)],
        C=[[1.0]],budget_rows_on_input_amplitudes=[[1.0,1.0]],
        state_normalization="One abstract relaxation coordinate, not a glacier state",
        input_normalization="Abstract dimensionless input; time-average L2 norm",
        output_normalization="Same abstract state coordinate",
        operator_provenance="Manufactured relaxation equation dx/dt=-x+f, no physical claim")
    regime=dict(temperature_C=0.2,depth_m=0.05,velocity_m_s=1.0,length_m=5.0,width_m=0.5,
        ice_position="bed",potential_drop_m=0.1,external_power_W=0.0,
        local_patch=dict(depth_m=0.05,speed_m_s=0.005,length_m=0.2))
    save(output/"forcing_manufactured.json",forcing)
    save(output/"regime_synthetic.json",regime)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    subs=parser.add_subparsers(dest="command",required=True)
    demo=subs.add_parser("demo"); demo.add_argument("--output",type=Path,required=True)
    for name in ("forcing","regime"):
        sub=subs.add_parser(name); sub.add_argument("--input",type=Path,required=True)
        sub.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    try:
        if args.command=="demo":
            demos(args.output); print(json.dumps({"generated":str(args.output)})); return
        raw=args.input.read_bytes(); data=json.loads(raw)
        result=optimal_forcing(data) if args.command=="forcing" else regime_report(data)
        report=dict(input_sha256=hashlib.sha256(raw).hexdigest(),command=args.command,result=result,
            scope="Offline mathematical/physical screening, no solver integration or environmental forecast")
        save(args.output,report); print(json.dumps(report,indent=2,allow_nan=False))
    except (ValueError,KeyError,TypeError,IndexError,OverflowError,OSError) as exc:
        parser.exit(2,f"Analysis failed: {exc}\n")


if __name__=="__main__":
    main()
