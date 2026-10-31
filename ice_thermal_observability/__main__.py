import argparse
import json
from pathlib import Path
from .forecast import demo_case,dump,write_forecast


def main():
    parser=argparse.ArgumentParser(description="Conditional regional-incision forecast from incomplete thermal information")
    sub=parser.add_subparsers(dest="command",required=True)
    d=sub.add_parser("demo");d.add_argument("--output",required=True)
    p=sub.add_parser("predict");p.add_argument("--case",required=True);p.add_argument("--output",required=True)
    args=parser.parse_args()
    if args.command=="demo":
        if Path(args.output).exists():
            raise FileExistsError("Choose a new case output")
        dump(args.output,demo_case())
        print(args.output)
    else:
        result=write_forecast(json.loads(Path(args.case).read_text()),args.output)
        print(json.dumps({"output":args.output,"scope_flags":result["scope_flags"],
                          "regions":[dict(name=r["region"]["name"],lowering_m=r["lowering_m"]) for r in result["regions"]]},indent=2))


if __name__=="__main__":
    main()
