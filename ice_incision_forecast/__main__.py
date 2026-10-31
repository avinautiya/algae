import argparse
import json
from pathlib import Path

from .forecast import write_forecast


def demo_case():
    return dict(schema_version=1,name="synthetic_temperate_channel_1h",
        cells=[dict(length_m=1.0,width_m=1.0,initial_bed_m=-0.01*i,initial_depth_m=0.2,
            conductance_lower_W_K=200.0,conductance_upper_W_K=1800.0,conductance_nominal_W_K=1000.0) for i in range(20)],
        prior=dict(total_lower_W_K=20000.0,total_upper_W_K=20000.0,exact_prefix_W_K=[],
            provenance="Synthetic prior with exact total conductance; no future measurements used"),
        forcing=[dict(duration_s=1200.0,inlet_Q_m3_s=q,inlet_excess_K=t,conductance_multiplier=s)
                 for q,t,s in ((0.08,0.20,1.0),(0.12,0.30,1.1),(0.10,0.18,0.95))],
        regions=[dict(name="declared_control_region",start=8,stop=12,critical_lowering_m=0.002)],
        properties=dict(rho_water=1000.0,rho_ice=917.0,cp=4217.0,latent_heat=334000.0),
        scope_checks=dict(max_fractional_depth_change=0.05,max_transit_to_interval_ratio=0.1))


def main():
    parser=argparse.ArgumentParser(description="Forecast conditional incision bounds without future outlet observations")
    subs=parser.add_subparsers(dest="command",required=True)
    demo=subs.add_parser("demo"); demo.add_argument("--output",type=Path,required=True)
    predict=subs.add_parser("predict"); predict.add_argument("--case",type=Path,required=True)
    predict.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    try:
        if args.command=="demo":
            with args.output.open("x") as stream:
                json.dump(demo_case(),stream,indent=2); stream.write("\n")
            print(str(args.output)); return
        result=write_forecast(json.loads(args.case.read_text()),args.output)
        print(json.dumps(result["summary"],indent=2,allow_nan=False))
    except (ValueError,KeyError,TypeError,OSError,ArithmeticError) as exc:
        parser.exit(2,f"Forecast failed: {exc}\n")


if __name__=="__main__":
    main()
