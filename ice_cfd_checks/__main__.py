"""Generate known references or score exported SI-unit CSV/JSON data."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

from .diagnostics import convergence, integrate_faces, score_dam, score_jump, score_rest, score_thermal
from .references import DamBreak, StefanReference, rectangular_jump
from .linear_response import analyze


def write_json(path,data):
    path=Path(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("x") as stream:
        json.dump(data,stream,indent=2,allow_nan=False); stream.write("\n")


def write_csv(path,rows):
    rows=list(rows)
    with Path(path).open("x",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator="\n"); writer.writeheader(); writer.writerows(rows)


def generate(output):
    output.mkdir(parents=True,exist_ok=False)
    g=9.81; h1=0.1; width=0.5; froude=2.0
    q=width*h1*froude*(g*h1)**0.5
    jump=rectangular_jump(h1,q,width)
    example=dict(width_m=width,upstream=dict(depth_m=h1,discharge_m3_s=q),
                 downstream=dict(depth_m=jump["depth_down_m"],discharge_m3_s=q))
    write_json(output/"jump_reference.json",jump)
    write_json(output/"jump_exact_example.json",example)
    write_csv(output/"hydrostatic_exact_example.csv",[
        dict(x_m=x/10,bed_m=bed,depth_m=max(1-bed,0.0),velocity_m_s=0.0)
        for x,bed in enumerate((0.1,0.5,1.2,0.8,0.0))])
    for name,depth in (("wet",0.1),("dry",0.0)):
        ref=DamBreak(1.0,depth)
        rows=[]
        for i in range(200):
            a,b=-5+i*0.05,-5+(i+1)*0.05
            h,q=ref.cell_average(a,b,0.3)
            rows.append(dict(left_m=a,right_m=b,depth_mean_m=h,hu_mean_m2_s=q))
        write_csv(output/f"dam_break_{name}_exact_example.csv",rows)
    stefan=[]
    for ste in (0.001,0.003,0.013,0.2):
        ref=StefanReference(ste)
        for time in (1.0,60.0,3600.0):
            stefan.append(dict(stefan_sensible_over_latent=ste,time_s=time,
                wall_excess_K=ref.wall_excess_K,lambda_value=ref.lambda_value,
                front_m=ref.front(time),interface_flux_W_m2=ref.interface_flux(time),
                **ref.ledger_per_area(time)))
    write_csv(output/"stefan_equal_density_reference.csv",stefan)
    ledger=[dict(time_s=t,water_sensible_J=100+2*t,advective_in_J=10*t,
        advective_out_J=5*t,conductive_in_J=t,conductive_out_J=t,
        external_J=t,to_ice_J=4*t) for t in (0,1,2)]
    write_csv(output/"thermal_exact_example.csv",ledger)
    write_csv(output/"face_flux_example.csv",[
        dict(area_m2=1,liquid_fraction=1,normal_velocity_m_s=2,
             liquid_volume_flux_m3_s=2,temperature_excess_K=0.2),
        dict(area_m2=1,liquid_fraction=0.5,normal_velocity_m_s=-1,
             liquid_volume_flux_m3_s=-0.5,temperature_excess_K=0.6)])
    write_json(output/"convergence_exact_example.json",[
        dict(spacing_m=h,value=2+3*h*h) for h in (0.2,0.1,0.05)])
    write_json(output/"linear_response_manufactured.json",dict(
        A_per_s=[[-1,10],[0,-2]],times_s=[0,0.1,0.3,0.6931471805599453,1,2,5],
        state_normalization="Two abstract dimensionless coordinates; no glacier states",
        input_normalization="Euclidean unit perturbation in the two abstract coordinates",
        output_normalization="Euclidean norm in the same coordinates",
        operator_provenance="Manufactured stable non-normal matrix, not a CFD operator or scientific result"))
    write_json(output/"manifest.json",dict(
        scope="Known analytical references and manufactured examples, no CFD results",
        dam_break=dict(left_depth_m=1.0,time_s=0.3,gravity_m_s2=g,
                       right_depths_m=dict(wet=0.1,dry=0.0)),
        stefan="Equal-density one-phase conduction limit; densities are not independently varied",
        source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                       for p in sorted(Path(__file__).resolve().parent.glob("*.py"))}))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    subs=parser.add_subparsers(dest="command",required=True)
    ref=subs.add_parser("reference"); ref.add_argument("--output",type=Path,required=True)
    for name in ("score-jump","score-rest","score-dam","score-thermal","integrate-faces","convergence","linear-response"):
        sub=subs.add_parser(name)
        sub.add_argument("--input",type=Path,required=True)
        sub.add_argument("--output",type=Path,required=True)
        if name=="score-rest":
            sub.add_argument("--surface",type=float,required=True)
            sub.add_argument("--depth-scale",type=float,required=True)
        if name=="score-dam":
            sub.add_argument("--left-depth",type=float,required=True)
            sub.add_argument("--right-depth",type=float,required=True)
            sub.add_argument("--time",type=float,required=True)
    args=parser.parse_args()
    try:
        if args.command=="reference":
            generate(args.output)
            print(json.dumps({"generated":str(args.output),"scope":"Analytical/manufactured examples only"}))
            return
        if args.command in ("score-jump","convergence","linear-response"):
            data=json.loads(args.input.read_text())
        else:
            with args.input.open(newline="") as stream:
                data=list(csv.DictReader(stream))
        functions={"score-jump":lambda:score_jump(data),
            "score-rest":lambda:score_rest(data,args.surface,args.depth_scale),
            "score-dam":lambda:score_dam(data,DamBreak(args.left_depth,args.right_depth),args.time),
            "score-thermal":lambda:score_thermal(data),"integrate-faces":lambda:integrate_faces(data),
            "convergence":lambda:convergence(data),"linear-response":lambda:analyze(data)}
        result=functions[args.command]()
        report=dict(command=args.command,input_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),
                    diagnostics=result,acceptance="No automatic physical validation or novelty verdict")
        write_json(args.output,report)
        print(json.dumps(report,indent=2,allow_nan=False))
    except (ValueError,KeyError,TypeError,OSError,OverflowError) as exc:
        parser.exit(2,f"Diagnostic failed: {exc}\n")


if __name__=="__main__":
    main()
