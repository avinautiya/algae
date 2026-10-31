"""Independent continuous-polytope vertex checks of the sharp formulas.

This validates implementation against a separate small optimization construction;
the analytic proof, not random testing, establishes the theorem.
"""
import argparse
import itertools
import json
import random
from pathlib import Path

from .core import ConductancePrior,ThermalMap,region_bounds


def solve_three(rows,values):
    a=[list(row)+[value] for row,value in zip(rows,values)]
    for j in range(3):
        pivot=max(range(j,3),key=lambda k:abs(a[k][j]))
        if abs(a[pivot][j])<1e-12:
            return None
        a[j],a[pivot]=a[pivot],a[j]
        scale=a[j][j]; a[j]=[x/scale for x in a[j]]
        for k in range(3):
            if k!=j:
                factor=a[k][j]; a[k]=[x-factor*y for x,y in zip(a[k],a[j])]
    return [a[i][3] for i in range(3)]


def vertices(lower,upper,total_lower,total_upper):
    # Enumerate intersections of three independent faces, rather than using
    # the monotonicity formulas or the production witness-construction routine.
    planes=[]
    for j in range(3):
        row=[float(i==j) for i in range(3)]
        planes.extend([(row,lower[j]),(row,upper[j])])
    planes.extend([([1,1,1],total_lower),([1,1,1],total_upper)])
    points=[]
    tol=1e-12*max(1,sum(upper))
    for selected in itertools.combinations(planes,3):
        point=solve_three([r for r,v in selected],[v for r,v in selected])
        if point is None:
            continue
        if all(lo-tol<=v<=hi+tol for lo,v,hi in zip(lower,point,upper)) and total_lower-tol<=sum(point)<=total_upper+tol:
            points.append(point)
    if not points:
        raise ArithmeticError("No feasible reference vertices")
    return points


def validate_randomized(cases=500,seed=3102026):
    if type(cases) is not int or cases<1:
        raise ValueError("cases must be a positive integer")
    rng=random.Random(seed); worst=0.0; rows=[]
    for i in range(cases):
        lower=[rng.uniform(100,10000) for _ in range(3)]
        upper=[l+rng.uniform(500,20000) for l in lower]
        center=[rng.uniform(l,u) for l,u in zip(lower,upper)]
        width=sum(u-l for l,u in zip(lower,upper))
        gl=max(sum(lower),sum(center)-rng.random()*width/3)
        gu=min(sum(upper),sum(center)+rng.random()*width/3)
        if i%3==0:
            gl=gu=sum(center)
        prior=ConductancePrior(tuple(lower),tuple(upper),gl,gu)
        model=ThermalMap(rng.uniform(.01,.2),rng.uniform(.01,.5))
        values=[model.region_power(p[0],p[1]) for p in vertices(lower,upper,gl,gu)]
        predicted=region_bounds(prior,1,2)
        errors=[]
        for label,reference in (("minimum",min(values)),("maximum",max(values))):
            w=predicted[label]; actual=model.region_power(w["prefix_W_K"],w["region_W_K"])
            errors.append(abs(actual-reference)/max(model.inlet_power_W,1.0))
        scaled=max(errors); worst=max(worst,scaled)
        if scaled>1e-10:
            raise AssertionError(f"Continuous vertex comparison failed in case {i}: {scaled}")
        rows.append(dict(case=i,scaled_extreme_error=scaled,reference_vertices=len(values)))
    return dict(status="PASS",cases=cases,seed=seed,max_power_scaled_error=worst,
        method="Independent enumeration of intersections of three polytope faces",
        scope="Numerical implementation check, not a substitute for proof or physical validation",rows=rows)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True); parser.add_argument("--cases",type=int,default=500)
    args=parser.parse_args()
    if args.cases<1:
        parser.error("cases must be positive")
    if args.output.exists():
        parser.error("Choose a new output file")
    result=validate_randomized(args.cases)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open("x") as stream:
        json.dump(result,stream,indent=2); stream.write("\n")
    print(json.dumps({k:v for k,v in result.items() if k!="rows"},indent=2))


if __name__=="__main__":
    main()
