"""Small analytical/numerical verification record, not a discovery campaign."""

import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import platform

import numpy as np
import scipy
from scipy.integrate import quad,solve_ivp

from .model import BaseConfig,HydraulicMap,SealModel,build_base,hydraulic_coefficients
from .matched import MatchedSealModel


def verify():
    # Exact uniform-coefficient shape derivative versus independent quadrature.
    q=.9; v=.7; alpha=beta=.5; fr=.2; length=2.
    d,hq,r,rv,rq=hydraulic_coefficients(v,q,alpha,beta,fr)
    k=float(rv/d); ell=(1-alpha)*v/((3-alpha)*q); z=float(hq-d*ell)
    bump=lambda x:x*x*(length-x)**2/length**4
    exact=-k/z*quad(lambda x:np.exp(-k*x)*bump(x),0,length,epsabs=1e-13)[0]
    quadrature=[]
    for n in [17,33,65,129,257]:
        x=np.linspace(0,length,n); mapping=HydraulicMap(x,np.full(n,v),q,alpha,beta,fr)
        j=mapping.solve(0,bump(x))[0]
        quadrature.append({"points":n,"spacing":length/(n-1),"prediction":j,
                           "independent_quadrature":exact,"absolute_error":abs(j-exact)})
    # Scalar response ODE uses exact dense base, not the discrete map recurrence.
    base=build_base(BaseConfig(left=-64,right=6,points=129))
    c=base.config
    def coefficients(x):
        velocity=base.evaluate([x])[0,0]
        d,hq,r,rv,rq=hydraulic_coefficients(velocity,1.,.5,.5,.3)
        return float(rv/d),float(rq+rv*hq/d)
    end_v=base.velocity[-1]
    d,hq,_,_,_=hydraulic_coefficients(end_v,1.,.5,.5,.3)
    terminal=float(hq-d*end_v/5)
    run=solve_ivp(lambda x,y:[coefficients(x)[0]*y[0]-coefficients(x)[1]],
                  (c.right,c.left),[terminal],dense_output=True,
                  rtol=2e-11,atol=2e-13,first_step=.01,max_step=.15)
    if not run.success: raise RuntimeError(run.message)
    # This is a tail-corrected finite reference, not a certified exact d_minus.
    reference=float(run.sol(-64)[0]+.5*64**-4)
    matching=[]
    for cutoff in [4,8,16,32,64]:
        value=float(run.sol(-cutoff)[0])
        matching.append({"left_cutoff":cutoff,"coefficient":value,
                         "tail_corrected_finite_reference":reference,
                         "reference_difference":reference-value,
                         "leading_asymptotic_term":.5*cutoff**-4})
    small_base=build_base(BaseConfig(points=65)); model=SealModel(small_base)
    forecast=model.forecast([(.5,.001),(.5,0)],sample_step=.1)
    matched=MatchedSealModel(small_base)
    kernel_times=np.linspace(0,2,21)
    matched_record={"kernel_scope":"Discretized outer-gauge renewal kernel C exp(L0 T) g_s; not an environmental timescale",
                    "kernel_times":kernel_times.tolist(),"renewal_kernel":matched.renewal_kernel(kernel_times).tolist(),
                    "diagnostics":matched.matching_diagnostics,
                    "forecast":matched.forecast([(.5,.001),(.5,0)],sample_step=.1)}
    source={str(p.name):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(Path(__file__).parent.glob('*.py'))}
    return {"scope":"Dimensionless equation/implementation verification only; no held-out physical skill test",
            "python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__,
            "source_sha256":source,"uniform_shape_quadrature":quadrature,
            "canonical_hydraulic_matching":matching,
            "matching_reference_scope":"Exact-base scalar ODE plus leading tail correction; not certified PDE error bounds",
            "forecast":forecast,"matched_outer_gauge":matched_record}


def main():
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument('--output',required=True)
    args=parser.parse_args(); target=Path(args.output)
    if target.exists(): raise FileExistsError(target)
    result=verify(); target.parent.mkdir(parents=True,exist_ok=True)
    with target.open('x') as f: json.dump(result,f,indent=2,allow_nan=False); f.write('\n')
    print(json.dumps({"record":str(target.resolve()),"scope":result["scope"],
                      "surrogate_head_residual":result['forecast']['maximum_total_head_ledger_residual'],
                      "surface_stage_residual":result['forecast']['maximum_surface_stage_storage_residual']}))


if __name__=='__main__': main()
