"""Reproducible synthetic study; protocol is committed before production runs."""
import argparse
import copy
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
from .channel import Channel,CrossStream,plug_cell_bulk_kernel,modal_witness
from .forecast import demo_case,dump,write_forecast
from .recovery import sharp_interval


def enumerate_vertices(c,e,d,eps,lo,hi):
    """Independent tiny LP reference: enumerate active faces, solve linear systems."""
    n=len(c); a=np.vstack([np.eye(n),-np.eye(n),e,-e]); b=np.r_[hi,-lo,d+eps,-d+eps]
    values=[]
    for ids in itertools.combinations(range(len(b)),n):
        aa=a[list(ids)]
        if abs(np.linalg.det(aa))<1e-12:
            continue
        x=np.linalg.solve(aa,b[list(ids)])
        if np.max(a@x-b)<1e-9:
            values.append(float(c@x))
    if not values:
        raise ArithmeticError("Independent vertex enumeration found no feasible point")
    return min(values),max(values)


def interval(model,stations,region,profile,cap,noise=0):
    e=np.array([model.bulk_kernel(x) for x in stations])
    c=model.heat_kernel_W_K(*region)
    return sharp_interval(c,e,e@profile,np.full(len(stations),noise),
                          np.zeros(model.cells),np.full(model.cells,cap))


def run(output):
    out=Path(output)
    if out.exists():
        raise FileExistsError("Choose a new study directory")
    protocol_path=Path(__file__).with_name("protocol.json")
    p=json.loads(protocol_path.read_text())
    out.mkdir(parents=True)
    dump(out/"protocol.json",p)
    dump(out/"source_snapshot.json",{f.name:hashlib.sha256(f.read_bytes()).hexdigest()
                                      for f in Path(__file__).parent.glob("*.py")})
    rng=np.random.default_rng(p["seed"])
    cap=p["temperature_cap_K"];stations=p["stations_m"];region=p["target_m"]
    ch=Channel();n=p["primary_cells"];m=CrossStream(ch,n)
    reference=np.full(n,p["reference_inlet_K"])
    # Forecast and source extension: actual CLI contract, no future observations.
    base_case=demo_case()
    base=write_forecast(base_case,out/"forecast")
    source_case=copy.deepcopy(base_case)
    source_case["surface_heat_segments"]=[dict(start_m=.04,stop_m=.12,lower_W_m2=0.,upper_W_m2=100.,
                                                    provenance="Synthetic declared nonnegative net surface-heating interval")]
    source=write_forecast(source_case,out/"surface_heating_forecast")
    print("Forecast and uncertain-surface-heating extension complete",flush=True)
    # An analytical continuum witness is independent of FV and optimization.
    modal=modal_witness(Channel(velocity="plug"),stations,*region,cap)
    dump(out/"continuum_witness.json",modal)
    refinement=[]
    for velocity in ("plug","shear_free"):
        for cells in p["transverse_refinement"]:
            mm=CrossStream(Channel(velocity=velocity),cells)
            f=np.full(cells,p["reference_inlet_K"])
            r=interval(mm,stations,region,f,cap)
            row=dict(velocity=velocity,cells=cells,lower_W=r["lower"],upper_W=r["upper"],
                     minimax_W=r["minimax"],max_dual_gap_W=max(abs(r["upper_certificate"]["dual_gap"]),
                                                               abs(r["negative_lower_certificate"]["dual_gap"])))
            if velocity=="plug":
                exact_a=plug_cell_bulk_kernel(mm.channel,mm.edges,region[0],1024)
                exact_b=plug_cell_bulk_kernel(mm.channel,mm.edges,region[1],1024)
                analytic=mm.channel.capacity_W_K*float((exact_a-exact_b)@f)
                fv=float(mm.heat_kernel_W_K(*region)@f)
                row.update(analytic_uniform_power_W=analytic,fv_uniform_power_W=fv,
                           relative_uniform_power_error=abs(fv-analytic)/analytic)
            refinement.append(row)
    dump(out/"refinement.json",refinement)
    print("Plug analytical and shear-flow refinement complete",flush=True)
    # Independent active-face enumeration, including noisy information.
    checks=[]
    for i in range(p["independent_lp_vertex_cases"]):
        size=3 if i%2 else 4; count=1 if size==3 else 2
        e=rng.uniform(0,1,(count,size));e/=e.sum(axis=1)[:,None]
        f=rng.uniform(.02,.18,size);d=e@f;eps=np.full(count,.001 if i%3 else 0)
        c=rng.normal(size=size);lo=np.zeros(size);hi=np.full(size,cap)
        r=sharp_interval(c,e,d,eps,lo,hi);a,b=enumerate_vertices(c,e,d,eps,lo,hi)
        error=max(abs(a-r["lower"]),abs(b-r["upper"]))
        if error>p["acceptance"]["independent_vertex_absolute_error"]:
            raise ArithmeticError("Independent LP reference disagreement")
        checks.append(dict(case=i,max_error=error,max_dual_gap=max(abs(r["upper_certificate"]["dual_gap"]),
                                                                abs(r["negative_lower_certificate"]["dual_gap"]))))
    dump(out/"independent_vertex_checks.json",checks)
    print("Independent continuous LP checks complete",flush=True)
    # Same information for all baselines; no training or test-set fitting.
    e=np.array([m.bulk_kernel(x) for x in stations]);c=m.heat_kernel_W_K(*region)
    uniform=np.ones(n);uniform_data=e@uniform;uniform_power=float(c@uniform)
    basis=m.eigenvectors[:,:2]/m.root_mass[:,None]
    basis=basis/np.max(np.abs(basis),axis=0)[None,:]
    two_mode_observation=e@basis
    cases=[]
    for i in range(p["smooth_profile_cases"]):
        eta=m.y/ch.depth_m;offset=rng.uniform(-.5,.5);amplitudes=rng.normal(size=5)
        signal=np.full(n,offset)
        for j,a in enumerate(amplitudes,1):
            signal+=a*np.cos(j*np.pi*eta)/j
        f=cap*np.sin(np.pi*eta/2)/(1+np.exp(-signal))
        d=e@f;truth=float(c@f)
        r=sharp_interval(c,e,d,np.zeros(len(stations)),np.zeros(n),np.full(n,cap))
        coefficient=float(np.clip(uniform_data@d/(uniform_data@uniform_data),0,cap))
        modal_coefficients=np.linalg.solve(two_mode_observation,d)
        modal_profile=basis@modal_coefficients
        predictions=dict(uniform_entrance=coefficient*uniform_power,
                         two_mode_graetz=float(c@modal_profile),minimax=r["minimax"])
        if not r["lower"]-1e-9<=truth<=r["upper"]+1e-9:
            raise ArithmeticError("Admitted true profile outside forecast interval")
        cases.append(dict(case=i,offset=float(offset),amplitudes=amplitudes.tolist(),
                          observations_K=d.tolist(),true_power_W=truth,
                          interval_W=[r["lower"],r["upper"]],predictions_W=predictions,
                          signed_errors_W={k:v-truth for k,v in predictions.items()},
                          uniform_observation_max_residual_K=float(np.max(np.abs(coefficient*uniform_data-d))),
                          two_mode_initial_bound_violation_K=float(max(0,-modal_profile.min(),modal_profile.max()-cap))))
        if (i+1)%30==0:
            dump(out/"smooth_profile_cases.json",cases)
            print("Smooth-profile comparison",i+1,"/",p["smooth_profile_cases"],flush=True)
    dump(out/"smooth_profile_cases.json",cases)
    metrics={}
    for name in ("uniform_entrance","two_mode_graetz","minimax"):
        errors=np.array([v["signed_errors_W"][name] for v in cases])
        metrics[name]=dict(mean_abs_error_W=float(np.mean(np.abs(errors))),
                           median_abs_error_W=float(np.median(np.abs(errors))),
                           max_abs_error_W=float(np.max(np.abs(errors))),mean_signed_error_W=float(errors.mean()))
    # Noise and station information, fixed reference: monotonicity is expected.
    noise=[]
    for count,ids in ((1,[1]),(2,[0,1])):
        for sigma in p["measurement_noise_K"]:
            ss=[stations[j] for j in ids]
            r=interval(m,ss,region,reference,cap,sigma)
            noise.append(dict(stations=count,noise_bound_K=sigma,lower_W=r["lower"],upper_W=r["upper"]))
    dump(out/"measurement_noise.json",noise)
    # Thermal-age curve with nested past station sets.
    age=[]
    length=ch.speed_m_s*ch.depth_m**2/ch.diffusivity
    for beta in np.geomspace(.005,.16,13):
        probe=2*beta*length;target=[3*beta*length,6*beta*length]
        all_stations=np.linspace(probe/4,probe,4)
        for count,ids in ((1,[3]),(2,[1,3]),(4,[0,1,2,3])):
            ss=all_stations[ids]
            r=interval(m,ss,target,reference,cap)
            age.append(dict(beta=float(beta),stations=count,target_start_m=float(target[0]),
                            target_stop_m=float(target[1]),lower_W=r["lower"],upper_W=r["upper"],
                            half_width_over_midpoint=float((r["upper"]-r["lower"])/(r["upper"]+r["lower"]))))
    dump(out/"thermal_age.json",age)
    base_r=base["regions"][0];source_r=source["regions"][0]
    certs=[base_r["power_interval"]["upper_certificate"],base_r["power_interval"]["negative_lower_certificate"]]
    summary=dict(status="PASS",protocol_sha256=hashlib.sha256(protocol_path.read_bytes()).hexdigest(),
                 scope="Synthetic conditional model study; not new CFD or field validation",test_profiles=len(cases),
                 independently_enumerated_cases=len(checks),max_vertex_error=max(v["max_error"] for v in checks),
                 continuum_witness_power_difference_W=modal["regional_power_difference_W"],
                 continuum_witness_max_observation_difference_K=max(abs(v) for v in modal["calibration_difference_K"]),
                 primary_heat_power_W={k:base_r["power_interval"][k] for k in ("lower","upper","minimax","worst_case_absolute_error")},
                 primary_lowering_m=base_r["lowering_m"],
                 upper_over_lower=base_r["power_interval"]["upper"]/base_r["power_interval"]["lower"],
                 max_primary_dual_gap_W=max(abs(v["dual_gap"]) for v in certs),
                 source_extended_heat_power_W={k:source_r["power_interval"][k] for k in ("lower","upper","minimax")},
                 smooth_profile_metrics=metrics,
                 smooth_cases_with_infeasible_two_mode_initial_profile=sum(v["two_mode_initial_bound_violation_K"]>1e-9 for v in cases),
                 smooth_cases_consistent_with_uniform_entrance_at_1e_6K=sum(v["uniform_observation_max_residual_K"]<=1e-6 for v in cases),
                 all_admitted_profiles_inside_conditional_intervals=True,
                 no_statistical_coverage_claim=True)
    dump(out/"summary.json",summary)
    print(json.dumps(summary,indent=2),flush=True)
    return summary


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--output",required=True)
    run(parser.parse_args().output)
