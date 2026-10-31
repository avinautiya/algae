"""Forward steady-scenario forecast; only upstream/current observations allowed."""
from dataclasses import asdict
from pathlib import Path
import csv
import hashlib
import json
import platform
import numpy as np
import scipy
from .channel import Channel, CrossStream, positive
from .recovery import sharp_interval, identifiable


def dump(path,data):
    Path(path).write_text(json.dumps(data,indent=2,allow_nan=False)+"\n")


def demo_case():
    ch = Channel()
    m = CrossStream(ch,96)
    station = [0.02,0.04]
    return dict(schema="thermal-observability-0.1",label="SYNTHETIC METHOD DEMONSTRATION",
                channel=asdict(ch),cells=96,stretch=1.0,forecast_origin_m=0.04,
                horizon_s=60.0,initial_lower_K=0.0,initial_upper_K=0.2,
                observations=[dict(x_m=x,bulk_excess_K=m.bulk(np.full(m.cells,0.1),x),
                                   error_bound_K=0.0,provenance="Synthetic uniform 0.1 K inlet; exact for declared FV grid") for x in station],
                surface_heat_segments=[],
                regions=[dict(name="downstream_control",start_m=0.06,stop_m=0.12)],
                scope=dict(max_declared_lowering_fraction_of_depth=0.05,
                           provenance="User screening limit; not a proved applicability threshold"))


def build(case):
    allowed={"schema","label","channel","cells","stretch","forecast_origin_m","horizon_s",
             "initial_lower_K","initial_upper_K","observations","surface_heat_segments","regions","scope"}
    if set(case)-allowed or case.get("schema")!="thermal-observability-0.1":
        raise ValueError("Unknown fields/schema; future observations are not accepted")
    ch=Channel(**case["channel"])
    m=CrossStream(ch,case["cells"],case.get("stretch",1.0))
    origin=positive(case["forecast_origin_m"],"forecast origin",True)
    positive(case["horizon_s"],"horizon")
    lo=np.broadcast_to(np.asarray(case["initial_lower_K"],float),(m.cells,)).copy()
    hi=np.broadcast_to(np.asarray(case["initial_upper_K"],float),(m.cells,)).copy()
    if not np.isfinite(lo).all() or not np.isfinite(hi).all() or np.any(lo<0) or np.any(lo>hi) or np.max(hi)>4:
        raise ValueError("Require 0 <= initial lower <= upper <= 4 K")
    observations=case["observations"]
    if not observations:
        raise ValueError("Supply at least one upstream/current observation")
    positions=[]
    for o in observations:
        if set(o)!={"x_m","bulk_excess_K","error_bound_K","provenance"}:
            raise ValueError("Unknown observation fields")
        x=positive(o["x_m"],"observation position",True)
        if x>origin:
            raise ValueError("Future/downstream observations are excluded from a forward forecast")
        positive(o["bulk_excess_K"],"observed excess",True)
        positive(o["error_bound_K"],"measurement error",True)
        if not isinstance(o["provenance"],str) or not o["provenance"].strip():
            raise ValueError("Observation provenance required")
        positions.append(x)
    if len(set(positions))!=len(positions):
        raise ValueError("Duplicate observation positions")
    sources=case.get("surface_heat_segments",[])
    for s in sources:
        if set(s)!={"start_m","stop_m","lower_W_m2","upper_W_m2","provenance"}:
            raise ValueError("Unknown surface-source fields")
        if not origin<=s["start_m"]<s["stop_m"]:
            raise ValueError("Forecast sources must lie at/downstream of origin")
        positive(s["lower_W_m2"],"surface lower",True)
        positive(s["upper_W_m2"],"surface upper",True)
        if s["lower_W_m2"]>s["upper_W_m2"] or not str(s["provenance"]).strip():
            raise ValueError("Invalid surface-source interval/provenance")
    if any(sources[i]["stop_m"]>sources[i+1]["start_m"] for i in range(len(sources)-1)):
        raise ValueError("Surface segments must be ordered and disjoint")
    for r in case["regions"]:
        if set(r)!={"name","start_m","stop_m"} or not origin<=r["start_m"]<r["stop_m"]:
            raise ValueError("Named forecast regions must lie at/downstream of origin")
    if not case["regions"] or len({r["name"] for r in case["regions"]})!=len(case["regions"]):
        raise ValueError("Nonempty unique regions required")
    e=np.array([np.r_[m.bulk_kernel(x),np.zeros(len(sources))] for x in positions])
    d=np.array([o["bulk_excess_K"] for o in observations])
    eps=np.array([o["error_bound_K"] for o in observations])
    return m,e,d,eps,np.r_[lo,[s["lower_W_m2"] for s in sources]],np.r_[hi,[s["upper_W_m2"] for s in sources]]


def witness_audit(m,vector,case,region):
    f=np.asarray(vector[:m.cells]); source_values=vector[m.cells:]
    def state(x):
        t=m.propagate(f,x)
        for s,q in zip(case.get("surface_heat_segments",[]),source_values):
            t=t+q*m.surface_response(s["start_m"],s["stop_m"],x)
        return t
    a,b=region["start_m"],region["stop_m"]
    ta,tb=state(a),state(b)
    pa=float(m.mass@ta*m.channel.width_m); pb=float(m.mass@tb*m.channel.width_m)
    supplied=sum(q*m.channel.width_m*max(0,min(b,s["stop_m"])-max(a,s["start_m"]))
                 for s,q in zip(case.get("surface_heat_segments",[]),source_values))
    heat=pa-pb+supplied
    rho_i,latent=m.channel.rho_ice,m.channel.latent_heat
    return dict(inlet_region_sensible_W=pa,outlet_region_sensible_W=pb,
                supplied_surface_W=float(supplied),bed_heat_W=float(heat),
                latent_ice_volume_m3_s=float(heat/(rho_i*latent)),
                bed_lowering_m=float(heat*m.lowering_factor(a,b,case["horizon_s"])),
                min_reported_temperature_K=float(min(ta.min(),tb.min())),
                max_reported_temperature_K=float(max(ta.max(),tb.max())),
                initial_profile_K=f.tolist(),surface_flux_W_m2=list(source_values),
                meaning="Fixed-flow thermal-to-latent diagnostic; no moving boundary or meltwater feedback")


def predict(case):
    m,e,d,eps,lo,hi=build(case)
    sources=case.get("surface_heat_segments",[])
    results=[]
    for region in case["regions"]:
        a,b=region["start_m"],region["stop_m"]
        objective=np.r_[m.heat_kernel_W_K(a,b),
                        [m.surface_heat_coefficient(s["start_m"],s["stop_m"],a,b) for s in sources]]
        interval=sharp_interval(objective,e,d,eps,lo,hi)
        factor=m.lowering_factor(a,b,case["horizon_s"])
        item=dict(region=region,power_interval=interval,
                  lowering_m={key:float(interval[key]*factor) for key in ("lower","upper","minimax","worst_case_absolute_error")},
                  linear_information=identifiable(objective,e),
                  lower_witness=witness_audit(m,interval["lower_witness_variables"],case,region),
                  upper_witness=witness_audit(m,interval["upper_witness_variables"],case,region))
        for name,expected in (("lower_witness",interval["lower"]),("upper_witness",interval["upper"])):
            if abs(item[name]["bed_heat_W"]-expected)>1e-9*max(abs(expected),1):
                raise ArithmeticError("Conservative replay disagrees with LP functional")
        results.append(item)
    max_fraction=max(r["lowering_m"]["upper"] for r in results)/m.channel.depth_m
    limit=case.get("scope",{}).get("max_declared_lowering_fraction_of_depth",0.05)
    positive(limit,"geometry screening limit")
    flags=[]
    if max_fraction>limit:
        flags.append("Declared geometry-change screen exceeded; fixed geometry needs re-evaluation")
    return dict(schema=case["schema"],label=case.get("label","USER SCENARIO"),
                model="Fixed-velocity transverse advection-diffusion; spatially steady scenario held over declared horizon",
                dimensionless=m.channel.dimensionless(float(np.max(hi[:m.cells]))),
                regions=results,scope_flags=flags,max_region_mean_lowering_fraction_of_depth=max_fraction,
                scope=["Sharp bounds are for the declared finite-volume model and profile/source class",
                       "No certified continuum discretization or environmental-model error enclosure",
                       "No evolving hydraulics, free surface, turbulence, buoyancy, ice conduction or meltwater mass feedback",
                       "Sources are prescribed nonnegative net surface heat; evaporation/cooling/freezing are not included",
                       "Steady forcing persistence is an input assumption, not a weather forecast",
                       "Intervals are deterministic feasible ranges, not confidence intervals"])


def write_forecast(case,output):
    path=Path(output)
    if path.exists():
        raise FileExistsError("Choose a new output directory")
    result=predict(case)
    path.mkdir(parents=True)
    dump(path/"case.json",case)
    dump(path/"forecast.json",result)
    rows=[]
    for r in result["regions"]:
        rows.append(dict(region=r["region"]["name"],**{k+"_m":v for k,v in r["lowering_m"].items()}))
    with (path/"regions.csv").open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    sources={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in Path(__file__).parent.glob("*.py")}
    dump(path/"manifest.json",dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
         case_sha256=hashlib.sha256((path/"case.json").read_bytes()).hexdigest(),source_sha256=sources))
    return result
