"""Forward conditional incision envelopes, with no future observations as inputs."""
import csv
import hashlib
import json
from pathlib import Path

from .core import ConductancePrior,ThermalMap,finite,region_bounds


MODEL="Quasi-steady temperate-bed heat exchange with meltwater mass addition"
ASSUMPTIONS=[
    "Prescribed inlet scenarios; this model does not predict the inlet weather or discharge",
    "Frozen channel geometry and effective conductance shape over the forecast horizon",
    "Known common conductance multiplier may vary between intervals",
    "Temperate bed-only melting; all interface heat becomes latent heat",
    "No additional atmospheric, frictional or cold-ice heat terms and no retained thermal storage",
    "Exact internal prefix constraints are idealized information, not uncertain field observations",
    "Bounds are deterministic conditional envelopes, not statistical prediction intervals",
]


def keys(d,required,allowed,name):
    if not isinstance(d,dict) or not set(required).issubset(d) or set(d)-set(allowed):
        raise ValueError(f"Missing or unknown {name} fields")


def validate(case):
    required={"schema_version","name","cells","prior","forcing","regions","properties","scope_checks"}
    keys(case,required,required,"case")
    if case["schema_version"]!=1 or not isinstance(case["name"],str):
        raise ValueError("Expected schema_version=1 and name")
    if not case["cells"] or not case["forcing"]:
        raise ValueError("Need cells and forcing intervals")
    ck={"length_m","width_m","initial_bed_m","initial_depth_m","conductance_lower_W_K",
        "conductance_upper_W_K","conductance_nominal_W_K"}
    for c in case["cells"]:
        keys(c,ck,ck,"cell")
        for k in ("length_m","width_m","initial_depth_m"):
            finite(c[k],k,0,True)
        finite(c["initial_bed_m"],"bed")
        finite(c["conductance_nominal_W_K"],"nominal conductance",0)
    prior=case["prior"]
    pk={"total_lower_W_K","total_upper_W_K","exact_prefix_W_K","provenance"}
    keys(prior,pk,pk,"prior")
    if not isinstance(prior["provenance"],str) or not prior["provenance"].strip():
        raise ValueError("Declare prior provenance")
    p=ConductancePrior(tuple(c["conductance_lower_W_K"] for c in case["cells"]),
        tuple(c["conductance_upper_W_K"] for c in case["cells"]),prior["total_lower_W_K"],
        prior["total_upper_W_K"],tuple(prior["exact_prefix_W_K"]))
    nominal=[c["conductance_nominal_W_K"] for c in case["cells"]]
    if not p.feasible(nominal):
        raise ValueError("Nominal conductance profile violates the declared prior")
    propkeys={"rho_water","rho_ice","cp","latent_heat"}
    keys(case["properties"],propkeys,propkeys,"properties")
    for f in case["forcing"]:
        fk={"duration_s","inlet_Q_m3_s","inlet_excess_K","conductance_multiplier"}
        keys(f,fk,fk,"forcing")
        finite(f["duration_s"],"duration",0,True)
        finite(f["conductance_multiplier"],"conductance multiplier",0)
        if f["inlet_excess_K"]>4:
            raise ValueError("Forecast interface is restricted to 0..4 K excess; core mathematics has broader scope")
        ThermalMap(f["inlet_Q_m3_s"],f["inlet_excess_K"],**case["properties"])
    if not case["regions"]:
        raise ValueError("Declare at least one forecast region")
    names=[]
    for r in case["regions"]:
        keys(r,{"name","start","stop"},{"name","start","stop","critical_lowering_m"},"region")
        if not isinstance(r["name"],str) or not r["name"]:
            raise ValueError("Region name must be nonempty")
        names.append(r["name"]); region_bounds(p,r["start"],r["stop"])
        if "critical_lowering_m" in r:
            finite(r["critical_lowering_m"],"critical lowering",0,True)
    if len(names)!=len(set(names)):
        raise ValueError("Region names must be unique")
    sk={"max_fractional_depth_change","max_transit_to_interval_ratio"}
    keys(case["scope_checks"],sk,sk,"scope_checks")
    for k in sk:
        finite(case["scope_checks"][k],k,0,True)
    return p


def profile_powers(model,profile,multiplier):
    prefix=0.0; powers=[]
    for g in profile:
        powers.append(model.region_power(prefix*multiplier,g*multiplier))
        prefix+=g
    return powers


def forecast_case(case):
    prior=validate(case); n=len(case["cells"])
    nominal=[c["conductance_nominal_W_K"] for c in case["cells"]]
    areas=[c["width_m"]*c["length_m"] for c in case["cells"]]
    witnesses=[region_bounds(prior,i,i+1) for i in range(n)]
    regional={r["name"]:region_bounds(prior,r["start"],r["stop"]) for r in case["regions"]}
    heat_lo=[0.0]*n; heat_hi=[0.0]*n; heat_nom=[0.0]*n
    region_heat={r["name"]:dict(minimum=0.0,maximum=0.0,nominal=0.0) for r in case["regions"]}
    bed_rows=[]; region_rows=[]; outlet_rows=[]; time=0.0
    cumulative=dict(inlet_water_m3=0.0,outlet_water_m3=0.0,meltwater_m3=0.0,
                    inlet_energy_J=0.0,outlet_energy_J=0.0,latent_energy_J=0.0)
    all_bounds=region_bounds(prior,0,n)
    volume=sum(c["width_m"]*c["length_m"]*c["initial_depth_m"] for c in case["cells"])
    max_transit_ratio=0.0
    for interval,f in enumerate(case["forcing"]):
        model=ThermalMap(f["inlet_Q_m3_s"],f["inlet_excess_K"],**case["properties"])
        dt=f["duration_s"]; multiplier=f["conductance_multiplier"]; time+=dt
        max_transit_ratio=max(max_transit_ratio,volume/f["inlet_Q_m3_s"]/dt)
        powers=profile_powers(model,nominal,multiplier)
        for i,w in enumerate(witnesses):
            heat_nom[i]+=powers[i]*dt
            for label,acc in (("minimum",heat_lo),("maximum",heat_hi)):
                acc[i]+=model.region_power(w[label]["prefix_W_K"]*multiplier,w[label]["region_W_K"]*multiplier)*dt
            factor=model.rho_ice*model.latent_heat*areas[i]
            bed_rows.append(dict(time_s=time,cell=i,lowering_min_m=heat_lo[i]/factor,
                lowering_nominal_m=heat_nom[i]/factor,lowering_max_m=heat_hi[i]/factor,
                bed_nominal_m=case["cells"][i]["initial_bed_m"]-heat_nom[i]/factor))
        for r in case["regions"]:
            w=regional[r["name"]]; h=region_heat[r["name"]]
            for label in ("minimum","maximum"):
                h[label]+=model.region_power(w[label]["prefix_W_K"]*multiplier,w[label]["region_W_K"]*multiplier)*dt
            h["nominal"]+=sum(powers[r["start"]:r["stop"]])*dt
            area=sum(areas[r["start"]:r["stop"]]); scale=model.rho_ice*model.latent_heat*area
            region_rows.append(dict(time_s=time,region=r["name"],lowering_min_m=h["minimum"]/scale,
                lowering_nominal_m=h["nominal"]/scale,lowering_max_m=h["maximum"]/scale,
                lowering_minimax_m=(h["minimum"]+h["maximum"])/(2*scale),
                minimax_worst_case_error_m=(h["maximum"]-h["minimum"])/(2*scale),
                nominal_worst_case_error_m=max(h["nominal"]-h["minimum"],h["maximum"]-h["nominal"])/scale,
                regional_heat_min_J=h["minimum"],regional_heat_max_J=h["maximum"]))
        state=model.state(sum(nominal)*multiplier)
        max_total=all_bounds["maximum"]["total_W_K"]
        min_total=all_bounds["minimum"]["total_W_K"]
        cold=model.state(max_total*multiplier); warm=model.state(min_total*multiplier)
        cumulative["inlet_water_m3"]+=model.inlet_Q_m3_s*dt
        cumulative["outlet_water_m3"]+=state["Q_m3_s"]*dt
        cumulative["meltwater_m3"]+=state["removed_W"]*dt/(model.rho_water*model.latent_heat)
        cumulative["inlet_energy_J"]+=model.inlet_power_W*dt
        cumulative["outlet_energy_J"]+=state["power_W"]*dt
        cumulative["latent_energy_J"]+=state["removed_W"]*dt
        outlet_rows.append(dict(time_s=time,interval=interval,inlet_Q_m3_s=model.inlet_Q_m3_s,
            outlet_Q_nominal_m3_s=state["Q_m3_s"],outlet_temperature_nominal_K=state["temperature_excess_K"],
            outlet_temperature_min_K=cold["temperature_excess_K"],outlet_temperature_max_K=warm["temperature_excess_K"],
            total_heat_nominal_W=state["removed_W"],total_heat_min_W=warm["removed_W"],total_heat_max_W=cold["removed_W"],
            water_residual_m3=cumulative["outlet_water_m3"]-cumulative["inlet_water_m3"]-cumulative["meltwater_m3"],
            energy_residual_J=cumulative["inlet_energy_J"]-cumulative["outlet_energy_J"]-cumulative["latent_energy_J"]))
    rho_i=case["properties"]["rho_ice"]; latent=case["properties"]["latent_heat"]
    max_fraction=max(heat_hi[i]/(rho_i*latent*areas[i]*case["cells"][i]["initial_depth_m"]) for i in range(n))
    flags=[]
    if max_fraction>case["scope_checks"]["max_fractional_depth_change"]:
        flags.append("FROZEN_GEOMETRY_CHANGE_THRESHOLD_EXCEEDED")
    if max_transit_ratio>case["scope_checks"]["max_transit_to_interval_ratio"]:
        flags.append("QUASI_STEADY_TRANSIT_THRESHOLD_EXCEEDED")
    summary=dict(model=MODEL,name=case["name"],assumptions=ASSUMPTIONS,scope_flags=flags,
        scope_threshold_note="User-declared screening limits, not universal validity theorems",
        horizon_s=time,cumulative_nominal=cumulative,total_ice_nominal_m3=cumulative["latent_energy_J"]/(rho_i*latent),
        max_fractional_depth_change_bound=max_fraction,max_transit_to_interval_ratio=max_transit_ratio,
        max_abs_water_residual_m3=max(abs(r["water_residual_m3"]) for r in outlet_rows),
        max_abs_energy_residual_J=max(abs(r["energy_residual_J"]) for r in outlet_rows),
        marginal_bound_warning="Different cell/region bounds have different witnesses; never sum marginal maxima as a jointly realizable forecast")
    summary["minimax_note"]="Regional scalar minimax estimates minimize worst-case absolute error within the assumed class; they are not probability means or a jointly conservative spatial profile"
    def crossing(rows,column,threshold):
        last_t=last_value=0.0
        for row in rows:
            value=row[column]
            if value>=threshold:
                return last_t+(row["time_s"]-last_t)*(threshold-last_value)/(value-last_value)
            last_t,last_value=row["time_s"],value
        return None
    threshold_forecasts={}
    for region in case["regions"]:
        if "critical_lowering_m" not in region:
            continue
        rows=[r for r in region_rows if r["region"]==region["name"]]
        threshold=region["critical_lowering_m"]
        earliest=crossing(rows,"lowering_max_m",threshold)
        latest=crossing(rows,"lowering_min_m",threshold)
        threshold_forecasts[region["name"]]=dict(threshold_m=threshold,
            earliest_possible_s=earliest,latest_guaranteed_s=latest,
            nominal_crossing_s=crossing(rows,"lowering_nominal_m",threshold),
            status="NO_PROFILE_CROSSES_WITHIN_HORIZON" if earliest is None else
                   "ALL_PROFILES_CROSS_WITHIN_HORIZON" if latest is not None else "OCCURRENCE_UNRESOLVED_WITHIN_HORIZON",
            meaning="Declared mean-incision threshold only, not an independently established hydraulic-failure or flood threshold")
    summary["threshold_forecasts"]=threshold_forecasts
    return dict(summary=summary,bed=bed_rows,regions=region_rows,outlet=outlet_rows,
                regional_witnesses=regional,cell_witnesses=witnesses)


def write_forecast(case,output):
    output=Path(output)
    if output.exists():
        raise FileExistsError("Choose a new output directory")
    result=forecast_case(case); output.mkdir(parents=True)
    (output/"case.json").write_text(json.dumps(case,indent=2,allow_nan=False)+"\n")
    for name in ("summary","regional_witnesses","cell_witnesses"):
        (output/(name+".json")).write_text(json.dumps(result[name],indent=2,allow_nan=False)+"\n")
    for name in ("bed","regions","outlet"):
        with (output/(name+".csv")).open("w",newline="") as stream:
            writer=csv.DictWriter(stream,fieldnames=list(result[name][0]),lineterminator="\n")
            writer.writeheader(); writer.writerows(result[name])
    code=Path(__file__).resolve().parent
    manifest=dict(case_sha256=hashlib.sha256(json.dumps(case,sort_keys=True,separators=(",",":"),allow_nan=False).encode()).hexdigest(),
        source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(code.glob("*.py"))},
        scope="Conditional model forecast; all provided demonstration data are synthetic")
    (output/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    return result
