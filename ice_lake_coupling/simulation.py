"""Executable synthetic lake/control comparison, deliberately not a CFD solver."""
import csv
import hashlib
import json
from dataclasses import asdict
from pathlib import Path

from .kernel import (
    Constants, DomainError, LakeForcing, LakeState, StageAreaTable,
    advance_lake, critical_discharge, finite, thermal_control,
)

MODEL = "Well-mixed lake + ideal critical outlet + prescribed thermal exchange"
LIMITS = [
    "Synthetic scenario, not a calibrated or validated glacier forecast",
    "No resolved CFD, channel storage, backwater or hydrofracture",
    "No debris, stratification, lake freezing, glacier mechanics or inundation routing",
    "Heat-transfer coefficient and contact area are prescribed, not learned or resolved",
    "Fixed lake hypsometry; channel meltwater leaves downstream, not back into lake",
    "First-order time integration and explicit sill update require refinement",
]


def require_keys(data, allowed, required, name):
    if not isinstance(data, dict):
        raise DomainError(f"{name} must be an object")
    extra, missing = set(data) - set(allowed), set(required) - set(data)
    if extra or missing:
        raise DomainError(f"{name}: unknown fields {sorted(extra)}, missing fields {sorted(missing)}")


def validate_case(case):
    keys = {"schema_version", "name", "duration_s", "max_dt_s", "max_lowering_per_step_m",
            "storage", "initial", "control", "forcing"}
    require_keys(case, keys, keys, "case")
    if case["schema_version"] != 1 or not isinstance(case["name"], str):
        raise DomainError("Expected schema_version=1 and a string name")
    for key in ("duration_s", "max_dt_s", "max_lowering_per_step_m"):
        finite(case[key], key, 0.0, strictly=True)
    require_keys(case["storage"], {"stages_m", "areas_m2"}, {"stages_m", "areas_m2"}, "storage")
    storage = StageAreaTable(**case["storage"])
    require_keys(case["initial"], {"stage_m", "temperature_K"}, {"stage_m", "temperature_K"}, "initial")
    storage.volume(case["initial"]["stage_m"])
    finite(case["initial"]["temperature_K"], "initial temperature", 0.0)
    ck = {"sill_m", "width_m", "contact_length_m", "htc_W_m2_K", "discharge_coefficient"}
    require_keys(case["control"], ck, ck, "control")
    c = case["control"]
    finite(c["contact_length_m"], "contact length", 0.0, strictly=True)
    finite(c["htc_W_m2_K"], "htc", 0.0)
    critical_discharge(case["initial"]["stage_m"], c["sill_m"], c["width_m"], c["discharge_coefficient"])
    if c["sill_m"] < storage.stages_m[0]:
        raise DomainError("Control below basin floor is not supported")
    if not isinstance(case["forcing"], list) or not case["forcing"]:
        raise DomainError("Forcing must contain one or more piecewise-constant intervals")
    previous = 0.0
    for f in case["forcing"]:
        fk = {"until_s", "inflow_m3_s", "inflow_temperature_K", "external_power_W"}
        require_keys(f, fk, fk, "forcing interval")
        finite(f["until_s"], "forcing interval end", previous, strictly=True)
        LakeForcing(**{k: v for k, v in f.items() if k != "until_s"})
        previous = f["until_s"]
    if previous != case["duration_s"]:
        raise DomainError("Final forcing interval must end exactly at duration_s")
    return storage


def run_case(case, melting, constants=Constants()):
    storage = validate_case(case)
    if not isinstance(melting, bool):
        raise DomainError("melting must be boolean")
    state = LakeState.from_temperature(storage.volume(case["initial"]["stage_m"]),
                                       case["initial"]["temperature_K"], constants)
    initial = state
    c = case["control"]
    sill = initial_sill = c["sill_m"]
    area = c["width_m"] * c["contact_length_m"]
    totals = dict(inflow_m3=0.0, inflow_J=0.0, external_J=0.0, port_m3=0.0,
                  port_J=0.0, delivered_m3=0.0, exported_J=0.0, latent_J=0.0,
                  meltwater_m3=0.0, ice_loss_m3=0.0)
    rows, time, segment, rejected = [], 0.0, 0, 0
    peak_q, peak_time = 0.0, 0.0
    def snapshot(dt, qlake, qout, tout, latent_power):
        residual_v = state.volume_m3 + totals["delivered_m3"] - initial.volume_m3 - totals["inflow_m3"] - totals["meltwater_m3"]
        residual_e = state.sensible_energy_J + totals["exported_J"] + totals["latent_J"] - initial.sensible_energy_J - totals["inflow_J"] - totals["external_J"]
        return dict(time_s=time, dt_s=dt, lake_stage_m=storage.stage(state.volume_m3),
                    lake_volume_m3=state.volume_m3, lake_temperature_K=state.temperature(constants),
                    sill_m=sill, lake_to_control_m3_s=qlake, delivered_m3_s=qout,
                    outlet_temperature_K=tout, latent_power_W=latent_power,
                    cumulative_delivered_m3=totals["delivered_m3"],
                    cumulative_exported_J=totals["exported_J"], cumulative_latent_J=totals["latent_J"],
                    water_residual_m3=residual_v, energy_residual_J=residual_e)
    q0 = critical_discharge(storage.stage(state.volume_m3), sill, c["width_m"], c["discharge_coefficient"], constants)
    peak_q, peak_time = q0, 0.0
    rows.append(snapshot(0.0, q0, None, None, None))
    while time < case["duration_s"]:
        f = case["forcing"][segment]
        forcing = LakeForcing(**{k: v for k, v in f.items() if k != "until_s"})
        endpoint = min(time + case["max_dt_s"], f["until_s"], case["duration_s"])
        dt = endpoint - time
        if dt <= 0.0:
            raise DomainError("Nonpositive step; forcing interval bookkeeping error")
        while True:
            lake = advance_lake(storage, state, forcing, dt,
                lambda h: critical_discharge(h, sill, c["width_m"], c["discharge_coefficient"], constants), constants)
            heat = thermal_control(lake.port, dt, area, c["htc_W_m2_K"], melting, constants)
            if heat.lowering_m <= case["max_lowering_per_step_m"]:
                break
            dt *= 0.5
            endpoint = time + dt
            rejected += 1
            if dt < 1e-6:
                raise DomainError("Geometry-step restriction requires dt below 1 microsecond")
        if sill - heat.lowering_m < storage.stages_m[0]:
            raise DomainError("Sill reaches lake floor; this demonstration cannot continue")
        state = lake.state
        sill -= heat.lowering_m
        time = endpoint
        for key, value in (
            ("inflow_m3", lake.inflow_volume_m3), ("inflow_J", lake.inflow_energy_J),
            ("external_J", lake.external_energy_J), ("port_m3", lake.port.volume_m3),
            ("port_J", lake.port.sensible_energy_J), ("delivered_m3", heat.outflow.volume_m3),
            ("exported_J", heat.outflow.sensible_energy_J), ("latent_J", heat.latent_energy_J),
            ("meltwater_m3", heat.meltwater_m3), ("ice_loss_m3", heat.ice_loss_m3)):
            totals[key] += value
        qlake, qout = lake.port.volume_m3/dt, heat.outflow.volume_m3/dt
        if qlake > peak_q:
            peak_q, peak_time = qlake, time
        tout = heat.outflow.sensible_energy_J / (constants.rc * heat.outflow.volume_m3) if heat.outflow.volume_m3 else None
        rows.append(snapshot(dt, qlake, qout, tout, heat.latent_energy_J/dt))
        if time == f["until_s"] and time < case["duration_s"]:
            segment += 1
    water_scale = max(initial.volume_m3 + totals["inflow_m3"] + totals["meltwater_m3"], 1.0)
    energy_scale = max(initial.sensible_energy_J + abs(totals["inflow_J"]) + abs(totals["external_J"]), 1.0)
    summary = dict(model=MODEL, scenario=case["name"], melting=melting, limits=LIMITS,
        constants=asdict(constants), accepted_steps=len(rows)-1, geometry_step_rejections=rejected,
        final_stage_m=rows[-1]["lake_stage_m"], final_temperature_K=state.temperature(constants),
        final_sill_m=sill, sill_lowering_m=initial_sill-sill, totals=totals,
        peak_lake_release_m3_s=peak_q, sampled_peak_time_s=peak_time,
        peak_definition="Maximum initial instantaneous or accepted-step mean lake-to-control discharge; requires refinement",
        max_water_residual_m3=max(abs(x["water_residual_m3"]) for x in rows),
        max_energy_residual_J=max(abs(x["energy_residual_J"]) for x in rows),
        water_scale_m3=water_scale, energy_scale_J=energy_scale,
        max_water_residual_relative=max(abs(x["water_residual_m3"]) for x in rows)/water_scale,
        max_energy_residual_relative=max(abs(x["energy_residual_J"]) for x in rows)/energy_scale,
        geometry_latent_volume_residual_m3=area*(initial_sill-sill)-totals["ice_loss_m3"])
    return rows, summary


def write_run(path, rows, summary):
    path = Path(path)
    path.mkdir(parents=True, exist_ok=False)
    with (path/"timeseries.csv").open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (path/"summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False)+"\n")


def case_digest(case):
    return hashlib.sha256(json.dumps(case, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()
