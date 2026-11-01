"""Synthetic two-route dose-collapse falsifier; not glacier validation.

Run: python3 thermal_clock_check.py --output thermal_clock_results.json
The auxiliary stock-export law is illustrative, not calibrated or inferred.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

G = 9.81
Q = 0.07
WIDTH = 0.5
AREA = 1.0
RHO_CP = 1000.0 * 4217.0
RHO_L = 917.0 * 334000.0
HORIZON = 48 * 3600.0


def derivative(_, state, temperature, export_rate=0.0, discharge_total=Q):
    """Flat ideal parallel controls, plug-flow cooling, optional stock clock."""
    z1, z2, stock = state
    beds = (z1 + stock, z2)
    coefficient = WIDTH * math.sqrt(G) * (2.0 / 3.0) ** 1.5

    def discharge(stage, bed):
        return coefficient * max(stage - bed, 0.0) ** 1.5

    stage = brentq(lambda level: sum(discharge(level, b) for b in beds) - discharge_total,
                   min(beds), max(beds) + 1.0, xtol=1e-14)
    heat = []
    for index, bed in enumerate(beds):
        flow = discharge(stage, bed)
        conductance = 200.0 / (1.0 + stock / 0.02) if index == 0 else 200.0
        power = (RHO_CP * flow * temperature *
                 (-math.expm1(-conductance * AREA / (RHO_CP * flow)))) if flow > 0 else 0.0
        heat.append(power)
    return [-heat[0] / (RHO_L * AREA), -heat[1] / (RHO_L * AREA),
            -export_rate * stock]


def run(schedule, export_rate=0.0):
    state = np.array([-0.04, 0.0, 0.03])
    time = 0.0
    evaluations = 0
    for duration, temperature in schedule:
        result = solve_ivp(lambda t, x: derivative(t, x, temperature, export_rate),
                           (time, time + duration), state, method="DOP853",
                           rtol=2e-11, atol=2e-13, max_step=1800)
        if not result.success:
            raise RuntimeError(result.message)
        state = result.y[:, -1]
        evaluations += result.nfev
        time += duration
    if abs(time - HORIZON) > 1e-8:
        raise ValueError("Schedule must span the declared horizon")
    return {"final_state_m": state.tolist(), "rhs_evaluations": evaluations,
            "dose_K_s": sum(d * t for d, t in schedule)}


def experiment():
    schedules = {
        "constant": [(HORIZON, 0.5)],
        "early_heat": [(HORIZON / 2, 1.0), (HORIZON / 2, 0.0)],
        "late_heat": [(HORIZON / 2, 0.0), (HORIZON / 2, 1.0)],
        "four_pulses": [(HORIZON / 8, t) for t in [1., 0., 1., 0., 1., 0., 1., 0.]],
    }
    fixed = {name: run(schedule) for name, schedule in schedules.items()}
    exported = {name: run(schedule, 1.0 / 86400.0) for name, schedule in schedules.items()}
    fixed_error = max(np.max(np.abs(np.array(item["final_state_m"]) -
                                          fixed["constant"]["final_state_m"]))
                      for item in fixed.values())
    early = np.array(exported["early_heat"]["final_state_m"])
    late = np.array(exported["late_heat"]["final_state_m"])
    return {
        "evidence": "Synthetic reduced model only; known time-reparameterization test",
        "assumptions": "Fixed Q, ideal parallel hydraulic controls, homogeneous steady passive heat",
        "counterexample": "Prescribed stock decay in real time; illustrative export rate 1/day",
        "fixed_stock": fixed, "real_time_stock_export": exported,
        "fixed_stock_max_schedule_difference_m": float(fixed_error),
        "export_early_minus_late_geometry_m": (early[:2] - late[:2]).tolist(),
        "assertions": {
            "dose_collapse": bool(fixed_error < 1e-9),
            "extra_clock_breaks_collapse": bool(np.max(np.abs(early[:2] - late[:2])) > 1e-4),
            "export_final_stock_schedule_invariant": bool(abs(early[2] - late[2]) < 1e-10),
        },
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    record = experiment()
    if not all(record["assertions"].values()):
        raise RuntimeError(record)
    content = json.dumps(record, indent=2) + "\n"
    if args.output:
        args.output.write_text(content)
    print(content)
