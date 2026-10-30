"""Compare Claude's actual thermal kernel against 60-digit segment solutions.

Usage: python -m ice_network_verification.audit_heat_kernel SNAPSHOT_ROOT --report OUT
numpy is required for the production function's arrays. No Numba installation is
needed: the AST-extracted function bodies run in Python without JIT decorators.
This tests Python arithmetic, not Numba compilation or physical closure validity.
"""

import argparse
import ast
from decimal import Decimal, localcontext
import hashlib
import json
from pathlib import Path

from .thermal import temperate_segment


def load_kernel(root):
    import numpy as np
    source = root / "ice_network/src/ice_network/thermal/channel_heat.py"
    constants = root / "ice_network/src/ice_network/constants.py"
    names = {}
    for statement in ast.parse(constants.read_text()).body:
        if isinstance(statement, ast.Assign) and len(statement.targets) == 1:
            names[statement.targets[0].id] = ast.literal_eval(statement.value)
    functions = []
    for statement in ast.parse(source.read_text()).body:
        if isinstance(statement, ast.FunctionDef) and statement.name in (
                "heat_transfer_coeff", "reach_heat"):
            statement.decorator_list = []
            functions.append(statement)
    if len(functions) != 2:
        raise ValueError("Expected exactly two known thermal functions")
    namespace = dict(names, np=np)
    module = ast.Module(body=functions, type_ignores=[])
    exec(compile(ast.fix_missing_locations(module), str(source), "exec"), namespace)
    return namespace, {str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
                       for path in (source, constants)}


def high_precision_segment(q, tin, tm, width, depth, length, phi, factor, constants):
    # Independent evaluation of the stated ODE, with no use of reach_heat or its
    # float equilibrium-temperature subtraction. Constants use decimal strings.
    with localcontext() as context:
        context.prec = 60
        d = lambda value: Decimal(str(value))
        q, tin, tm, w, h, length, phi, factor = map(
            d, (q, tin, tm, width, depth, length, phi, factor))
        rho, cp, kw, nuw, pr = (d(constants[key]) for key in
                                 ("RHO_W", "CP_W", "K_W", "NU_W", "PR_W"))
        perimeter = w + 2 * h
        dh = 4 * w * h / perimeter
        reynolds = max(q * dh / (w * h * nuw), Decimal(1))
        htc = factor * d(.023) * reynolds ** d(.8) * pr ** d(.3) * kw / dh
        capacity_rate = rho * cp * q
        a = htc * perimeter * length / capacity_rate
        decay = (-a).exp()
        transfer = 1 - decay
        tout = tm + (tin - tm) * decay + phi / capacity_rate * transfer / a
        heat = capacity_rate * (tin - tm) * transfer + phi * (1 - transfer / a)
        return float(tout), float(heat), float(a), float(reynolds)


def run(root):
    import numpy as np
    kernel, hashes = load_kernel(root)
    scenarios = [
        ("warm_no_dissipation", .001, 2., 0., 1.),
        ("cold_with_dissipation", .001, 0., 20., 1.),
        ("warm_with_dissipation", .001, 2., 20., 1.),
        ("high_discharge", .1, 2., 20., 1.),
        ("low_discharge", 1e-5, .1, .01, 1.),
        ("weak_transfer_stress", .001, 2., 20., 1e-12),
        ("weak_transfer_cold_stress", .001, 0., 20., 1e-12),
    ]
    results = []
    for name, q, tin, total_phi, factor in scenarios:
        for segments in (1, 10, 100):
            width, depth, length, tm = .5, .1, 10., 0.
            h = np.full(segments + 1, depth)
            w = np.full(segments + 1, width)
            temperatures = np.empty(segments + 1)
            heat = np.empty(segments)
            phi = np.full(segments, total_phi / segments)
            actual_tout = kernel["reach_heat"](q, tin, tm, h, w, length / segments,
                                                phi, 0., factor, temperatures, heat)
            actual_heat = float(heat.sum())
            expected_tout, expected_heat, a, re = high_precision_segment(
                q, tin, tm, width, depth, length, total_phi, factor, kernel)
            capacity = kernel["RHO_W"] * kernel["CP_W"] * q
            stable_tout, stable_heat = temperate_segment(capacity, a * capacity,
                                                        tin, total_phi, tm)
            energy_residual = (kernel["RHO_W"] * kernel["CP_W"] * q *
                               (tin - actual_tout) + total_phi - actual_heat)
            incoming = kernel["RHO_W"] * kernel["CP_W"] * q * tin + total_phi
            results.append({"scenario": name, "segments": segments,
                            "q_m3_s": q, "inlet_temperature_c": tin,
                            "dissipation_w": total_phi, "nu_factor": factor,
                            "thermal_exposure": a, "reynolds": re,
                            "below_nominal_pipe_turbulence_re_10000": re < 10000,
                            "stress_test_not_physical_calibration": factor != 1,
                            "actual_outlet_temperature_c": float(actual_tout),
                            "reference_outlet_temperature_c": expected_tout,
                            "actual_heat_to_ice_w": actual_heat,
                            "reference_heat_to_ice_w": expected_heat,
                            "stable_heat_to_ice_w": stable_heat,
                            "stable_heat_abs_error_w": abs(stable_heat - expected_heat),
                            "stable_temperature_abs_error_c": abs(stable_tout - expected_tout),
                            "temperature_abs_error_c": abs(actual_tout - expected_tout),
                            "heat_abs_error_w": abs(actual_heat - expected_heat),
                            "heat_relative_error": (abs(actual_heat - expected_heat) /
                                                    abs(expected_heat)) if expected_heat else None,
                            "heat_error_fraction_of_input_power":
                                abs(actual_heat - expected_heat) / incoming if incoming else None,
                            "thermal_budget_residual_w": float(energy_residual),
                            "negative_segment_heat": bool((heat < 0).any())})
    return {"evidence_type": "independent_kernel_numerical_audit",
            "source_sha256": hashes, "execution": "AST-extracted original Python bodies; no JIT",
            "not_validated": ["Dittus-Boelter for open channels", "geometry evolution",
                              "hydraulics", "field melt", "Numba behavior"],
            "cases": results}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("snapshot_root", type=Path)
    parser.add_argument("--report", required=True, type=Path)
    args = parser.parse_args()
    report = run(args.snapshot_root)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    ordinary = [row for row in report["cases"] if not row["stress_test_not_physical_calibration"]]
    stress = [row for row in report["cases"] if row["stress_test_not_physical_calibration"]]
    print(json.dumps({"ordinary_max_heat_error_w": max(r["heat_abs_error_w"] for r in ordinary),
                      "stress_max_heat_error_w": max(r["heat_abs_error_w"] for r in stress),
                      "negative_segment_heat_cases": sum(r["negative_segment_heat"] for r in report["cases"]),
                      "report": str(args.report)}, indent=2))


if __name__ == "__main__":
    main()
