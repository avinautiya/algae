"""Small synthetic forcing-order/flow-commutator diagnostic; no field claim."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from thermal_clock_check import G, WIDTH, derivative

INITIAL = np.array([-0.03, 0.0, 0.03])
QA, QB, QORDINARY, TEMPERATURE = 0.0002, 0.003, 0.001, 0.5


def field(state, q):
    return np.asarray(derivative(0.0, state, TEMPERATURE, discharge_total=q))


def jacobian(state, q, epsilon):
    columns = []
    for j in range(len(state)):
        shift = np.zeros_like(state)
        shift[j] = epsilon
        columns.append((field(state + shift, q) - field(state - shift, q)) / (2 * epsilon))
    return np.column_stack(columns)


def partition(state, q):
    beds = np.array([state[0] + state[2], state[1]])
    c = WIDTH * np.sqrt(G) * (2.0 / 3.0) ** 1.5
    fluxes = lambda head: c * np.maximum(head - beds, 0) ** 1.5
    head = brentq(lambda head: sum(fluxes(head)) - q, min(beds), max(beds) + 1, xtol=1e-14)
    flow = fluxes(head)
    return flow / q


def advance(state, q, duration):
    result = solve_ivp(lambda t, x: field(x, q), (0, duration), state,
                       method="DOP853", rtol=1e-12, atol=1e-15, max_step=duration / 16)
    if not result.success:
        raise RuntimeError(result.message)
    fractions = [partition(x, q) for x in result.y.T]
    return result.y[:, -1], float(np.min(fractions))


def experiment(initial=INITIAL, qa=QA, qb=QB, ordinary=QORDINARY):
    brackets = {str(eps): (jacobian(initial, qb, eps) @ field(initial, qa) -
                           jacobian(initial, qa, eps) @ field(initial, qb)).tolist()
                for eps in [1e-5, 1e-6, 1e-7]}
    bracket = np.array(brackets["1e-07"])
    cases = []
    for h in [600.0, 300.0, 150.0, 75.0]:
        a, amin = advance(initial, qa, h)
        ab, abmin = advance(a, qb, h)
        b, bmin = advance(initial, qb, h)
        ba, bamin = advance(b, qa, h)
        difference = ab - ba
        predicted = h * h * bracket
        ordinary_ab, min_ab = advance(ab, ordinary, 600.0)
        ordinary_ba, min_ba = advance(ba, ordinary, 600.0)
        cases.append({"block_seconds": h,
                      "AB_minus_BA_m": difference.tolist(),
                      "bracket_prediction_m": predicted.tolist(),
                      "relative_bracket_error": float(np.linalg.norm(difference-predicted) /
                                                       np.linalg.norm(predicted)),
                      "difference_fraction_of_total_incision": float(np.linalg.norm(difference[:2]) /
                                                                     np.linalg.norm((ab-initial)[:2])),
                      "difference_divided_by_h_squared": (difference / h**2).tolist(),
                      "minimum_route_fraction_during_protocol": min(amin, abmin, bmin, bamin),
                      "ordinary_flow_initial_toe_fraction_difference": float(
                          partition(ab, ordinary)[0] - partition(ba, ordinary)[0]),
                      "after_common_600s_geometry_difference_m": (ordinary_ab-ordinary_ba).tolist(),
                      "after_common_600s_toe_fraction_difference": float(
                          partition(ordinary_ab, ordinary)[0]-partition(ordinary_ba, ordinary)[0]),
                      "minimum_route_fraction_during_common_return": min(min_ab, min_ba)})
    errors = [item["relative_bracket_error"] for item in cases]
    differences = [np.linalg.norm(item["AB_minus_BA_m"]) for item in cases]
    return {"evidence": "Synthetic smooth two-route implementation mechanism; not environmental novelty",
            "protocol": "AB vs BA, each state same duration, common inlet superheat 0.5 K",
            "QA_m3_s": qa, "QB_m3_s": qb, "ordinary_Q_m3_s": ordinary,
            "initial_state_m": initial.tolist(),
            "budget_note": "AB and BA have identical water, inlet heat and joint Q,T distribution",
            "ordinary_return_note": "Quasi-steady repartition immediately at common Q, then 600 s common evolution",
            "finite_difference_brackets": brackets, "cases": cases,
            "observed_halving_orders": [float(np.log2(a/b)) for a,b in zip(differences[:-1],differences[1:])],
            "assertions": {"all_routes_remain_wet": all(c["minimum_route_fraction_during_protocol"]>0.01 for c in cases),
                           "bracket_error_decreases": all(a>b for a,b in zip(errors[:-1],errors[1:])),
                           "smallest_block_bracket_error_below_5percent": errors[-1]<0.05}}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    record = {
        "primary_environmental_scale_baseline": experiment(np.array([-0.04, 0.0, 0.03]),
                                                          0.035, 0.105, 0.07),
        "additional_low_Q_sensitivity_not_environmental_validation": experiment(),
        "interpretation": "Baseline first, fixed h=200 W/m2/K. Lower-Q case probes finite heat-depletion sensitivity, not an improved environmental prediction.",
    }
    text = json.dumps(record, indent=2) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text)
    if not all(value for group in record.values() if isinstance(group, dict)
               for value in group["assertions"].values()):
        raise SystemExit("Diagnostic assumption or convergence assertion failed")
