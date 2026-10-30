"""Strict export comparison and independently recomputed inventory budgets."""

import math
from typing import Dict, List, Mapping, Sequence


class InvalidExport(ValueError):
    """Missing, non-finite, misaligned or physically invalid exported values."""


def validate_times(rows: Sequence[Mapping[str, float]]) -> None:
    if len(rows) < 2:
        raise InvalidExport("At least two output rows are required")
    previous = -math.inf
    for row in rows:
        if "time_s" not in row or not math.isfinite(row["time_s"]):
            raise InvalidExport("Every row needs a finite time_s")
        t = row["time_s"]
        if t < 0 or t <= previous:
            raise InvalidExport("Times must be nonnegative and strictly increasing")
        previous = t


def compare_rows(reference: Sequence[Mapping[str, float]],
                 candidate: Sequence[Mapping[str, float]],
                 rtol: float = 0.01, atol: float = 1e-8) -> Dict:
    """Pointwise tolerance: abs(error) <= atol + rtol * abs(reference).

    atol has each field's SI units. Reported errors remain available independently
    of this configurable numerical gate. No interpolation hides time alignment.
    """
    if not math.isfinite(rtol) or not math.isfinite(atol) or min(rtol, atol) < 0:
        raise ValueError("Tolerances must be finite and nonnegative")
    validate_times(reference)
    validate_times(candidate)
    if len(reference) != len(candidate):
        raise InvalidExport("Export must contain exactly the requested output times")
    fields = set(reference[0]) - {"time_s"}
    for expected, actual in zip(reference, candidate):
        if not math.isclose(expected["time_s"], actual["time_s"], rel_tol=0, abs_tol=1e-9):
            raise InvalidExport("Export time_s differs from the requested reference grid")
        if not fields.issubset(actual):
            raise InvalidExport("Missing fields: " + ", ".join(sorted(fields - set(actual))))
        if any(not math.isfinite(actual[key]) for key in fields):
            raise InvalidExport("Non-finite candidate value")
        for key in fields:
            # These quantities cannot be negative in the supplied test cases.
            if actual[key] < 0:
                raise InvalidExport("Negative value in nonnegative reference field: " + key)
    result = {}
    for key in sorted(fields):
        errors = [a[key] - e[key] for e, a in zip(reference, candidate)]
        scale = max(abs(e[key]) for e in reference)
        rmse = math.sqrt(math.fsum(error * error for error in errors) / len(errors))
        passed = all(abs(error) <= atol + rtol * abs(e[key])
                     for error, e in zip(errors, reference))
        result[key] = {"rmse": rmse, "bias": math.fsum(errors) / len(errors),
                       "max_abs_error": max(map(abs, errors)),
                       "normalized_rmse": rmse / scale if scale else None,
                       "reference_scale": scale, "passed": passed}
    return {"passed": all(metric["passed"] for metric in result.values()),
            "rows": len(reference), "rtol": rtol, "atol": atol, "fields": result}


MASS_FIELDS = ("liquid_mass_kg", "mass_in_kg", "mass_out_kg", "melt_mass_kg")
ENERGY_FIELDS = ("sensible_energy_j", "latent_energy_j", "boundary_heat_j",
                 "advected_energy_in_j", "advected_energy_out_j", "dissipated_to_heat_j")


def _budget(rows, fields, signs, atol, rtol):
    if not math.isfinite(atol) or not math.isfinite(rtol) or min(atol, rtol) < 0:
        raise ValueError("Budget tolerances must be finite and nonnegative")
    if not fields or not rows:
        raise InvalidExport("Empty budget")
    validate_times(rows)
    for row in rows:
        if any(key not in row or not math.isfinite(row[key]) for key in fields):
            raise InvalidExport("Budget requires all finite fields: " + ", ".join(fields))
    first = rows[0]
    residuals, scales, passed = [], [], True
    for row in rows:
        deltas = [row[key] - first[key] for key in fields]
        residual = math.fsum(sign * delta for sign, delta in zip(signs, deltas))
        scale = max(map(abs, deltas))
        residuals.append(residual)
        scales.append(scale)
        passed = passed and abs(residual) <= atol + rtol * scale
    peak = max(scales)
    return {"passed": passed, "max_abs_residual": max(map(abs, residuals)),
            "max_relative_residual": max((abs(r) / s for r, s in zip(residuals, scales)
                                           if s > 0), default=None),
            "peak_inventory_change": peak,
            "nonzero_inventory_change": peak > 0,
            "atol": atol, "rtol": rtol}


def audit_mass(rows: List[Mapping[str, float]], atol=1e-8, rtol=1e-8) -> Dict:
    """Delta liquid = boundary in - boundary out + temperate ice melted.

    No refreezing in this contract. All boundary/melt entries are cumulative kg.
    """
    result = _budget(rows, MASS_FIELDS, (1, -1, 1, -1), atol, rtol)
    if any(row["liquid_mass_kg"] < 0 for row in rows):
        raise InvalidExport("Liquid storage cannot be negative")
    for key in MASS_FIELDS[1:]:
        if any(rows[i][key] < rows[i - 1][key] for i in range(1, len(rows))):
            raise InvalidExport("Cumulative mass fields must not decrease: " + key)
    return result


def audit_energy(rows: List[Mapping[str, float]], atol=1e-6, rtol=1e-8) -> Dict:
    """Delta sensible + latent = net boundary heat + advected in - out + dissipation.

    This is a thermal-energy ledger, not a full mechanical-energy conservation
    test. Each source must be counted once; the checker cannot infer its origin.
    """
    return _budget(rows, ENERGY_FIELDS, (1, 1, -1, -1, 1, -1), atol, rtol)


def observed_order(errors: Sequence[float], refinement_ratio=2.0) -> List[float]:
    """Orders for successive errors against an independent known solution.

    Caller must establish asymptotic behavior; this is not a Richardson estimator
    from three numerical solutions. Zero errors have no defined logarithmic order.
    """
    if not math.isfinite(refinement_ratio) or refinement_ratio <= 1:
        raise ValueError("Refinement ratio must exceed one")
    if len(errors) < 2 or any(not math.isfinite(e) or e <= 0 for e in errors):
        raise ValueError("At least two finite positive reference errors are required")
    return [math.log(coarse / fine) / math.log(refinement_ratio)
            for coarse, fine in zip(errors[:-1], errors[1:])]
