"""Conditional flood-volume bounds; no field calibration or CFD implied."""
from dataclasses import dataclass, asdict
import argparse
import json
import math
from scipy.integrate import quad


@dataclass(frozen=True)
class Bound:
    critical_peak_m3_s: float
    minimum_time_s: float
    minimum_volume_m3: float
    reachable_at_peak: bool
    scope: str = "static toe and bypass; first zero-depth contact only"


def static_bound(*, area_m2, coefficient, exponent, bypass_m, toe_m,
                 initial_stage_m, peak_m3_s):
    """Exact fastest/least-total-water first-contact bound for declared ODE.

    A H' = I - C max(H-zb,0)**p, 0 <= I <= peak.
    Inputs use a common vertical datum. Water means TOTAL inflow, not excess
    above a baseline. Existing pond storage is represented by initial_stage.
    """
    vals = (area_m2, coefficient, exponent, bypass_m, toe_m,
            initial_stage_m, peak_m3_s)
    if not all(math.isfinite(v) and not isinstance(v, bool) for v in vals):
        raise ValueError("all inputs must be finite numbers, not bools")
    if min(area_m2, coefficient, exponent, peak_m3_s) <= 0:
        raise ValueError("area, coefficient, exponent and peak must be positive")
    critical = coefficient * max(toe_m-bypass_m, 0.)**exponent
    if initial_stage_m >= toe_m:
        return Bound(critical, 0., 0., True)
    if peak_m3_s <= critical:
        return Bound(critical, math.inf, math.inf, False)
    t, _ = quad(lambda h: area_m2 / (
        peak_m3_s-coefficient*max(h-bypass_m, 0.)**exponent),
        initial_stage_m, toe_m, epsabs=1e-9, epsrel=1e-10)
    return Bound(critical, t, peak_m3_s*t, True)


def deadline_certificate(bound, *, available_volume_m3, available_time_s):
    """Necessary and sufficient static-ODE condition, not sustained wetness."""
    if not all(math.isfinite(x) and x >= 0 for x in (
            available_volume_m3, available_time_s)):
        raise ValueError("budgets must be finite and nonnegative")
    return (bound.reachable_at_peak and
            available_volume_m3 >= bound.minimum_volume_m3 and
            available_time_s >= bound.minimum_time_s)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--area", type=float, required=True)
    p.add_argument("--coefficient", type=float, required=True)
    p.add_argument("--exponent", type=float, default=1.5)
    p.add_argument("--bypass", type=float, required=True)
    p.add_argument("--toe", type=float, required=True)
    p.add_argument("--stage", type=float, required=True)
    p.add_argument("--peak", type=float, required=True)
    a = p.parse_args()
    b = static_bound(area_m2=a.area, coefficient=a.coefficient,
                     exponent=a.exponent, bypass_m=a.bypass, toe_m=a.toe,
                     initial_stage_m=a.stage, peak_m3_s=a.peak)
    # Standard JSON: infinity is represented explicitly as null.
    print(json.dumps({k: (None if isinstance(v, float) and not math.isfinite(v)
                          else v) for k, v in asdict(b).items()}, indent=2))


if __name__ == "__main__":
    main()
