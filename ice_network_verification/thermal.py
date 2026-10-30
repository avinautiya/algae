"""Numerically stable constant-coefficient heat-exchanger segment solution.

This is an independently derived utility, not a change to Claude's production
solver and not validation of any empirical heat-transfer coefficient.
"""

import math
from typing import Tuple


def temperate_segment(capacity_rate_w_k: float, wall_conductance_w_k: float,
                      inlet_temperature_c: float, dissipation_w: float = 0.0,
                      melting_temperature_c: float = 0.0) -> Tuple[float, float]:
    """Return outlet temperature (degC) and heat delivered to ice (W).

    capacity_rate = rho_w * cp * Q; wall_conductance = h_T * P * dx.
    Dissipation is uniform along this segment. Temperate ice; no refreezing.
    Exact ODE solution with the cancellation-prone heat term evaluated by series.
    """
    values = (capacity_rate_w_k, wall_conductance_w_k, inlet_temperature_c,
              dissipation_w, melting_temperature_c)
    if any(not math.isfinite(value) for value in values):
        raise ValueError("Segment inputs must be finite")
    if capacity_rate_w_k <= 0 or wall_conductance_w_k < 0 or dissipation_w < 0:
        raise ValueError("Positive capacity rate and nonnegative conductance/source required")
    if inlet_temperature_c < melting_temperature_c:
        raise ValueError("Subfreezing inlet requires a model supporting refreezing")
    c, ua, phi = capacity_rate_w_k, wall_conductance_w_k, dissipation_w
    a = ua / c
    if not math.isfinite(a):
        raise ValueError("Thermal exposure overflows float representation")
    transfer = -math.expm1(-a)
    if a == 0:
        source_out_fraction, source_ice_fraction = 1.0, 0.0
    elif a < 1e-4:
        # 1 - (1-exp(-a))/a = a/2 - a^2/6 + a^3/24 - a^4/120 + ...
        source_ice_fraction = a * (.5 + a * (-1 / 6 + a * (1 / 24 - a / 120)))
        source_out_fraction = 1 - source_ice_fraction
    else:
        source_out_fraction = transfer / a
        source_ice_fraction = 1 - source_out_fraction
    excess = inlet_temperature_c - melting_temperature_c
    outlet = (melting_temperature_c + excess * math.exp(-a)
              + phi / c * source_out_fraction)
    heat_to_ice = c * excess * transfer + phi * source_ice_fraction
    return outlet, heat_to_ice
