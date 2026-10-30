"""Closed-form limiting cases derived independently of the production solver.

All times are seconds and all quantities use SI units. These cases verify selected
equations under restrictive assumptions; they do not validate glacier predictions.
"""

import math
from dataclasses import dataclass
from typing import Callable, Dict, List


@dataclass(frozen=True)
class Case:
    name: str
    description: str
    parameters: Dict[str, float]
    units: Dict[str, str]
    duration_s: float
    samples: int
    solution: Callable[[float], Dict[str, float]]

    def rows(self) -> List[Dict[str, float]]:
        return [dict(time_s=t, **self.solution(t)) for t in
                (self.duration_s * i / (self.samples - 1) for i in range(self.samples))]


def analytical_cases() -> List[Case]:
    """Fixed, reproducible cases; parameter values are idealized, not field data."""
    cases = []

    # A dh/dt = -G h. This is a linear hydraulic element, not Manning flow.
    area, h0, conductance = 10.0, 1.0, 0.02
    def linear(t):
        h = h0 * math.exp(-conductance * t / area)
        return {"level_m": h, "storage_m3": area * h,
                "discharge_m3_s": conductance * h,
                "outflow_m3": area * (h0 - h)}
    cases.append(Case("linear_reservoir", "Constant-area reservoir with Q=G h",
                      {"area_m2": area, "initial_head_m": h0,
                       "conductance_m2_s": conductance},
                      {"level_m": "m", "storage_m3": "m3", "discharge_m3_s": "m3/s",
                       "outflow_m3": "m3"}, 2000.0, 101, linear))

    # A dh/dt = -C sqrt(h); after emptying, h and Q remain zero.
    area_t, head_t, coeff_t = 4.0, 1.0, 0.04
    dry_time = 2 * area_t * math.sqrt(head_t) / coeff_t
    def torricelli(t):
        root = max(0.0, math.sqrt(head_t) - coeff_t * t / (2 * area_t))
        h = root * root
        return {"level_m": h, "storage_m3": area_t * h,
                "discharge_m3_s": coeff_t * root,
                "outflow_m3": area_t * (head_t - h)}
    cases.append(Case("finite_drainage", "Ideal reservoir outlet Q=C sqrt(h), including emptying",
                      {"area_m2": area_t, "initial_head_m": head_t,
                       "coefficient_m5_2_s": coeff_t, "emptying_time_s": dry_time},
                      {"level_m": "m", "storage_m3": "m3", "discharge_m3_s": "m3/s",
                       "outflow_m3": "m3"}, 300.0, 121, torricelli))

    # Equal outlet elevations; independently known conductances. Symmetric and
    # asymmetric cases prevent a hard-coded 50/50 split from passing both.
    def split_case(name, g1, g2):
        a, initial = 8.0, 0.8
        total = g1 + g2
        def solution(t):
            h = initial * math.exp(-total * t / a)
            drained = a * (initial - h)
            return {"level_m": h, "storage_m3": a * h,
                    "outlet_1_m3_s": g1 * h, "outlet_2_m3_s": g2 * h,
                    "outlet_1_volume_m3": drained * g1 / total,
                    "outlet_2_volume_m3": drained * g2 / total}
        return Case(name, "Parallel linear outlets at a common datum; no channel storage",
                    {"area_m2": a, "initial_head_m": initial,
                     "conductance_1_m2_s": g1, "conductance_2_m2_s": g2},
                    {"level_m": "m", "storage_m3": "m3", "outlet_1_m3_s": "m3/s",
                     "outlet_2_m3_s": "m3/s", "outlet_1_volume_m3": "m3",
                     "outlet_2_volume_m3": "m3"}, 1600.0, 101, solution)
    cases.extend([split_case("symmetric_outlets", 0.01, 0.01),
                  split_case("unequal_outlets", 0.01, 0.03)])

    # Up to first overtopping ONLY. No discharge or post-breach claim.
    a_s, initial_s, sill, inflow = 5.0, 0.1, 0.5, 0.002
    overtopping = a_s * (sill - initial_s) / inflow
    def filling(t):
        return {"level_m": initial_s + inflow * t / a_s,
                "storage_m3": a_s * initial_s + inflow * t,
                "inflow_volume_m3": inflow * t, "discharge_m3_s": 0.0}
    cases.append(Case("pre_overtopping", "Filling a constant-area basin up to first contact with a sill",
                      {"area_m2": a_s, "initial_head_m": initial_s,
                       "sill_height_m": sill, "inflow_m3_s": inflow,
                       "first_sill_contact_s": overtopping},
                      {"level_m": "m", "storage_m3": "m3", "inflow_volume_m3": "m3",
                       "discharge_m3_s": "m3/s"}, overtopping, 101, filling))

    # T_t + u T_x = -beta T, prescribed uniform velocity, T(0,t)=Tin,
    # T(x,0)=0 for x>0. No diffusion. Evaluate inlet state at the front.
    speed, x, beta, inlet = 0.2, 5.0, 0.01, 2.0
    arrival = x / speed
    def thermal(t):
        temp = inlet * math.exp(-beta * x / speed) if t >= arrival else 0.0
        return {"temperature_c": temp}
    cases.append(Case("thermal_front", "Plug-flow temperature step with linear wall cooling",
                      {"velocity_m_s": speed, "probe_x_m": x,
                       "cooling_rate_s_1": beta, "inlet_temperature_c": inlet,
                       "arrival_time_s": arrival}, {"temperature_c": "degC"},
                      60.0, 121, thermal))

    # Fixed-geometry heat exchanger, temperate ice, constant effective wall
    # conductance. This verifies the integrated heat-to-melt conversion, not
    # geometry evolution or the validity of a convection closure.
    flow, rho_w, cp, tin = 0.001, 1000.0, 4180.0, 2.0
    ua, latent, rho_i = 200.0, 334000.0, 917.0
    tout = tin * math.exp(-ua / (rho_w * flow * cp))
    power = rho_w * flow * cp * (tin - tout)
    def melting(t):
        energy = power * t
        mass = energy / latent
        return {"outlet_temperature_c": tout, "heat_to_ice_j": energy,
                "melt_mass_kg": mass, "ice_volume_lost_m3": mass / rho_i,
                "meltwater_volume_m3": mass / rho_w}
    cases.append(Case("heat_to_melt", "Fixed-geometry warm-water heat exchanger against temperate ice",
                      {"discharge_m3_s": flow, "water_density_kg_m3": rho_w,
                       "water_cp_j_kg_k": cp, "inlet_temperature_c": tin,
                       "wall_conductance_w_k": ua, "latent_heat_j_kg": latent,
                       "ice_density_kg_m3": rho_i},
                      {"outlet_temperature_c": "degC", "heat_to_ice_j": "J",
                       "melt_mass_kg": "kg", "ice_volume_lost_m3": "m3",
                       "meltwater_volume_m3": "m3"}, 3600.0, 61, melting))
    return cases
