"""Monotone passive parallel-route ratings; all inputs and outputs use SI units."""
from dataclasses import dataclass
import math

G = 9.80665


def _number(value):
    return not isinstance(value,bool) and isinstance(value,(int,float)) and math.isfinite(value)


@dataclass(frozen=True)
class Route:
    name: str
    floor_z_m: float
    debris_height_m: float
    width_m: float
    weir_coefficient_sqrt_m_per_s: float
    aperture_height_m: float = 0.0
    discharge_coefficient: float = 0.0
    darcy_conductance_m2_s: float = 0.0

    def __post_init__(self):
        if not isinstance(self.name,str) or not self.name.strip():
            raise ValueError('nonempty route name required')
        for name, value in vars(self).items():
            if name != 'name' and not _number(value):
                raise ValueError(name + ' must be finite')
        for name in ('debris_height_m', 'width_m', 'weir_coefficient_sqrt_m_per_s',
                     'aperture_height_m', 'discharge_coefficient', 'darcy_conductance_m2_s'):
            if getattr(self, name) < 0:
                raise ValueError(name + ' must be nonnegative')
        if self.aperture_height_m > self.debris_height_m:
            raise ValueError('aperture must fit below debris top')
        if self.aperture_height_m and self.darcy_conductance_m2_s:
            raise ValueError('choose aperture OR Darcy transmission, not both')


def route_discharge_m3_s(route, stage_m, downstream_stage_m):
    """One-way passive rating, zero for head below the actual floor or receiver.

    Aperture area is width * wetted aperture height. Darcy conductance is a
    whole-route effective conductance (m²/s), already including area/length.
    """
    if not all(map(_number, (stage_m, downstream_stage_m))):
        raise ValueError('stages must be finite')
    head = max(0.0, stage_m - max(route.floor_z_m, downstream_stage_m))
    crest = route.floor_z_m + route.debris_height_m
    over_head = max(0.0, stage_m - max(crest, downstream_stage_m))
    overtopping = route.weir_coefficient_sqrt_m_per_s * route.width_m * over_head ** 1.5
    # Integrate the hydrostatic head over a vertical aperture rather than
    # assigning its bottom pressure to every point in the opening.
    H=max(0.,stage_m-route.floor_z_m)
    a=min(route.aperture_height_m,H)
    receiver=max(0.,downstream_stage_m-route.floor_z_m)
    submerged=min(a,receiver)
    integral=submerged*math.sqrt(head)
    free_span=max(0.,a-receiver)
    if free_span>0 and head>0:
        # Integral sqrt(H-z) dz on z in [receiver,a], stable for tiny apertures.
        fraction=min(1.,free_span/head)
        drop=1. if fraction==1. else -math.expm1(1.5*math.log1p(-fraction))
        integral+=(2/3)*head**1.5*drop
    under=(route.discharge_coefficient*route.width_m*math.sqrt(2*G)*integral+
           route.darcy_conductance_m2_s*head)
    return overtopping + under


@dataclass(frozen=True)
class StageSolution:
    stage_m: float
    route_flows_m3_s: dict
    residual_m3_s: float
    iterations: int


def solve_common_stage(routes, inflow_m3_s, downstream_stage_m, *,
                       flow_tolerance_m3_s=1e-9, max_iterations=256,
                       max_stage_m=None):
    """Bisection of summed monotone ratings for >=2 disjoint parallel routes.

    max_stage_m is an optional physical stage bound; an unreachable inflow
    raises ValueError, never an invented transmission or overflow path.
    """
    routes = tuple(routes)
    if len(routes) < 2 or len({r.name for r in routes}) != len(routes):
        raise ValueError('supply at least two routes with unique names')
    if not _number(inflow_m3_s) or inflow_m3_s < 0:
        raise ValueError('inflow must be finite and nonnegative')
    if not _number(downstream_stage_m):
        raise ValueError('downstream stage must be finite')
    if not _number(flow_tolerance_m3_s) or flow_tolerance_m3_s <= 0 or type(max_iterations) is not int or max_iterations < 1:
        raise ValueError('positive tolerance and iteration limit required')
    def result(h, iterations):
        flows = {r.name: route_discharge_m3_s(r, h, downstream_stage_m) for r in routes}
        return StageSolution(h, flows, sum(flows.values())-inflow_m3_s, iterations)
    low = min(downstream_stage_m, min(r.floor_z_m for r in routes))
    if inflow_m3_s == 0:
        return result(low, 0)
    if not any((r.width_m and (r.weir_coefficient_sqrt_m_per_s or
                r.aperture_height_m*r.discharge_coefficient)) or r.darcy_conductance_m2_s for r in routes):
        raise ValueError('positive inflow cannot pass zero-transmission routes')
    if max_stage_m is not None:
        if not _number(max_stage_m) or max_stage_m < low:
            raise ValueError('invalid maximum stage')
        high = max_stage_m
    else:
        high = max(downstream_stage_m, max(r.floor_z_m+r.debris_height_m for r in routes))+1
        for _ in range(max_iterations):
            if result(high, 0).residual_m3_s >= 0:
                break
            high = low + 2*(high-low)
        else:
            raise ValueError('unable to bracket inflow')
    if result(high, 0).residual_m3_s < -flow_tolerance_m3_s:
        raise ValueError('inflow exceeds discharge at maximum stage')
    for iteration in range(1, max_iterations+1):
        mid = (low+high)/2
        answer = result(mid, iteration)
        if abs(answer.residual_m3_s) <= flow_tolerance_m3_s:
            return answer
        if mid == low or mid == high:
            break
        if answer.residual_m3_s < 0:
            low = mid
        else:
            high = mid
    raise RuntimeError('bisection did not achieve requested flow tolerance')


def hydrostatic_diagnostics(depth_m, length_m, speed_m_s, bed_slope,
                            curvature_radius_m=None):
    """Dimensionless indicators only; no indicator validates hydrostatic flow."""
    values = (depth_m, length_m, speed_m_s, bed_slope)
    if not all(map(_number, values)) or depth_m <= 0 or length_m <= 0 or speed_m_s < 0:
        raise ValueError('finite inputs, positive depth/length, nonnegative speed required')
    if curvature_radius_m is not None and (not _number(curvature_radius_m) or curvature_radius_m <= 0):
        raise ValueError('radius must be positive finite or None (unknown)')
    warnings = ['Indicators are warnings for review; they do not validate a hydrostatic approximation.']
    if curvature_radius_m is None:
        warnings.append('Curvature radius unknown: vertical acceleration cannot be assessed.')
    return {
        'depth_to_length': depth_m/length_m,
        'curvature_acceleration_over_g': None if curvature_radius_m is None else speed_m_s**2/(G*curvature_radius_m),
        'froude_number': speed_m_s/math.sqrt(G*depth_m),
        'bed_slope': bed_slope,
        'warnings': warnings,
    }
