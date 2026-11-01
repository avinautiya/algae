"""Unfitted thermal-component comparators; these are not new models or field scores.

All inputs SI. Fully developed tube correlations are only conditional
comparators for a rough, free-surface ice channel. No target-based tuning here.
"""
from dataclasses import dataclass
import math


def _finite(value, name, positive=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f'{name} must be finite numeric')
    if positive and value <= 0:
        raise ValueError(f'{name} must be positive')
    return float(value)


@dataclass(frozen=True)
class WaterIceProperties:
    nu_m2_s: float
    prandtl: float
    water_k_W_m_K: float
    ice_density_kg_m3: float
    latent_J_kg: float

    def __post_init__(self):
        for name in self.__dataclass_fields__:
            _finite(getattr(self, name), name, positive=True)


def nusselt(re, pr, method, darcy_factor=None):
    re = _finite(re, 'Re', True)
    pr = _finite(pr, 'Pr', True)
    if method in ('db_published_04', 'db_cooling_03'):
        if re < 10000 or not .6 <= pr <= 160:
            raise ValueError('Dittus-Boelter comparator outside numerical Re/Pr gate')
        if darcy_factor is not None:
            raise ValueError('Dittus-Boelter does not take a friction factor')
        exponent = .4 if method == 'db_published_04' else .3
        return .023 * re**.8 * pr**exponent
    if method not in ('gnielinski_smooth', 'gnielinski_supplied_f'):
        raise ValueError('unknown comparator')
    # Conservative lower gate 3000, stricter than the MOOSE implementation's 2300.
    if not 3000 <= re <= 5e6 or not .5 <= pr <= 2000:
        raise ValueError('Gnielinski comparator outside numerical Re/Pr gate')
    if method == 'gnielinski_smooth':
        if darcy_factor is not None:
            raise ValueError('smooth comparator calculates its own friction factor')
        f = (1.82 * math.log10(re) - 1.64)**-2
    else:
        f = _finite(darcy_factor, 'supplied Darcy friction factor', True)
    return (f/8) * (re-1000) * pr / (1 + 12.7 * math.sqrt(f/8) * (pr**(2/3)-1))


def heat_transfer(velocity_m_s, hydraulic_diameter_m, properties, method, darcy_factor=None):
    velocity_m_s = _finite(velocity_m_s, 'velocity', True)
    diameter = _finite(hydraulic_diameter_m, 'hydraulic diameter', True)
    re = velocity_m_s * diameter / properties.nu_m2_s
    nu = nusselt(re, properties.prandtl, method, darcy_factor)
    return {'Re': re, 'Nu': nu, 'h_W_m2_K': nu * properties.water_k_W_m_K / diameter,
            'applicability': 'numerical_Re_Pr_gate_only; open_channel_transfer_unvalidated'}


def integrate_constant_cells(cells, start_s, end_s, properties, method):
    """Integrate documented, contiguous imposed cell means over exact overlap.

    Each cell: start_s,end_s,temperature_C,velocity_m_s,hydraulic_diameter_m,
    optionally darcy_factor. No extrapolation, gap filling, or temperature offset.
    Input cells must partition the requested interval. A cell's Tw and velocity
    are constant assumptions, even if constructed from measured averages.
    """
    a = _finite(start_s, 'interval start'); b = _finite(end_s, 'interval end')
    if b <= a:
        raise ValueError('positive duration required')
    cursor = a; energy = 0.; used = []
    previous_end = None
    for cell in cells:
        left = _finite(cell['start_s'], 'cell start'); right = _finite(cell['end_s'], 'cell end')
        if right <= left or (previous_end is not None and left < previous_end):
            raise ValueError('unordered, overlapping or nonpositive cells')
        previous_end = right
        lo, hi = max(left, a), min(right, b)
        if hi <= lo:
            continue
        if not math.isclose(lo, cursor, abs_tol=1e-6, rel_tol=0):
            raise ValueError('forcing gap or unsupported interval endpoint')
        temp = _finite(cell['temperature_C'], 'temperature')
        if temp < 0:
            raise ValueError('negative water temperature/refreezing unsupported')
        h = heat_transfer(cell['velocity_m_s'], cell['hydraulic_diameter_m'], properties,
                          method, cell.get('darcy_factor'))
        supply = h['h_W_m2_K'] * temp * (hi-lo)
        energy += supply
        used.append(dict(start_s=lo, end_s=hi, energy_J_m2=supply, **h))
        cursor = hi
    if not math.isclose(cursor, b, abs_tol=1e-6, rel_tol=0):
        raise ValueError('forcing does not cover interval')
    return {'thermal_energy_J_m2': energy,
            'predicted_lowering_m': energy / (properties.ice_density_kg_m3 * properties.latent_J_kg),
            'cells': used, 'scope': 'conditional local heat-transfer hindcast; ice at melting point'}
