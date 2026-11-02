"""Local, explicit heat-transfer choices for imposed rectangular open-channel flow.

Tube correlations are conditional engineering closures here, NOT validated
rough/free-surface ice-channel laws. No target-derived gain is supplied.
"""
from __future__ import annotations
from dataclasses import asdict, dataclass
import math

SCHEMA_VERSION = 1
METHODS = ('fixed', 'db04', 'db03', 'gnielinski_smooth')
G = 9.81


class ApplicabilityRefusal(ValueError):
    """A requested closure cannot be evaluated under its declared contract."""


def numeric(x, name, minimum=None, strictly_positive=False):
    if isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x):
        raise ApplicabilityRefusal(name + ' must be a finite numeric value')
    if strictly_positive and x <= 0 or minimum is not None and x < minimum:
        raise ApplicabilityRefusal(name + ' outside its permitted sign range')
    return float(x)


@dataclass(frozen=True)
class ClosureConfig:
    schema_version: int = SCHEMA_VERSION
    method: str = 'fixed'
    fixed_h_W_m2_K: float = 200.
    nu_m2_s: float = 1.8e-6
    water_k_W_m_K: float = .56
    prandtl: float = 13.5
    gain: float = 1.
    calibration_provenance: str | None = None
    domain_policy: str = 'strict'
    extrapolation_reason: str | None = None

    def __post_init__(self):
        if isinstance(self.schema_version, bool) or self.schema_version != SCHEMA_VERSION:
            raise ApplicabilityRefusal('unknown closure schema version')
        if self.method not in METHODS:
            raise ApplicabilityRefusal('unknown heat-transfer method; no automatic fallback')
        for name in ['nu_m2_s', 'water_k_W_m_K', 'prandtl', 'gain']:
            numeric(getattr(self, name), name, strictly_positive=True)
        numeric(self.fixed_h_W_m2_K, 'fixed h', minimum=0.)
        if self.gain != 1. and (not isinstance(self.calibration_provenance, str) or not self.calibration_provenance.strip()):
            raise ApplicabilityRefusal('nonunit gain requires explicit calibration provenance')
        if self.calibration_provenance is not None and (not isinstance(self.calibration_provenance, str) or not self.calibration_provenance.strip()):
            raise ApplicabilityRefusal('calibration provenance must be nonempty text')
        if self.domain_policy not in ('strict', 'allow_high_re_extrapolation'):
            raise ApplicabilityRefusal('unknown domain policy')
        if self.domain_policy != 'strict' and (not isinstance(self.extrapolation_reason, str) or not self.extrapolation_reason.strip()):
            raise ApplicabilityRefusal('high-Re extrapolation requires an explicit reason')
        if self.domain_policy == 'strict' and self.extrapolation_reason is not None:
            raise ApplicabilityRefusal('extrapolation reason supplied for strict policy')


def config_from_dict(value=None, *, legacy_fixed_h_W_m2_K=200.):
    """Absent configuration preserves the legacy supplied fixed coefficient."""
    if value is None:
        return ClosureConfig(fixed_h_W_m2_K=legacy_fixed_h_W_m2_K)
    if not isinstance(value, dict):
        raise ApplicabilityRefusal('closure configuration must be an object')
    allowed = set(ClosureConfig.__dataclass_fields__)
    if set(value) - allowed:
        raise ApplicabilityRefusal('unknown closure fields: ' + ','.join(sorted(set(value)-allowed)))
    value = dict(value)
    if value.get('method', 'fixed') == 'fixed' and 'fixed_h_W_m2_K' not in value:
        value['fixed_h_W_m2_K'] = legacy_fixed_h_W_m2_K
    try:
        return ClosureConfig(**value)
    except TypeError as exc:
        raise ApplicabilityRefusal('invalid closure configuration') from exc


def rectangular_flow(Q_m3_s, width_m, depth_m, nu_m2_s=1.8e-6):
    """Compute local geometry and dimensionless diagnostics; no hydraulic solve."""
    q = numeric(Q_m3_s, 'Q', minimum=0.)
    w = numeric(width_m, 'width', strictly_positive=True)
    d = numeric(depth_m, 'depth', minimum=0.)
    nu = numeric(nu_m2_s, 'kinematic viscosity', strictly_positive=True)
    if q > 0 and d == 0:
        raise ApplicabilityRefusal('positive prescribed flow needs positive water depth')
    if d == 0:
        return dict(Q_m3_s=q, width_m=w, depth_m=d, area_m2=0., wetted_perimeter_m=0.,
                    hydraulic_diameter_m=0., velocity_m_s=0., Re=0., Fr=0.,
                    flowing=False)
    area = numeric(w*d, 'derived area', strictly_positive=True)
    perimeter = numeric(w+2*d, 'derived perimeter', strictly_positive=True)
    diameter = numeric(4*(area/perimeter), 'derived hydraulic diameter', strictly_positive=True)
    u = numeric(q/area, 'derived velocity', minimum=0., strictly_positive=q > 0)
    # Algebraically identical to u*Dh/nu without a cancellation-prone large u*Dh.
    re = numeric(4*(q/perimeter)/nu, 'derived Reynolds number', minimum=0., strictly_positive=q > 0)
    fr = numeric(u/math.sqrt(G*d), 'derived Froude number', minimum=0.)
    return dict(Q_m3_s=q, width_m=w, depth_m=d, area_m2=area,
                wetted_perimeter_m=perimeter, hydraulic_diameter_m=diameter,
                velocity_m_s=u, Re=re, Fr=fr, flowing=q > 0)


def evaluate(Q_m3_s, width_m, depth_m, config=None, *, legacy_fixed_h_W_m2_K=200.):
    if not isinstance(config, ClosureConfig):
        config = config_from_dict(config, legacy_fixed_h_W_m2_K=legacy_fixed_h_W_m2_K)
    flow = rectangular_flow(Q_m3_s, width_m, depth_m, config.nu_m2_s)
    status = ('prescribed_fixed_coefficient_no_fit_claim' if config.method == 'fixed' and config.gain == 1.
              else 'unfitted_correlation' if config.gain == 1. else 'explicit_user_calibration')
    result = dict(flow, closure=asdict(config), method=config.method, warnings=[],
                  calibration_status=status,
                  applicability='imposed-flow engineering closure; rough free-surface transfer not field-validated')
    if not flow['flowing']:
        return dict(result, Nu=None, h_W_m2_K=0., numerical_domain='zero_discharge_no_advective_exchange',
                    standing_water_convection='not represented')
    if config.method == 'fixed':
        h = numeric(config.fixed_h_W_m2_K*config.gain, 'effective fixed h', minimum=0.)
        prescribed_nu = numeric(h*(flow['hydraulic_diameter_m']/config.water_k_W_m_K), 'prescribed Nusselt', minimum=0.)
        return dict(result, Nu=prescribed_nu,
                    h_W_m2_K=h, numerical_domain='prescribed_coefficient_no_correlation_domain_claim')
    re, pr = flow['Re'], config.prandtl
    lower, upper, pl, ph = (1e4, 1.2e5, .7, 160.) if config.method.startswith('db') else (3e3, 5e6, .5, 2000.)
    if re < lower or not pl <= pr <= ph:
        raise ApplicabilityRefusal(f'{config.method}: Re={re:g}, Pr={pr:g} outside lower-Re/Pr domain; no laminar fallback')
    extrapolated = re > upper
    if extrapolated and config.domain_policy == 'strict':
        raise ApplicabilityRefusal(f'{config.method}: Re={re:g} exceeds declared upper-Re gate {upper:g}')
    if extrapolated:
        result['warnings'].append('explicit high-Re extrapolation; not within cited numerical domain')
    if config.method.startswith('db'):
        exponent = .4 if config.method == 'db04' else .3
        nusselt = .023*re**.8*pr**exponent
    else:
        f = (1.82*math.log10(re)-1.64)**-2
        denominator = 1+12.7*math.sqrt(f/8)*(pr**(2/3)-1)
        if denominator <= 0:
            raise ApplicabilityRefusal('nonpositive Gnielinski denominator')
        nusselt = (f/8)*(re-1000)*pr/denominator
        result['smooth_tube_darcy_factor'] = f
    nusselt = numeric(nusselt*config.gain, 'effective Nusselt', strictly_positive=True)
    h = numeric(nusselt*(config.water_k_W_m_K/flow['hydraulic_diameter_m']), 'derived h', strictly_positive=True)
    return dict(result, Nu=nusselt, h_W_m2_K=h,
                numerical_domain='explicit_high_Re_extrapolation' if extrapolated else 'inside_cited_Re_Pr_gate_only')


def finite_heat_chain(Q_m3_s, inlet_superheat_K, conductances_W_K, *, rho_water_kg_m3=1000., cp_water_J_kg_K=4217.):
    """Exact serial zero-C surfaces, with no heat supplied twice.

    No storage or natural convection at Q=0. No signed/refreezing exchange here.
    Use the separate wet-debris kernel for a nonzero-temperature solid surface.
    """
    q = numeric(Q_m3_s, 'Q', minimum=0.)
    tin = numeric(inlet_superheat_K, 'inlet superheat', minimum=0.)
    rho = numeric(rho_water_kg_m3, 'water density', strictly_positive=True)
    cp = numeric(cp_water_J_kg_K, 'water heat capacity', strictly_positive=True)
    if not isinstance(conductances_W_K, (list, tuple)):
        raise ApplicabilityRefusal('ordered conductance list required')
    gs = [numeric(x, 'conductance', minimum=0.) for x in conductances_W_K]
    capacity = numeric(q*rho*cp, 'water heat-capacity rate', minimum=0.)
    if q > 0 and capacity == 0:
        raise ApplicabilityRefusal('positive-flow heat-capacity rate underflow')
    initial = numeric(capacity*tin, 'incoming sensible heat', minimum=0.)
    remaining = initial; rows=[]
    for conductance in gs:
        ntu = conductance/capacity if capacity > 0 else 0.
        removed = remaining*(-math.expm1(-ntu)) if capacity > 0 else 0.
        after = remaining*math.exp(-ntu) if capacity > 0 else 0.
        rows.append(dict(conductance_W_K=conductance, inlet_heat_W=remaining,
                         heat_W=removed, outlet_heat_W=after))
        remaining=after
    total=math.fsum(row['heat_W'] for row in rows)
    return dict(inlet_heat_W=initial, delivered_heat_W=total, outlet_heat_W=remaining,
                outlet_superheat_K=remaining/capacity if capacity > 0 else None,
                residual_W=initial-total-remaining, sections=rows,
                scope='steady serial zero-C contact; no flow/storage inferred at zero discharge')
