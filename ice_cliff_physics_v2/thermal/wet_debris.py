"""Exact, passive lumped wet-debris / bare-ice thermal contact.

A prescribed flow exchanges heat with coextensive bare and debris patches. Water
is quasi-steady and has finite advective heat capacity. Debris is a single
uniform-temperature reservoir; the attached ice is temperate (0 C). This is a
closure model, not resolved water/debris CFD or a field-validated transfer law.
"""
from __future__ import annotations
from dataclasses import dataclass, asdict
import math


class UnsupportedThermalState(ValueError):
    """Unsupported phase/domain; no silently clipped heat or temperature."""


def _number(value, name, nonnegative=False, positive=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(name + ' must be finite numeric')
    if (nonnegative and value < 0) or (positive and value <= 0):
        raise ValueError(name + ' outside supported sign range')
    return float(value)


def _phi(x):
    """(1-exp(-x))/x, including its exact zero limit."""
    if x == 0:
        return 1.
    if x < 1e-4:
        return 1. + x*(-.5 + x*(1./6 + x*(-1./24 + x/120)))
    return -math.expm1(-x)/x


def _one_minus_phi(x):
    if x < 1e-4:
        return x*(.5 + x*(-1./6 + x*(1./24 + x*(-1./120 + x/720))))
    return 1.-_phi(x)


def _psi(x):
    """(x + expm1(-x))/x**2, with stable zero limit."""
    if x < 1e-4:
        return .5 + x*(-1./6 + x*(1./24 + x*(-1./120 + x/720)))
    return (x+math.expm1(-x))/(x*x)


@dataclass(frozen=True)
class ThermalState:
    temperature_C: float

    def validate(self):
        _number(self.temperature_C, 'debris temperature', nonnegative=True)
        return self


@dataclass(frozen=True)
class ThermalParameters:
    capacity_J_K: float
    ice_conductance_W_K: float
    bare_water_conductance_W_K: float
    debris_water_conductance_W_K: float

    def validate(self):
        for k, v in asdict(self).items():
            _number(v, k, nonnegative=True)
        return self


@dataclass(frozen=True)
class ThermalResult:
    state: ThermalState
    integrated_debris_temperature_K_s: float
    mean_debris_temperature_C: float
    debris_storage_change_J: float
    atmospheric_heat_J: float
    water_inlet_heat_J: float
    water_outlet_heat_J: float
    water_to_debris_heat_J: float
    water_to_bare_heat_J: float
    covered_ice_heat_J: float
    bare_ice_heat_J: float
    total_ice_heat_J: float
    global_energy_residual_J: float
    water_exchange_residual_J: float
    passivity_margin_J: float
    outlet_mean_temperature_C: float | None
    zero_capacity_equilibrium: bool


def advance(state, parameters, dt_s, inlet_temperature_C, water_capacity_W_K,
            atmospheric_power_W=0., ice_temperature_C=0.):
    """Advance one constant-forcing interval exactly within this closure.

    Conductances and capacity are EXTENSIVE, not per-area. For actual flow Q,
    water_capacity_W_K = rho_w * cp_w * Q. Atmospheric power acts on the debris
    reservoir only; bare atmospheric heat is a separate caller-owned port.

    Outputs at the ice are energy available for latent melt of temperate ice;
    conversion to volume is the caller's responsibility. There is no ice-cold
    debt and no refreezing. Any trajectory crossing below zero is refused.
    Atmospheric power can be negative provided the entire solution remains
    nonnegative. Zero capacity is an explicit algebraic steady limit.
    """
    if not isinstance(state, ThermalState) or not isinstance(parameters, ThermalParameters):
        raise TypeError('ThermalState and ThermalParameters required')
    state.validate(); parameters.validate()
    dt = _number(dt_s, 'duration', positive=True)
    tw = _number(inlet_temperature_C, 'inlet temperature', nonnegative=True)
    flow = _number(water_capacity_W_K, 'water capacity', nonnegative=True)
    air = _number(atmospheric_power_W, 'atmospheric power')
    ti = _number(ice_temperature_C, 'ice temperature')
    if ti != 0:
        raise UnsupportedThermalState('temperate 0 C substrate required; cold ice/refreezing is not implemented')
    C, K, Gb, Gd = (parameters.capacity_J_K, parameters.ice_conductance_W_K,
                    parameters.bare_water_conductance_W_K, parameters.debris_water_conductance_W_K)
    t0 = state.temperature_C
    G = Gb+Gd
    if not math.isfinite(G):
        raise ValueError('conductance sum overflow')
    if flow == 0 or G == 0:
        # Zero flow is the caller's dry-contact state; there is no shared water.
        a = b = bare_in = bare_debris = 0.
        alpha = H = 0.
        water_decay = 1.
        water_decay_complement = 0.
    else:
        ntu = G/flow
        phi = _phi(ntu)
        alpha = Gd/G
        a = Gd*phi
        # Stable positive form of Gd*(1-(1-phi)*alpha).
        b = Gd*(Gb/G + phi*alpha)
        bare_in = Gb*phi
        bare_debris = Gb*_one_minus_phi(ntu)*alpha
        water_decay = math.exp(-ntu)
        water_decay_complement = -math.expm1(-ntu)
        H = flow*water_decay_complement
    beta = b+K
    forcing = air+a*tw
    if not math.isfinite(beta) or not math.isfinite(forcing):
        raise ValueError('derived thermal coefficients overflow')
    if C == 0:
        if beta == 0:
            if forcing != 0:
                raise UnsupportedThermalState('zero capacity and zero loss cannot absorb nonzero forcing')
            tend = t0
        else:
            tend = forcing/beta
        integral = tend*dt
    elif beta == 0:
        increment = (forcing/C)*dt
        tend = t0+increment
        integral = (t0+.5*increment)*dt
    else:
        x = (beta/C)*dt
        if x < 1e-4:
            derivative = (forcing-beta*t0)/C
            tend = t0+derivative*dt*_phi(x)
            integral = t0*dt+derivative*dt*dt*_psi(x)
        else:
            steady = forcing/beta
            decay = math.exp(-x)
            tend = steady+(t0-steady)*decay
            integral = steady*dt+(t0-steady)*(C/beta)*(-math.expm1(-x))
    if not all(math.isfinite(x) for x in [tend, integral]):
        raise ValueError('derived solution overflow')
    if tend < 0 or integral < 0:
        raise UnsupportedThermalState('debris crosses below 0 C; freezing/cold-substrate closure required')
    # Scalar constant-forcing solution is monotone: endpoint positivity certifies
    # positivity for the complete interval (zero-C limit is algebraic).
    qd = a*tw*dt-b*integral
    qb = bare_in*tw*dt+bare_debris*integral
    qi = K*integral
    storage = C*(tend-t0)
    water_in = flow*tw*dt
    water_exchange = qd+qb
    # Convex integral formula stays nonnegative even when tiny throughflow
    # heat is dwarfed by debris-to-bare redistribution. Do not subtract two
    # much larger contact transfers to recover this small outlet quantity.
    water_out = flow*(water_decay*tw*dt+water_decay_complement*alpha*integral)
    independent_water_exchange = H*(tw*dt-alpha*integral)
    residual = air*dt+independent_water_exchange-storage-qi-qb
    water_residual = water_exchange-independent_water_exchange
    # Negative atmospheric power removes energy; positive atmospheric power and
    # initial debris enthalpy are the only non-water resources admitted here.
    available = C*t0+water_in+max(air, 0.)*dt
    ice = qi+qb
    margin = available-ice
    values = [qd,qb,qi,storage,water_in,water_out,residual,water_residual,margin]
    if not all(math.isfinite(x) for x in values):
        raise ValueError('derived heat inventory overflow')
    # Do not clamp small roundoff; expose signed residual and passivity margin.
    scale=max(1.,abs(water_in),abs(water_out),abs(air*dt),abs(C*t0),abs(ice))
    if water_out < -1e-12*scale or margin < -1e-12*scale:
        raise ArithmeticError('passivity failure')
    return ThermalResult(ThermalState(tend), integral, integral/dt, storage, air*dt,
        water_in, water_out, qd, qb, qi, qb, ice, residual, water_residual, margin,
        water_out/(flow*dt) if flow>0 else None, C==0)
