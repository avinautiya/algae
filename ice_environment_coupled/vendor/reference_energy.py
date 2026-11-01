"""Fixed-area, per-unit-area energy closure; SI units throughout.

This is a lumped cold-content reservoir, not a cold-ice temperature PDE.
Meltwater leaves immediately; refreezing and changing column area are unsupported.
"""
from dataclasses import dataclass
import math


def _finite(name, value, minimum=None):
    if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or (minimum is not None and value < minimum):
        raise ValueError(f"invalid {name}: {value}")


@dataclass(frozen=True)
class ColumnParameters:
    debris_capacity_j_m2_k: float
    conductance_w_m2_k: float
    latent_heat_j_kg: float

    def __post_init__(self):
        for name, value in vars(self).items():
            _finite(name, value, 0)
        if self.debris_capacity_j_m2_k == 0 or self.latent_heat_j_kg == 0:
            raise ValueError("capacity and latent heat must be positive")


@dataclass(frozen=True)
class ColumnState:
    debris_temperature_c: float
    cold_reserve_j_m2: float
    melted_kg_m2: float = 0.0

    def __post_init__(self):
        _finite("temperature", self.debris_temperature_c)
        _finite("cold reserve", self.cold_reserve_j_m2, 0)
        _finite("melt", self.melted_kg_m2, 0)


@dataclass(frozen=True)
class BareIceState:
    cold_reserve_j_m2: float
    melted_kg_m2: float = 0.0

    def __post_init__(self):
        _finite("cold reserve", self.cold_reserve_j_m2, 0)
        _finite("melt", self.melted_kg_m2, 0)


@dataclass(frozen=True)
class StepResult:
    state: object  # ColumnState or BareIceState
    atmospheric_j_m2: float
    water_j_m2: float
    debris_storage_j_m2: float
    cold_storage_j_m2: float
    latent_j_m2: float
    residual_j_m2: float


def step_bare_ice(state, latent_heat_j_kg, dt_s, atmospheric_w_m2,
                  water_w_m2, *, refreeze=False):
    """Direct fixed-area bare-ice ledger, with no debris temperature/storage.

    Both supplied fluxes act directly on lumped ice enthalpy. Cooling builds
    cold reserve, never refreezes departed meltwater. No surface temperature
    or subfreezing temperature profile is diagnosed by this closure.
    """
    if refreeze:
        raise NotImplementedError("refreezing is unsupported")
    if not isinstance(state, BareIceState):
        raise TypeError("bare ice requires BareIceState")
    _finite("latent heat", latent_heat_j_kg, 0)
    if latent_heat_j_kg == 0:
        raise ValueError("latent heat must be positive")
    _finite("dt", dt_s, 0)
    _finite("atmospheric flux", atmospheric_w_m2)
    _finite("water flux", water_w_m2, 0)
    atmosphere, water = dt_s * atmospheric_w_m2, dt_s * water_w_m2
    remaining = state.cold_reserve_j_m2 - atmosphere - water
    cold, latent = max(remaining, 0.0), max(-remaining, 0.0)
    new = BareIceState(cold, state.melted_kg_m2 + latent / latent_heat_j_kg)
    cold_storage = state.cold_reserve_j_m2 - cold
    return StepResult(new, atmosphere, water, 0.0, cold_storage, latent,
                      atmosphere + water - cold_storage - latent)


def step_column(state, parameters, dt_s, atmospheric_w_m2, water_w_m2,
                *, refreeze=False):
    """Backward-Euler debris storage coupled to a zero-C reference interface.

    Signed atmospheric heat enters debris; nonnegative water heat enters ice.
    The interface uses a *zero-C reference*, even with cold content outstanding:
    this closure cannot resolve subfreezing ice temperature or gradients.
    Conductive energy plus water heat first repays cold content, then melts ice.
    Negative interface heat increases cold content; it never freezes departed
    meltwater. No finite ice mass/exhaustion or evolving-area accounting is made.
    """
    if refreeze:
        raise NotImplementedError("refreezing is unsupported")
    _finite("dt", dt_s, 0)
    _finite("atmospheric flux", atmospheric_w_m2)
    _finite("water flux", water_w_m2, 0)
    c = parameters.debris_capacity_j_m2_k
    g = parameters.conductance_w_m2_k
    temperature = (c * state.debris_temperature_c + dt_s * atmospheric_w_m2) / (c + dt_s * g)
    delivered = dt_s * (g * temperature + water_w_m2)
    remaining = state.cold_reserve_j_m2 - delivered
    cold = max(remaining, 0.0)
    latent = max(-remaining, 0.0)
    new = ColumnState(temperature, cold, state.melted_kg_m2 + latent / parameters.latent_heat_j_kg)
    atmospheric = dt_s * atmospheric_w_m2
    water = dt_s * water_w_m2
    debris_storage = c * (temperature - state.debris_temperature_c)
    cold_storage = state.cold_reserve_j_m2 - cold
    residual = atmospheric + water - debris_storage - cold_storage - latent
    return StepResult(new, atmospheric, water, debris_storage, cold_storage, latent, residual)


def atmospheric_flux(*, face_shortwave_w_m2, incoming_longwave_w_m2,
                     surface_temperature_c, air_temperature_c,
                     air_specific_humidity_kg_kg, wind_speed_m_s, pressure_pa,
                     albedo, emissivity, sensible_exchange_coefficient,
                     moisture_exchange_coefficient, air_heat_capacity_j_kg_k,
                     latent_vaporization_j_kg, dry_air_gas_constant_j_kg_k,
                     stefan_boltzmann_w_m2_k4):
    """Return signed net W/m2 and components, positive toward the surface.

    Required shortwave is already projected/shaded onto the face. Neutral bulk
    transfer: H=rho cp Ch U(Ta-Ts), LE=rho Lv Ce U(qa-qs).
    Density uses dry ideal gas; qs uses liquid-water Tetens saturation pressure
    (610.78 exp(17.27 T/(T+237.3)) Pa) and epsilon=0.622. No stability correction,
    rain heat, terrain radiation, or sublimation is included. Surface T must be
    0..50 C; evaporation consumes latent heat, condensation supplies it.
    """
    values = locals().copy()
    for name, value in values.items():
        _finite(name, value)
    if not 0 <= surface_temperature_c <= 50 or not -50 <= air_temperature_c <= 50:
        raise ValueError("temperature outside closure range")
    for name in ("albedo", "emissivity", "air_specific_humidity_kg_kg"):
        if not 0 <= values[name] <= 1:
            raise ValueError(f"{name} must be in [0,1]")
    for name in ("face_shortwave_w_m2", "incoming_longwave_w_m2", "wind_speed_m_s",
                 "sensible_exchange_coefficient", "moisture_exchange_coefficient"):
        _finite(name, values[name], 0)
    for name in ("pressure_pa", "air_heat_capacity_j_kg_k", "latent_vaporization_j_kg",
                 "dry_air_gas_constant_j_kg_k", "stefan_boltzmann_w_m2_k4"):
        if values[name] <= 0:
            raise ValueError(f"{name} must be positive")
    ts = surface_temperature_c
    es = 610.78 * math.exp(17.27 * ts / (ts + 237.3))
    if es >= pressure_pa:
        raise ValueError("saturation vapor pressure must be below total pressure")
    qs = 0.622 * es / (pressure_pa - 0.378 * es)
    rho = pressure_pa / (dry_air_gas_constant_j_kg_k * (air_temperature_c + 273.15))
    parts = {
        "shortwave": (1 - albedo) * face_shortwave_w_m2,
        "longwave": emissivity * (incoming_longwave_w_m2 - stefan_boltzmann_w_m2_k4 * (ts + 273.15)**4),
        "sensible": rho * air_heat_capacity_j_kg_k * sensible_exchange_coefficient * wind_speed_m_s * (air_temperature_c - ts),
        "latent": rho * latent_vaporization_j_kg * moisture_exchange_coefficient * wind_speed_m_s * (air_specific_humidity_kg_kg - qs),
    }
    return {**parts, "net": sum(parts.values())}
