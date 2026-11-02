"""Conservative lumped zero-C interface closure; SI units, no cold-ice PDE."""
from dataclasses import dataclass
import math


def _number(name, value, minimum=None):
    if (isinstance(value, bool) or not isinstance(value, (int, float))
            or not math.isfinite(value) or (minimum is not None and value < minimum)):
        raise ValueError(f"invalid {name}: {value}")


@dataclass(frozen=True)
class Parameters:
    capacity_j_m2_k: float
    conductance_w_m2_k: float
    latent_heat_j_kg: float

    def __post_init__(self):
        for name, value in vars(self).items():
            _number(name, value, 0)
        if self.capacity_j_m2_k == 0 or self.latent_heat_j_kg == 0:
            raise ValueError("capacity and latent heat must be positive")


@dataclass(frozen=True)
class State:
    area_m2: float
    debris_energy_j: float
    cold_debt_j: float
    cumulative_latent_j: float = 0.0

    def __post_init__(self):
        _number("area", self.area_m2, 0)
        if self.area_m2 == 0:
            raise ValueError("area must be positive")
        _number("debris energy", self.debris_energy_j)
        _number("cold debt", self.cold_debt_j, 0)
        _number("cumulative latent", self.cumulative_latent_j, 0)

    @property
    def accounted_energy_j(self):
        return self.debris_energy_j - self.cold_debt_j + self.cumulative_latent_j


@dataclass(frozen=True)
class AreaContent:
    """Declared material carried per m2; latent is material-associated history."""
    debris_energy_j_m2: float
    cold_debt_j_m2: float
    cumulative_latent_j_m2: float

    def __post_init__(self):
        _number("debris content", self.debris_energy_j_m2)
        _number("cold content", self.cold_debt_j_m2, 0)
        _number("latent content", self.cumulative_latent_j_m2, 0)


@dataclass(frozen=True)
class RemapResult:
    state: State
    debris_transfer_j: float
    cold_debt_transfer_j: float
    latent_transfer_j: float
    residual_j: float


def remap_area(state, new_area_m2, *, incoming=None, outgoing=None):
    """Change area only by explicitly declared signed material transfers.

    Transfers are positive into the domain. Outgoing content must describe the
    removed material, not necessarily the domain average. Reject removal of
    more cold debt or latent history than exists. Debris enthalpy is signed.
    """
    if not isinstance(state, State):
        raise TypeError("state must be State")
    _number("new area", new_area_m2, 0)
    if new_area_m2 == 0:
        raise ValueError("new area must be positive")
    da = new_area_m2 - state.area_m2
    if da == 0:
        if incoming is not None or outgoing is not None:
            raise ValueError("unchanged area takes no transfer")
        return RemapResult(state, 0., 0., 0., 0.)
    content = incoming if da > 0 else outgoing
    if not isinstance(content, AreaContent) or (outgoing if da > 0 else incoming) is not None:
        raise ValueError("declare only incoming content for expansion or outgoing for contraction")
    de, dc, dl = (da * x for x in vars(content).values())
    new = State(new_area_m2, state.debris_energy_j + de,
                state.cold_debt_j + dc, state.cumulative_latent_j + dl)
    residual = (new.accounted_energy_j - state.accounted_energy_j) - (de - dc + dl)
    _number("remap residual", residual)
    return RemapResult(new, de, dc, dl, residual)


@dataclass(frozen=True)
class StepResult:
    state: State
    atmospheric_j: float
    water_j: float
    debris_storage_j: float
    cold_storage_j: float
    latent_j: float
    residual_j: float


def step_bare(state, dt_s, atmospheric_w_m2, water_w_m2):
    """Direct signed ports into bare ice; debris energy must be exactly zero."""
    if not isinstance(state, State):
        raise TypeError("expected State")
    if state.debris_energy_j != 0:
        raise ValueError("bare ice cannot carry debris energy")
    _number("dt", dt_s, 0)
    _number("atmospheric flux", atmospheric_w_m2)
    _number("water flux", water_w_m2)
    atmosphere = dt_s * atmospheric_w_m2 * state.area_m2
    water = dt_s * water_w_m2 * state.area_m2
    _number("atmospheric energy", atmosphere)
    _number("water energy", water)
    remaining = state.cold_debt_j - atmosphere - water
    cold, latent = max(remaining, 0.), max(-remaining, 0.)
    new = State(state.area_m2, 0., cold, state.cumulative_latent_j + latent)
    cold_storage = state.cold_debt_j - cold
    residual = atmosphere + water - cold_storage - latent
    _number("bare residual", residual)
    return StepResult(new, atmosphere, water, 0., cold_storage, latent, residual)


def step(state, parameters, dt_s, atmospheric_w_m2, water_w_m2):
    """Backward Euler, matching the fixed-area reference for positive water.

    Both external fluxes are signed positive inward. Atmosphere heats debris;
    water acts on the ice reservoir. Delivered heat repays debt before melt.
    Negative delivered heat grows debt, never reverses past melt. Parameters
    stay fixed during this step; area changes require a separate remap.
    """
    if not isinstance(state, State) or not isinstance(parameters, Parameters):
        raise TypeError("expected State and Parameters")
    _number("dt", dt_s, 0)
    _number("atmospheric flux", atmospheric_w_m2)
    _number("water flux", water_w_m2)
    if dt_s == 0:
        return StepResult(state, 0., 0., 0., 0., 0., 0.)
    a = state.area_m2
    c = parameters.capacity_j_m2_k
    g = parameters.conductance_w_m2_k
    atmosphere = dt_s * atmospheric_w_m2 * a
    water = dt_s * water_w_m2 * a
    _number("atmospheric energy", atmosphere)
    _number("water energy", water)
    denominator = 1.0 + dt_s * g / c
    _number("implicit denominator", denominator)
    debris = (state.debris_energy_j + atmosphere) / denominator
    _number("debris energy", debris)
    storage = debris - state.debris_energy_j
    delivered = atmosphere - storage + water
    remaining = state.cold_debt_j - delivered
    cold, latent = max(remaining, 0.), max(-remaining, 0.)
    new = State(a, debris, cold, state.cumulative_latent_j + latent)
    cold_storage = state.cold_debt_j - cold
    residual = atmosphere + water - storage - cold_storage - latent
    _number("step residual", residual)
    return StepResult(new, atmosphere, water, storage, cold_storage, latent, residual)
