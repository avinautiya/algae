"""SI-unit lake storage and conservative exchange with an external channel.

The implicit example step assumes a nonnegative, monotone stage-discharge law.
General CFD coupling should exchange integrated boundary fluxes via PortTransfer;
it must not impose both a lake stage and an unrelated discharge at the same port.
"""
from dataclasses import dataclass
from bisect import bisect_right
from math import expm1, isfinite, sqrt
from typing import Callable, Tuple


class DomainError(ValueError):
    """The requested state is outside the declared physical/numerical domain."""


def finite(value, name, minimum=None, strictly=False):
    if not isfinite(value):
        raise DomainError(f"{name} must be finite")
    if minimum is not None and (value < minimum or (strictly and value == minimum)):
        raise DomainError(f"{name} must be {'>' if strictly else '>='} {minimum}")
    return value


@dataclass(frozen=True)
class Constants:
    rho_water: float = 1000.0
    rho_ice: float = 917.0
    cp_water: float = 4217.0
    latent_heat: float = 334000.0
    gravity: float = 9.81

    def __post_init__(self):
        for key, value in vars(self).items():
            finite(value, key, 0.0, strictly=True)

    @property
    def rc(self):
        return self.rho_water * self.cp_water


@dataclass(frozen=True)
class StageAreaTable:
    """Piecewise-linear wetted surface area, exactly integrated for storage.

The first stage is the zero-volume floor. All areas must be positive. Tables
with a point-like zero-area bottom need explicit preprocessing, not silent
regularization. No stage or volume extrapolation is allowed.
"""
    stages_m: Tuple[float, ...]
    areas_m2: Tuple[float, ...]

    def __post_init__(self):
        object.__setattr__(self, "stages_m", tuple(self.stages_m))
        object.__setattr__(self, "areas_m2", tuple(self.areas_m2))
        if len(self.stages_m) < 2 or len(self.stages_m) != len(self.areas_m2):
            raise DomainError("Need at least two matched stage/area entries")
        for h in self.stages_m:
            finite(h, "stage")
        for a in self.areas_m2:
            finite(a, "area", 0.0, strictly=True)
        if any(b <= a for a, b in zip(self.stages_m, self.stages_m[1:])):
            raise DomainError("Stages must strictly increase")

    @property
    def cumulative_volumes(self):
        values = [0.0]
        for i in range(len(self.stages_m) - 1):
            values.append(values[-1] + 0.5 * (self.areas_m2[i] + self.areas_m2[i+1])
                          * (self.stages_m[i+1] - self.stages_m[i]))
        return tuple(values)

    @property
    def capacity_m3(self):
        return self.cumulative_volumes[-1]

    def _stage_segment(self, stage):
        finite(stage, "stage")
        if not self.stages_m[0] <= stage <= self.stages_m[-1]:
            raise DomainError("Stage outside measured storage table")
        return min(bisect_right(self.stages_m, stage) - 1, len(self.stages_m) - 2)

    def area(self, stage):
        i = self._stage_segment(stage)
        f = (stage - self.stages_m[i]) / (self.stages_m[i+1] - self.stages_m[i])
        return self.areas_m2[i] + f * (self.areas_m2[i+1] - self.areas_m2[i])

    def volume(self, stage):
        i = self._stage_segment(stage)
        dh = stage - self.stages_m[i]
        return self.cumulative_volumes[i] + 0.5 * (self.areas_m2[i] + self.area(stage)) * dh

    def stage(self, volume):
        finite(volume, "volume", 0.0)
        cv = self.cumulative_volumes
        if volume > cv[-1]:
            raise DomainError("Volume exceeds storage-table capacity")
        i = min(bisect_right(cv, volume) - 1, len(cv) - 2)
        dv = volume - cv[i]
        a = self.areas_m2[i]
        slope = (self.areas_m2[i+1] - a) / (self.stages_m[i+1] - self.stages_m[i])
        # Stable positive root of dv = a*dh + slope*dh^2/2; also works at slope=0.
        dh = 2.0 * dv / (a + sqrt(a*a + 2.0*slope*dv))
        return self.stages_m[i] + dh


@dataclass(frozen=True)
class LakeState:
    volume_m3: float
    sensible_energy_J: float

    def __post_init__(self):
        finite(self.volume_m3, "lake volume", 0.0)
        finite(self.sensible_energy_J, "sensible energy", 0.0)
        if self.volume_m3 == 0.0 and self.sensible_energy_J != 0.0:
            raise DomainError("An empty lake cannot store sensible water heat")

    def temperature(self, constants=Constants()):
        return self.sensible_energy_J / (constants.rc * self.volume_m3) if self.volume_m3 else 0.0

    @classmethod
    def from_temperature(cls, volume, temperature, constants=Constants()):
        finite(temperature, "temperature excess above melting", 0.0)
        return cls(volume, constants.rc * volume * temperature)


@dataclass(frozen=True)
class LakeForcing:
    inflow_m3_s: float = 0.0
    inflow_temperature_K: float = 0.0
    external_power_W: float = 0.0

    def __post_init__(self):
        finite(self.inflow_m3_s, "inflow", 0.0)
        finite(self.inflow_temperature_K, "inflow temperature", 0.0)
        finite(self.external_power_W, "external power")


@dataclass(frozen=True)
class PortTransfer:
    """Time-integrated advective flux; positive from lake into channel.

Use the SAME object with opposite signs on the two domains. Negative values
represent backflow; its heat must be computed from the channel donor state.
Thermal conduction across a zero-flow port belongs in a separate heat ledger.
"""
    volume_m3: float
    sensible_energy_J: float

    def __post_init__(self):
        finite(self.volume_m3, "port volume")
        finite(self.sensible_energy_J, "port energy")
        if self.volume_m3 * self.sensible_energy_J < 0.0:
            raise DomainError("Advective volume and sensible heat must have compatible signs")
        if self.volume_m3 == 0.0 and self.sensible_energy_J != 0.0:
            raise DomainError("A zero-flow advective port cannot transport sensible heat")


def exchange_port(lake, transfer):
    """Apply one integrated CFD exchange without clipping or silently limiting it.

The caller owns timestep rejection if a flux removes unavailable mass/heat.
Do not apply this a second time after advance_lake, which already applies it.
"""
    return LakeState(lake.volume_m3 - transfer.volume_m3,
                     lake.sensible_energy_J - transfer.sensible_energy_J)


@dataclass(frozen=True)
class LakeStep:
    state: LakeState
    port: PortTransfer
    inflow_volume_m3: float
    inflow_energy_J: float
    external_energy_J: float
    water_residual_m3: float
    energy_residual_J: float


def advance_lake(storage, state, forcing, dt_s, discharge: Callable[[float], float], constants=Constants()):
    """Backward-Euler well-mixed lake; fixed monotone rating during one step.

Solves V1 + dt*Q(H(V1)) = V0 + dt*I. Positivity is structural, without an
arbitrary wet/dry floor. This is first-order; timestep convergence is required.
The rating must give zero discharge at the basin floor and be nondecreasing.
This routine is NOT an implicit coupling solver for an unsteady CFD domain.
"""
    finite(dt_s, "dt", 0.0, strictly=True)
    storage.stage(state.volume_m3)
    def q(v):
        return finite(discharge(storage.stage(v)), "discharge", 0.0)
    if q(0.0) != 0.0:
        raise DomainError("Rating drains at the zero-storage floor")
    vin = forcing.inflow_m3_s * dt_s
    ein = constants.rc * vin * forcing.inflow_temperature_K
    eatm = forcing.external_power_W * dt_s
    available = state.volume_m3 + vin
    energy = state.sensible_energy_J + ein + eatm
    if energy < 0.0:
        raise DomainError("Cooling would freeze water; phase change in the lake is out of scope")
    if available == 0.0:
        if energy != 0.0:
            raise DomainError("External heating of a dry lake is not represented")
        return LakeStep(state, PortTransfer(0.0, 0.0), vin, ein, eatm, 0.0, 0.0)
    lo, hi = 0.0, min(available, storage.capacity_m3)
    if hi + dt_s*q(hi) < available:
        raise DomainError("Lake exceeds storage table; extend measured geometry or model overflow")
    if hi == available and q(hi) == 0.0:
        v1 = available  # exact closed-outlet storage; avoid cumulative root roundoff
    else:
        # Absolute+relative resolution in volume, with a bounded iteration count.
        tolerance = 1e-13 * max(1.0, available)
        for _ in range(100):
            mid = 0.5 * (lo + hi)
            if mid + dt_s*q(mid) > available:
                hi = mid
            else:
                lo = mid
            if hi - lo <= tolerance:
                break
        v1 = 0.5 * (lo + hi)
    out = dt_s*q(v1)
    total = v1 + out
    e1 = energy * (v1 / total)
    eout = energy * (out / total)
    result = LakeState(v1, e1)
    return LakeStep(result, PortTransfer(out, eout), vin, ein, eatm,
                    v1 + out - available, e1 + eout - energy)


def critical_discharge(stage_m, sill_m, width_m, coefficient=1.0, constants=Constants()):
    """Ideal rectangular critical control: Q=Cd*b*sqrt(g)*(2*head/3)^1.5.

Negligible approach velocity, unsubmerged control, no downstream backwater.
This is a demonstration rating, not a prediction from resolved CFD.
"""
    finite(stage_m, "stage")
    finite(sill_m, "sill")
    finite(width_m, "width", 0.0, strictly=True)
    finite(coefficient, "discharge coefficient", 0.0)
    if coefficient > 1.0:
        raise DomainError("Demo discharge coefficient must be within [0, 1]")
    head = max(stage_m - sill_m, 0.0)
    return coefficient * width_m * sqrt(constants.gravity) * (2.0*head/3.0)**1.5


@dataclass(frozen=True)
class ControlHeat:
    latent_energy_J: float
    outflow: PortTransfer
    ice_loss_m3: float
    meltwater_m3: float
    lowering_m: float


def thermal_control(port, dt_s, contact_area_m2, htc_W_m2_K, melting=True, constants=Constants()):
    """Steady plug-flow heat-exchange illustration for a temperate ice control.

No transient channel storage, freezing, friction heating or turbulent closure
prediction. The prescribed heat-transfer coefficient is an uncertain input.
The same heat loss enters latent melting; new meltwater is at melting point.
"""
    finite(dt_s, "dt", 0.0, strictly=True)
    finite(contact_area_m2, "contact area", 0.0, strictly=True)
    finite(htc_W_m2_K, "heat-transfer coefficient", 0.0)
    if port.volume_m3 < 0.0:
        raise DomainError("Demo control supports forward flow only; generic port supports backflow")
    heat = 0.0
    if melting and port.volume_m3 > 0.0:
        exponent = htc_W_m2_K * contact_area_m2 * dt_s / (constants.rc * port.volume_m3)
        heat = port.sensible_energy_J * (-expm1(-exponent))
    ice = heat / (constants.rho_ice * constants.latent_heat)
    meltwater = heat / (constants.rho_water * constants.latent_heat)
    return ControlHeat(heat, PortTransfer(port.volume_m3 + meltwater,
                       port.sensible_energy_J - heat), ice, meltwater, ice/contact_area_m2)
