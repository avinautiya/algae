"""Section hydraulics with explicit observation provenance.
Pressure head p/(rho g) = eta - z_s equals full depth only when sensor clearance above the bed is KNOWN: h = p/(rho g) + (z_s - z_b). Unknown clearance -> depth is NOT determined (refusal, or an explicitly labelled pressure PROXY).
Area, wetted perimeter, hydraulic diameter, velocity and Reynolds number are recomputed from the section's own geometry for every call; a station's Re can never be taken from another station's series.
"""
import math
from dataclasses import dataclass
NU = 1.8e-6; G = 9.81
class HydraulicRefusal(ValueError): pass
DEPTH_BASES = ("measured_bed_clearance", "pressure_proxy_zero_clearance_assumed", "sourced_assumption", "critical_control_verified")
Q_BASES = ("measured", "reconstructed", "imposed_scenario")
def depth_from_pressure(pressure_head_m, sensor_clearance_m=None, allow_proxy=False):
    if not (isinstance(pressure_head_m, (int, float)) and math.isfinite(pressure_head_m) and pressure_head_m > 0): raise HydraulicRefusal("positive finite pressure head required (missing/purged records are not zero depth)")
    if sensor_clearance_m is None:
        if not allow_proxy: raise HydraulicRefusal("sensor clearance unknown: pressure head is submergence, not bed clearance; depth undetermined")
        return pressure_head_m, "pressure_proxy_zero_clearance_assumed"
    if not (math.isfinite(sensor_clearance_m) and sensor_clearance_m >= 0): raise HydraulicRefusal("clearance must be finite and >= 0")
    return pressure_head_m + sensor_clearance_m, "measured_bed_clearance"
@dataclass(frozen=True)
class SectionHydraulics:
    station: str; width_m: float; depth_m: float; discharge_m3_s: float; depth_basis: str; discharge_basis: str; geometry_source: str
    def __post_init__(self):
        if self.depth_basis not in DEPTH_BASES or self.discharge_basis not in Q_BASES: raise HydraulicRefusal("depth/discharge basis must be labelled explicitly")
        for v, n in ((self.width_m, "width"), (self.depth_m, "depth"), (self.discharge_m3_s, "discharge")):
            if not (isinstance(v, (int, float)) and math.isfinite(v) and v > 0): raise HydraulicRefusal(n + " must be positive finite")
    def derived(self):
        a = self.width_m * self.depth_m; p = self.width_m + 2 * self.depth_m; dh = 4 * a / p; v = self.discharge_m3_s / a
        return dict(station=self.station, area_m2=a, wetted_perimeter_m=p, hydraulic_diameter_m=dh, velocity_m_s=v, Re=v * dh / NU, Re_alt=4 * self.discharge_m3_s / (NU * p), Fr=v / math.sqrt(G * self.depth_m),
                    labels=dict(depth=self.depth_basis, discharge=self.discharge_basis, geometry=self.geometry_source, reconstructed=("reconstructed" in self.discharge_basis) or ("proxy" in self.depth_basis)))
def require_station(section, expected_station):
    if section.station != expected_station: raise HydraulicRefusal(f"section belongs to {section.station}; {expected_station} requires its own geometry (no P3->P5 substitution)")
    return section
