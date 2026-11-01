"""Integration contract between a future coupled flow/heat/moving-ice solver and the grain/exposure layer: data structures, validation, flux export parsing, budget adapters, analytical fixtures.
NONE of this verifies a coupled solver. It fixes the numbers that must cross each interface and checks their internal consistency.
Per-step exchange record (all extensive, integrated over the reporting surface and the step; SI):
  from FLOW/HEAT to ICE      : q_wall_W        heat delivered to the ice interface (>= 0 means melting energy supply)             [W]  (same value used for water cooling)
                               wetted_area_m2  ice area in contact with flowing water                                              [m2]
  from ICE to FLOW           : melt_volume_ice_m3 (ice-volume lost), meltwater_volume_m3 = rho_i/rho_w * melt_volume_ice            [m3]
                               interface_speed_n m/s per patch (v_I.n, n from ice into water)
  from FLOW to GRAINS        : near-toe velocity v [m/s] and depth h [m] (clast-scale force balance inputs); discharge Q [m3/s]
  from GRAINS to FLOW/ICE    : covered_area_m2 (no stream-ice heat exchange there), roughness/blockage (not implemented: UnsupportedRun)
  from ICE to GRAINS         : S_eng = c_mix A r_n  (frozen grains released), face geometry
Time-step ledgers required from the solver: sensible_in_J, sensible_out_J, heat_to_ice_J, storage_change_J, dissipation_J (counted once), latent_J = rho_i Lf melt_volume_ice (must equal heat_to_ice_J up to the declared residual)."""
import csv, math
from dataclasses import dataclass
from .units import RHO_I, RHO_W, LF, require
@dataclass(frozen=True)
class ExchangeRecord:
    t0: float; t1: float; q_wall_W: float; wetted_area_m2: float; melt_volume_ice_m3: float; meltwater_volume_m3: float; sensible_in_J: float; sensible_out_J: float; storage_change_J: float; dissipation_J: float = 0.0
    def validate(self, rtol=1e-6):
        if not self.t1 > self.t0: raise ValueError("t1 must exceed t0")
        for k in ("wetted_area_m2", "melt_volume_ice_m3", "meltwater_volume_m3"): require(getattr(self, k), k, 0.0)
        if abs(self.meltwater_volume_m3 - RHO_I / RHO_W * self.melt_volume_ice_m3) > rtol * max(self.meltwater_volume_m3, 1e-30): raise ValueError("meltwater volume != rho_i/rho_w * ice volume (phase-change mass accounting)")
        heat = self.q_wall_W * (self.t1 - self.t0); latent = RHO_I * LF * self.melt_volume_ice_m3
        return dict(latent_minus_heat_J=latent - heat, latent_over_heat=(latent - heat) / max(abs(heat), 1e-300))
    def water_energy_residual(self): return self.sensible_in_J + self.dissipation_J - self.sensible_out_J - self.q_wall_W * (self.t1 - self.t0) - self.storage_change_J
def parse_exchange_csv(path):
    """columns: t0,t1,q_wall_W,wetted_area_m2,melt_volume_ice_m3,meltwater_volume_m3,sensible_in_J,sensible_out_J,storage_change_J[,dissipation_J]"""
    out = []
    with open(path) as f:
        for r in csv.DictReader(f): out.append(ExchangeRecord(**{k: float(v) for k, v in r.items()}))
    return out
def neumann_fixture(Ste, n=200, t_end=3600.0, rho_w=RHO_W, rho_i=RHO_I, cp=4217.0, k=0.561, Lf=LF):
    """analytic one-phase Stefan solution (equal-density limit, liquid velocity 0) as exchange records per unit area: s = 2 lam sqrt(kappa t), lam e^{lam^2} erf lam = Ste/sqrt(pi).
    Independent reference for adapters. Uses the exact similarity solution, not the discrete scheme."""
    from math import erf, exp, sqrt, pi
    lo, hi = 1e-9, 5.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if mid * exp(mid * mid) * erf(mid) < Ste / sqrt(pi): lo = mid
        else: hi = mid
    lam = 0.5 * (lo + hi); kap = k / (rho_w * cp); Tw = Ste * Lf / cp; recs = []; ts = [t_end * (i + 1) / n for i in range(n)]
    s = lambda t: 2 * lam * sqrt(kap * t)
    for i, t1 in enumerate(ts):
        t0 = 0.0 if i == 0 else ts[i - 1]; dV = s(t1) - s(t0); q = rho_i * Lf * dV / (t1 - t0)
        recs.append(ExchangeRecord(t0, t1, q, 1.0, dV, rho_i / rho_w * dV, 0.0, 0.0, 0.0))
    return lam, Tw, recs
def volume_to_grain_inputs(rec: ExchangeRecord, c_mix):
    """S_eng over the step from the ice volume lost: grains released = c_mix/(1-c_mix) per pure-ice volume melted... here melt_volume_ice is PURE-ICE volume, so solid released = c_ice * V_ice with c_ice = c_mix/(1-c_mix)."""
    return c_mix / (1.0 - c_mix) * rec.melt_volume_ice_m3
