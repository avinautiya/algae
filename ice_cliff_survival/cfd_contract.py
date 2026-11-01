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
_FIELDS = ("rho_i", "rho_w", "t0", "t1", "q_wall_W", "wetted_area_m2", "melt_volume_ice_m3", "meltwater_volume_m3", "sensible_in_J", "sensible_out_J", "storage_change_J", "dissipation_J")
def _check_tol(**tols):
    for k, v in tols.items():
        if not (isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v) and v >= 0.0): raise ValueError(f"tolerance {k} must be finite and >= 0, got {v!r}")
@dataclass(frozen=True)
class ExchangeRecord:
    """Signs: q_wall_W >0 = heat delivered to the ice (melting); it MAY be negative (refreezing/heat extraction) in a diagnostic record, but then melt volume must match
    the signed latent change only if melt_volume_ice_m3 >= 0 is kept; here melt volumes/areas are >= 0 so a record with q<0 fails the latent check (refreezing is unsupported).
    storage_change_J, sensible_in/out_J and dissipation_J are signed energies (storage may be negative); every field must be finite."""
    t0: float; t1: float; q_wall_W: float; wetted_area_m2: float; melt_volume_ice_m3: float; meltwater_volume_m3: float; sensible_in_J: float; sensible_out_J: float; storage_change_J: float; dissipation_J: float = 0.0; rho_i: float = RHO_I; rho_w: float = RHO_W      # densities used by the producing solver (equal-density fixtures set both)
    def check_finite(self):
        for k in _FIELDS:
            v = getattr(self, k)
            if not isinstance(v, (int, float)) or isinstance(v, bool) or not math.isfinite(v): raise ValueError(f"field {k} must be a finite number, got {v!r}")
    def validate(self, rtol=1e-6):
        """structural checks (finite, ordered, non-negative extensives, phase-change mass); returns DIAGNOSTIC latent-vs-heat residuals (no tolerance on them)."""
        _check_tol(rtol=rtol); self.check_finite()
        if not self.t1 > self.t0: raise ValueError("t1 must exceed t0")
        for k in ("wetted_area_m2", "melt_volume_ice_m3", "meltwater_volume_m3"): require(getattr(self, k), k, 0.0)
        require(self.rho_i, "rho_i", 0.0, lo_strict=True); require(self.rho_w, "rho_w", 0.0, lo_strict=True)
        if abs(self.meltwater_volume_m3 - self.rho_i / self.rho_w * self.melt_volume_ice_m3) > rtol * max(self.meltwater_volume_m3, 1e-30): raise ValueError("meltwater volume != rho_i/rho_w * ice volume (phase-change mass accounting)")
        return self.diagnostics()
    def heat_J(self): return self.q_wall_W * (self.t1 - self.t0)
    def diagnostics(self):
        heat = self.heat_J(); latent = self.rho_i * LF * self.melt_volume_ice_m3
        return dict(latent_minus_heat_J=latent - heat, latent_over_heat=(latent - heat) / max(abs(heat), abs(latent), 1e-300), water_energy_residual_J=self.water_energy_residual())
    def water_energy_residual(self): return self.sensible_in_J + self.dissipation_J - self.sensible_out_J - self.heat_J() - self.storage_change_J
def diagnose_records(records, rtol=1e-6):
    """DIAGNOSTIC pass: structural validation of each record (finite, ordered, contiguous) and the residuals; never applies an energy/latent tolerance. Returns a list of dicts."""
    out = []; prev = None
    for i, r in enumerate(records):
        d = r.validate(rtol)
        if prev is not None and abs(r.t0 - prev.t1) > 1e-9 * max(1.0, abs(r.t1)): raise ValueError(f"record {i}: gap/overlap in time")
        den = abs(r.sensible_in_J) + abs(r.sensible_out_J) + abs(r.heat_J()) + abs(r.storage_change_J) + abs(r.dissipation_J) + 1e-300
        d["water_energy_rel"] = d["water_energy_residual_J"] / den; out.append(d); prev = r
    return out
def validate_records(records, rtol=1e-6, latent_rtol=1e-3, energy_rtol=1e-6):
    """ACCEPTANCE GATE: diagnose_records() then enforce tolerances: |water energy residual| <= energy_rtol * sum|terms| and |latent-heat|/max(|heat|,|latent|) <= latent_rtol. Raises on any NaN/inf, bad tolerance, or violation."""
    _check_tol(rtol=rtol, latent_rtol=latent_rtol, energy_rtol=energy_rtol)
    for i, d in enumerate(diagnose_records(records, rtol)):
        if abs(d["water_energy_rel"]) > energy_rtol: raise ValueError(f"record {i}: water energy residual {d['water_energy_residual_J']:.3e} J exceeds {energy_rtol} of terms")
        if abs(d["latent_over_heat"]) > latent_rtol: raise ValueError(f"record {i}: latent/heat mismatch {d['latent_over_heat']:.3e}")
    return True
def parse_exchange_csv(path, validate=True):
    """columns: t0,t1,q_wall_W,wetted_area_m2,melt_volume_ice_m3,meltwater_volume_m3,sensible_in_J,sensible_out_J,storage_change_J[,dissipation_J]"""
    out = []
    with open(path) as f:
        for r in csv.DictReader(f): out.append(ExchangeRecord(**{k: float(v) for k, v in r.items()}))
    if validate: validate_records(out)
    return out
def _neumann_lambda(Ste):
    from math import erf, exp, sqrt, pi
    lo, hi = 1e-12, 8.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if mid * exp(mid * mid) * erf(mid) < Ste / sqrt(pi): lo = mid
        else: hi = mid
    return 0.5 * (lo + hi)
def neumann_fixture(Ste, n=200, t_end=3600.0, rho=RHO_W, cp=4217.0, k=0.561, Lf=LF, rho_w=None, rho_i=None):
    """analytic one-phase Stefan solution for an EQUAL-DENSITY ice/water pair (zero liquid velocity): s = 2 lam sqrt(kappa t), lam e^{lam^2} erf lam = Ste/sqrt(pi), Ste = cp Tw/Lf.
    M2 CORRECTION: the previous default (rho_w=1000, rho_i=917) paired the equal-density similarity solution with unequal densities, so emitted heat was rho_i/rho_w = 0.917 of its own analytic conductive heat. Now one density `rho` is used throughout;
    passing rho_w / rho_i that differ raises (a genuine unequal-density reference needs the density-jump similarity solution, see cfd/stefan2d.solve_lambda; not offered as a contract fixture).
    Records per unit area: ice volume lost dV=s(t1)-s(t0); heat from the ANALYTIC front conductive flux q_f(t)=k Tw e^{-lam^2}/(sqrt(pi kappa t) erf lam) integrated exactly (= rho Lf dV by the Stefan condition; see neumann_front_heat for the independent quadrature)."""
    from math import erf, exp, sqrt, pi
    for nm, v in (("rho_w", rho_w), ("rho_i", rho_i)):
        if v is not None and abs(v - rho) > 1e-12 * rho: raise ValueError(f"neumann_fixture is an equal-density solution: {nm}={v} != rho={rho}")
    require(Ste, "Ste", 0.0, lo_strict=True); require(t_end, "t_end", 0.0, lo_strict=True)
    if not (isinstance(n, int) and n >= 1): raise ValueError("n must be a positive int")
    lam = _neumann_lambda(Ste); kap = k / (rho * cp); Tw = Ste * Lf / cp; recs = []; ts = [t_end * (i + 1) / n for i in range(n)]
    s = lambda t: 2 * lam * sqrt(kap * t)
    for i, t1 in enumerate(ts):
        t0 = 0.0 if i == 0 else ts[i - 1]; dV = s(t1) - s(t0)
        heat = k * Tw * exp(-lam * lam) / (sqrt(pi * kap) * erf(lam)) * 2.0 * (sqrt(t1) - sqrt(t0))     # exact integral of the analytic front flux
        recs.append(ExchangeRecord(t0, t1, heat / (t1 - t0), 1.0, dV, dV, heat, 0.0, 0.0, 0.0, rho, rho))   # sensible_in = heat supplied at the hot wall; liquid-layer storage not modelled in this contract fixture
    return lam, Tw, recs
def volume_to_grain_inputs(rec: ExchangeRecord, c_mix):
    """S_eng over the step from the ice volume lost: grains released = c_mix/(1-c_mix) per pure-ice volume melted... here melt_volume_ice is PURE-ICE volume, so solid released = c_ice * V_ice with c_ice = c_mix/(1-c_mix)."""
    return c_mix / (1.0 - c_mix) * rec.melt_volume_ice_m3
