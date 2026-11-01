"""Conservative toe-stock accounting (solid volume). dVs/dt = S_up + S_mantle + S_eng + S_ext - Q_out.
Each source is recorded under its own name; exclusivity is enforced: a debris parcel released from the upper mantle may not also appear as external rockfall (caller declares `mantle_released_m3` already contained in S_ext to be subtracted, or must pass S_ext exclusive of it).
The export used in a step is the ACTUAL flux (capacity clipped so the stock cannot go negative) and is what the ledger records; clipping without recording the changed flux is not permitted here."""
from dataclasses import dataclass, field, asdict
from .units import require
SOURCES = ("S_up", "S_mantle", "S_eng", "S_ext")
@dataclass
class ToeStock:
    Vs: float = 0.0                                   # solid volume in the toe control volume [m3]
    t: float = 0.0
    cum: dict = field(default_factory=lambda: {k: 0.0 for k in SOURCES + ("Q_out_actual", "Q_out_capacity_unused", "roundoff_correction")})
    def __post_init__(self): require(self.Vs, "Vs", 0.0); require(self.t, "t")
    def step(self, dt, sources, capacity, mantle_double_counted_in_ext=0.0):
        """sources: dict of the four source rates [m3/s]; capacity: admissible export rate [m3/s] (>= 0, from a declared transport closure);
        mantle_double_counted_in_ext: rate [m3/s] of mantle material the caller has ALSO included in S_ext; it is removed from S_ext (and reported) so nothing is counted twice."""
        dt = require(dt, "dt", 0.0, lo_strict=True); unknown = set(sources) - set(SOURCES)
        if unknown: raise ValueError(f"unknown source(s) {sorted(unknown)}")
        s = {k: require(sources.get(k, 0.0), k, 0.0) for k in SOURCES}; dbl = require(mantle_double_counted_in_ext, "mantle_double_counted_in_ext", 0.0)
        if dbl > s["S_ext"] + 1e-300: raise ValueError("declared double count exceeds S_ext")
        s["S_ext"] -= dbl; cap = require(capacity, "capacity", 0.0); S = sum(s.values())
        q = min(cap, S + self.Vs / dt)                 # actual export: cannot remove more than is available during the step
        self.Vs = self.Vs + dt * (S - q)
        tol = 1e-12 * max(1.0, S * dt, abs(self.Vs))
        if self.Vs < -tol: raise ArithmeticError(f"negative stock {self.Vs} after step (export logic violated)")
        if self.Vs < 0: self.cum["roundoff_correction"] += -self.Vs; self.Vs = 0.0      # recorded, not silent (included in budget_residual)
        for k in SOURCES: self.cum[k] += s[k] * dt
        self.cum["Q_out_actual"] += q * dt; self.cum["Q_out_capacity_unused"] += (cap - q) * dt; self.t += dt
        return dict(Q_out=q, S=S, Vs=self.Vs, removed_double_count=dbl)
    def budget_residual(self, Vs0):
        """independent reconstruction: Vs(t) - Vs0 - (sum of recorded sources - recorded export)."""
        return self.Vs - Vs0 - (sum(self.cum[k] for k in SOURCES) - self.cum["Q_out_actual"] + self.cum["roundoff_correction"])
    def state(self): return dict(Vs=self.Vs, t=self.t, cum=dict(self.cum))
