"""Typed configuration with validation and deterministic fingerprints (canonical JSON + sha256)."""
import hashlib, json
from dataclasses import dataclass, asdict, field
from .units import require
def fingerprint(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(",", ":"), default=float).encode()).hexdigest()
@dataclass(frozen=True)
class ToeConfig:
    """Reduced toe-stock configuration (structural baseline). W [m] cliff-toe width; porosity of deposited debris; H_cliff [m]; theta_deg cliff face angle; beta_deg debris repose angle (< theta); Vs0 [m3 solid]."""
    W: float; porosity: float; H_cliff: float; theta_deg: float; beta_deg: float; Vs0: float = 0.0
    def __post_init__(self):
        require(self.W, "W", 0.0, lo_strict=True); require(self.porosity, "porosity", 0.0, 1.0, hi_strict=True); require(self.H_cliff, "H_cliff", 0.0, lo_strict=True)
        require(self.theta_deg, "theta_deg", 0.0, 90.0, lo_strict=True); require(self.beta_deg, "beta_deg", 0.0, 90.0, lo_strict=True); require(self.Vs0, "Vs0", 0.0)
        if self.beta_deg >= self.theta_deg: raise ValueError("beta_deg (repose) must be < theta_deg")
    def fingerprint(self): return fingerprint(asdict(self))
@dataclass(frozen=True)
class Provenance:
    """Every injectable closure declares where its numbers come from. kind in {'measured','published','structural','assumed','synthetic'}."""
    name: str; kind: str; source: str = ""
    def __post_init__(self):
        if self.kind not in ("measured", "published", "structural", "assumed", "synthetic"): raise ValueError("bad provenance kind")
