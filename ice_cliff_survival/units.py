"""SI constants and validation helpers. All quantities SI unless a name says otherwise (…_cm, …_deg)."""
import math
RHO_W = 1000.0; RHO_I = 917.0; CP_W = 4217.0; LF = 3.34e5; K_W = 0.561; G = 9.81; RHO_ROCK_DEFAULT = 2700.0
def require(x, name, lo=None, hi=None, lo_strict=False, hi_strict=False):
    if not isinstance(x, (int, float)) or isinstance(x, bool) or not math.isfinite(x): raise ValueError(f"{name} must be a finite number, got {x!r}")
    if lo is not None and (x < lo or (lo_strict and x == lo)): raise ValueError(f"{name}={x} below lower bound {lo}")
    if hi is not None and (x > hi or (hi_strict and x == hi)): raise ValueError(f"{name}={x} above upper bound {hi}")
    return float(x)
