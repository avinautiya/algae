"""Serializable, registry-backed closure definitions for ToeRun (M2).
A closure is an executable model definition that can state EXACTLY what it computes: `spec()` returns a JSON-able {"type", "params"} from which `from_spec()` rebuilds an equivalent closure through the registry.
Restart identity is derived from these specs (type + numeric parameters + actual forcing-table contents), never from function names, repr() or descriptive strings.
Forcing tables are piecewise constant on declared edges; `prefix_fingerprint(t)` hashes the forcing content that determines the closure on [0, t], so a horizon extension (adding intervals after t) is a permitted change and anything altering [0, t] is not.
Extension policy is declared per table: extension='forbid' (evaluation past the last edge raises) is the only policy implemented."""
import bisect, math
from .config import fingerprint
from .units import require
REGISTRY = {}
def register(cls): REGISTRY[cls.TYPE] = cls; return cls
class Closure:
    TYPE = None
    def spec(self): raise NotImplementedError
    def prefix_fingerprint(self, t): return fingerprint(self.spec())     # default: whole definition is the identity (no extendable content)
    def type_params_fingerprint(self): s = self.spec(); return fingerprint(dict(type=s["type"], params=self._static_params()))
    def _static_params(self): return self.spec()["params"]
@register
class ConstantSources(Closure):
    TYPE = "constant_sources"
    def __init__(self, rates): self.rates = {k: require(v, k, 0.0) for k, v in dict(rates).items()}
    def __call__(self, t): return dict(self.rates)
    def spec(self): return dict(type=self.TYPE, params=dict(rates=dict(sorted(self.rates.items()))))
    @classmethod
    def from_spec(cls, s): return cls(s["params"]["rates"])
@register
class ConstantCapacity(Closure):
    TYPE = "constant_capacity"
    def __init__(self, Q): self.Q = require(Q, "Q", 0.0)
    def __call__(self, t, Vs): return self.Q
    def spec(self): return dict(type=self.TYPE, params=dict(Q=self.Q))
    @classmethod
    def from_spec(cls, s): return cls(s["params"]["Q"])
@register
class LinearCapacity(Closure):
    """Q = max(0, c0 + c1 * Vs) [m3/s]."""
    TYPE = "linear_capacity"
    def __init__(self, c0, c1): self.c0, self.c1 = require(c0, "c0"), require(c1, "c1")
    def __call__(self, t, Vs): return max(0.0, self.c0 + self.c1 * Vs)
    def spec(self): return dict(type=self.TYPE, params=dict(c0=self.c0, c1=self.c1))
    @classmethod
    def from_spec(cls, s): return cls(**s["params"])
@register
class TableSources(Closure):
    """piecewise-constant source rates: edges t_0 < ... < t_n (t_0 = 0), rates[i] = dict valid on [t_i, t_{i+1}). extension='forbid'."""
    TYPE = "table_sources"
    def __init__(self, edges, rates, extension="forbid"):
        if extension != "forbid": raise ValueError("only extension='forbid' is implemented")
        e = [require(x, "edge") for x in edges]
        if len(e) < 2 or e[0] != 0.0 or any(b <= a for a, b in zip(e, e[1:])): raise ValueError("edges must start at 0 and increase strictly")
        if len(rates) != len(e) - 1: raise ValueError("need one rate dict per interval")
        self.edges, self.rates, self.extension = e, [{k: require(v, k, 0.0) for k, v in dict(r).items()} for r in rates], extension
    def __call__(self, t):
        if not (0.0 <= t <= self.edges[-1] * (1 + 1e-12)): raise ValueError(f"forcing table covers [0, {self.edges[-1]}], asked t={t}")
        i = min(max(bisect.bisect_right(self.edges, t) - 1, 0), len(self.rates) - 1); return dict(self.rates[i])
    def spec(self): return dict(type=self.TYPE, params=dict(edges=self.edges, rates=[dict(sorted(r.items())) for r in self.rates], extension=self.extension))
    def _static_params(self): return dict(extension=self.extension)
    def prefix_fingerprint(self, t):
        n = sum(1 for a in self.edges[:-1] if a < t) if t > 0 else 0
        if t > self.edges[-1]: raise ValueError("checkpoint time beyond table end")
        return fingerprint(dict(type=self.TYPE, intervals=[[self.edges[i], dict(sorted(self.rates[i].items()))] for i in range(n)]))
    @classmethod
    def from_spec(cls, s): p = s["params"]; return cls(p["edges"], p["rates"], p.get("extension", "forbid"))
def from_spec(s):
    if not isinstance(s, dict) or s.get("type") not in REGISTRY or "params" not in s: raise ValueError(f"unknown or malformed closure spec: {s!r}")
    return REGISTRY[s["type"]].from_spec(s)
