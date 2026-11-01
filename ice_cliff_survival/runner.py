"""Deterministic reduced toe-stock runner with restart. Structural baseline: closures are registry-backed definitions (closures.py) with explicit provenance; a run refuses to be labelled environmental.
Checkpoint identity (M2): version, config fingerprint, dt, label, closure TYPE+PARAMETERS (or extendable-forcing PREFIX CONTENT up to the checkpoint time), provenance. Opaque Python callables are not part of verified restarts:
they are accepted only with unsafe_opaque=True, can be saved/restored only with unsafe_restart=True, and such checkpoints are stamped restart_verified=False. Permitted change on restart: a longer horizon (forcing tables may be extended after the checkpoint time)."""
import json, math, os
from . import UnsupportedRun, VERSION
from .config import ToeConfig, Provenance, fingerprint
from .closures import Closure, from_spec
from .toe_stock import ToeStock
from .geometry import exposure
_IDENT_KEYS = ("version", "cfg", "dt", "label", "closures", "prov")
class ToeRun:
    def __init__(self, cfg: ToeConfig, sources, capacity, prov: dict, dt: float, label="synthetic", unsafe_opaque=False):
        if label not in ("synthetic", "component_test"): raise UnsupportedRun("environmental toe-stock runs are not supported in milestone 0 (no verified coupled flow/heat/ice/grain model); label must be 'synthetic' or 'component_test'")
        if isinstance(dt, bool) or not isinstance(dt, (int, float)) or not math.isfinite(dt) or dt <= 0.0: raise ValueError(f"dt must be a finite number > 0, got {dt!r}")
        for k in ("sources", "capacity"):
            if k not in prov or not isinstance(prov[k], Provenance): raise ValueError(f"closure '{k}' needs a Provenance")
            if prov[k].kind == "measured": raise UnsupportedRun(f"closure '{k}' claims measured provenance in a synthetic/component run: measured forcing may only enter a verified environmental model (not available)")
        self.opaque = [n for n, c in (("sources", sources), ("capacity", capacity)) if not isinstance(c, Closure)]
        if self.opaque and not unsafe_opaque: raise UnsupportedRun(f"opaque callable closure(s) {self.opaque} have no verifiable definition; use registry closures (closures.py) or pass unsafe_opaque=True (restart then needs unsafe_restart=True and is never verified)")
        self.cfg, self.sources, self.capacity, self.prov, self.dt, self.label = cfg, sources, capacity, prov, float(dt), label
        self.stock = ToeStock(Vs=cfg.Vs0); self.Vs0 = cfg.Vs0; self.history = []; self.restart_verified = not self.opaque
    def _closure_ident(self, t):
        out = {}
        for n, c in (("sources", self.sources), ("capacity", self.capacity)):
            out[n] = dict(opaque=True) if n in self.opaque else dict(type=c.spec()["type"], static=c.type_params_fingerprint(), prefix=c.prefix_fingerprint(t))
        return out
    def ident(self, t=None):
        """identity at checkpoint time t (default: current stock time)."""
        t = self.stock.t if t is None else t
        return dict(version=VERSION, cfg=self.cfg.fingerprint(), dt=self.dt, label=self.label, closures=self._closure_ident(t), prov={k: dict(v.__dict__, params=[list(p) for p in v.params]) for k, v in self.prov.items()})
    def advance(self, t_end):
        if isinstance(t_end, bool) or not isinstance(t_end, (int, float)) or not math.isfinite(t_end): raise ValueError(f"horizon must be finite, got {t_end!r}")
        if t_end < self.stock.t - 1e-12: raise ValueError(f"backward horizon {t_end} < current time {self.stock.t}")
        while self.stock.t < t_end - 1e-12:
            t = self.stock.t; dt = min(self.dt, t_end - t); src = self.sources(t + 0.5 * dt); cap = self.capacity(t + 0.5 * dt, self.stock.Vs); r = self.stock.step(dt, src, cap)
            ex = exposure(self.stock.Vs, self.cfg.W, self.cfg.porosity, self.cfg.H_cliff, self.cfg.theta_deg, self.cfg.beta_deg)
            self.history.append(dict(label=self.label, t=self.stock.t, Vs=self.stock.Vs, Q_out=r["Q_out"], S=r["S"], covered_fraction=ex["covered_fraction"], exposed_area=ex["exposed_area"], fully_buried=ex["fully_buried"]))
    def save(self, path, unsafe_restart=False):
        if self.opaque and not unsafe_restart: raise UnsupportedRun(f"cannot checkpoint: opaque closure(s) {self.opaque} have no verifiable identity (pass unsafe_restart=True for an UNVERIFIED checkpoint)")
        tmp = path + ".tmp"
        d = dict(ident=self.ident(), closure_specs={n: c.spec() for n, c in (("sources", self.sources), ("capacity", self.capacity)) if n not in self.opaque}, restart_verified=self.restart_verified, stock=self.stock.state(), Vs0=self.Vs0, history_len=len(self.history))
        with open(tmp, "w") as f: json.dump(d, f, indent=1)
        os.replace(tmp, path)
    def load(self, path, unsafe_restart=False):
        with open(path) as f: d = json.load(f)
        if not isinstance(d, dict) or not isinstance(d.get("ident"), dict) or not isinstance(d.get("stock"), dict): raise ValueError("malformed checkpoint: missing ident/stock")
        a = d["ident"]; miss = [k for k in _IDENT_KEYS if k not in a]
        if miss: raise ValueError(f"malformed checkpoint: identity fields missing {miss}")
        if not isinstance(a["closures"], dict) or set(a["closures"]) != {"sources", "capacity"}: raise ValueError("malformed checkpoint: closure identity")
        for k in ("Vs", "t", "cum"):
            if k not in d["stock"]: raise ValueError(f"malformed checkpoint: stock.{k} missing")
        if "Vs0" not in d or "restart_verified" not in d: raise ValueError("malformed checkpoint: Vs0/restart_verified missing")
        tchk = d["stock"]["t"]
        if isinstance(tchk, bool) or not isinstance(tchk, (int, float)) or not math.isfinite(tchk) or tchk < 0: raise ValueError("malformed checkpoint: bad time")
        if (self.opaque or not d["restart_verified"]) and not unsafe_restart: raise UnsupportedRun("checkpoint or run involves opaque closures: restart is unverifiable; pass unsafe_restart=True to accept an UNVERIFIED restart")
        try: b = self.ident(tchk)
        except ValueError as e: raise ValueError(f"restart refused: current definitions cannot cover the checkpoint time: {e}")
        bad = [k for k in _IDENT_KEYS if a[k] != b[k]]
        if bad: raise ValueError(f"restart refused: checkpoint differs in {bad}")
        if self.opaque: self.restart_verified = False
        self.stock.Vs = d["stock"]["Vs"]; self.stock.t = tchk; self.stock.cum = dict(d["stock"]["cum"]); self.Vs0 = d["Vs0"]
