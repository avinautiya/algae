"""Deterministic reduced toe-stock runner with restart. Structural baseline: closures are injected (sources(t) -> dict, capacity(t, Vs) -> m3/s) with explicit provenance; a run refuses to be labelled environmental.
State file: JSON with config fingerprint, forcing fingerprint, closure provenance, stock state, ledger totals; restart is refused on any fingerprint mismatch (permitted change: longer horizon)."""
import json, os
from . import UnsupportedRun, VERSION
from .config import ToeConfig, Provenance, fingerprint
from .toe_stock import ToeStock
from .geometry import exposure
class ToeRun:
    def __init__(self, cfg: ToeConfig, sources, capacity, prov: dict, forcing_id: str, dt: float, label="synthetic"):
        if label not in ("synthetic", "component_test"): raise UnsupportedRun("environmental toe-stock runs are not supported in milestone 0 (no verified coupled flow/heat/ice/grain model); label must be 'synthetic' or 'component_test'")
        for k in ("sources", "capacity"):
            if k not in prov or not isinstance(prov[k], Provenance): raise ValueError(f"closure '{k}' needs a Provenance")
        self.cfg, self.sources, self.capacity, self.prov, self.forcing_id, self.dt, self.label = cfg, sources, capacity, prov, forcing_id, dt, label
        self.stock = ToeStock(Vs=cfg.Vs0); self.Vs0 = cfg.Vs0; self.history = []
    def ident(self): return dict(version=VERSION, cfg=self.cfg.fingerprint(), forcing=self.forcing_id, dt=self.dt, prov={k: v.__dict__ for k, v in self.prov.items()}, label=self.label)
    def advance(self, t_end):
        while self.stock.t < t_end - 1e-12:
            t = self.stock.t; dt = min(self.dt, t_end - t); src = self.sources(t + 0.5 * dt); cap = self.capacity(t + 0.5 * dt, self.stock.Vs); r = self.stock.step(dt, src, cap)
            ex = exposure(self.stock.Vs, self.cfg.W, self.cfg.porosity, self.cfg.H_cliff, self.cfg.theta_deg, self.cfg.beta_deg)
            self.history.append(dict(t=self.stock.t, Vs=self.stock.Vs, Q_out=r["Q_out"], S=r["S"], covered_fraction=ex["covered_fraction"], exposed_area=ex["exposed_area"], fully_buried=ex["fully_buried"]))
    def save(self, path):
        tmp = path + ".tmp"; json.dump(dict(ident=self.ident(), stock=self.stock.state(), Vs0=self.Vs0, history_len=len(self.history)), open(tmp, "w"), indent=1); os.replace(tmp, path)
    def load(self, path):
        d = json.load(open(path)); a, b = d["ident"], self.ident()
        bad = [k for k in a if a[k] != b[k]]
        if bad: raise ValueError(f"restart refused: checkpoint differs in {bad}")
        self.stock.Vs = d["stock"]["Vs"]; self.stock.t = d["stock"]["t"]; self.stock.cum = dict(d["stock"]["cum"]); self.Vs0 = d["Vs0"]
