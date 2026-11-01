"""Explicit interface between the field-forcing adapter and the coupled plan-view model (reactivation/swe2d.Sim).
For each forcing interval: supported -> set Sim.Tw (time-weighted mean over supported time), Sim.Q, march hydraulics, then ONE bed update of length = supported seconds;
unsupported -> NO bed update, NO hydraulic inference (contact/flow are not set to zero); durations kept separate. Restart files carry the bed AND the flow state (so restart is exactly equivalent), counters and source fingerprints."""
import hashlib, json, os
import numpy as np
from .forcing import SupportAccount
def _sha(b): return hashlib.sha256(b).hexdigest()
class CoupledDriver:
    def __init__(self, sim, provenance, include_reconstructed=True, hyd_max=60.0):
        self.sim = sim; self.provenance = dict(provenance); self.include_reconstructed = include_reconstructed; self.hyd_max = hyd_max
        self.time = dict(observed_s=0.0, reconstructed_s=0.0, unsupported_s=0.0, advanced_s=0.0); self.log = []; self.melt = dict(toe=0.0, byp=0.0)
    def advance(self, account: SupportAccount, q_m3_s):
        a = account; self.time["observed_s"] += a.observed_s; self.time["unsupported_s"] += a.unsupported_s
        if self.include_reconstructed: self.time["reconstructed_s"] += a.reconstructed_s
        else: self.time["unsupported_s"] += a.reconstructed_s
        dur = a.supported_s if self.include_reconstructed else a.observed_s; tm = a.mean_temperature(self.include_reconstructed)
        if dur <= 0 or tm is None: self.log.append(dict(advanced=False, reason="unsupported interval: bed/contact state unchanged and not inferred")); return None
        s = self.sim; s.Tw = float(tm); s.Q = float(q_m3_s); s.hydraulics(tmax=self.hyd_max, cool=True); d = s.diagnostics(); vol, vt, vb = s.lower_bed(dur)
        self.time["advanced_s"] += dur; self.melt["toe"] += vt; self.melt["byp"] += vb; self.log.append(dict(advanced=True, seconds=dur, Tw=s.Tw, Q=s.Q, f_toe=d["f_toe"], ledger_gap=s.ledger_gap)); return d
    def checkpoint(self, path):
        meta = dict(provenance=self.provenance, time=self.time, melt=self.melt, include_reconstructed=self.include_reconstructed, n_log=len(self.log), dx=self.sim.dx, solver="reactivation.swe2d v4")
        tmp = path + ".tmp.npz"; np.savez(tmp, zb=self.sim.zb, h=self.sim.h, hu=self.sim.hu, hv=self.sim.hv, hT=self.sim.hT, meta=json.dumps(meta, sort_keys=True, default=float)); os.replace(tmp, path); return meta
    def restore(self, path):
        z = np.load(path, allow_pickle=False); meta = json.loads(str(z["meta"]))
        if meta["provenance"] != self.provenance or meta["dx"] != self.sim.dx or meta["include_reconstructed"] != self.include_reconstructed: raise ValueError("restart refused: provenance/choices differ")
        self.sim.zb = z["zb"].copy(); self.sim.h = z["h"].copy(); self.sim.hu = z["hu"].copy(); self.sim.hv = z["hv"].copy(); self.sim.hT = z["hT"].copy(); self.time = meta["time"]; self.melt = meta["melt"]; return meta
