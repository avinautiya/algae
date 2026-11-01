"""Finite incoming-water-heat accounting for supplied routes (wraps Codex predict_external.predict_interval; review notes in INTEGRATION_REPORT.md).
Each SUPPORTED interval uses the time-weighted mean positive temperature over supported time: because h depends on Re (hence on Q), not on T, delivered power is LINEAR in T and the mean is exact for the exponential-depletion law.
Unsupported time is neither zero-filled nor extrapolated: it is accumulated separately and no heat is allocated to it.
Heat is allocated ONCE along the ordered route (section outlet temperature becomes the next inlet); within a section bed/sidewall split by wetted-perimeter fraction.
Separate ledgers: floor, wall, outgoing water, latent, fresh-ice sensible; atmospheric and dissipation are NOT implemented (reported as such); no moving geometry, no cold-ice PDE.
"""
import hashlib, json, math
from . import _paths
import predict_external as PE
from .forcing import SupportAccount
LEDGERS = ("incoming_water_heat_J", "floor_heat_J", "wall_heat_J", "outgoing_water_heat_J", "latent_J", "fresh_ice_sensible_J", "floor_ice_volume_m3", "wall_ice_volume_m3")
class IntegratedHeatModel:
    def __init__(self, sections, method, ice_temperature_C=0.0, metadata=None):
        if method not in PE.PHYSICAL_MODELS: raise PE.ApplicabilityRefusal("physical model required")
        self.sections = list(sections); self.method = method; self.ice_temperature_C = ice_temperature_C; self.metadata = dict(metadata or {})
        self.ledger = {k: 0.0 for k in LEDGERS}; self.time = dict(requested_s=0.0, observed_s=0.0, reconstructed_s=0.0, unsupported_s=0.0); self.history = []
    def config_fingerprint(self):
        cfg = dict(method=self.method, ice_T=self.ice_temperature_C, sections=[(s.name, s.width_m, s.depth_m, s.length_m, s.point_control, s.source, s.geometry_basis) for s in self.sections], metadata=self.metadata)
        return hashlib.sha256(json.dumps(cfg, sort_keys=True, default=str).encode()).hexdigest()
    def advance(self, account: SupportAccount, q_m3_s, include_reconstructed=True):
        a = account; self.time["requested_s"] += a.requested_s; self.time["observed_s"] += a.observed_s; self.time["unsupported_s"] += a.unsupported_s
        use_rec = include_reconstructed
        if use_rec: self.time["reconstructed_s"] += a.reconstructed_s
        tmean = a.mean_temperature(include_reconstructed=use_rec); dur = a.supported_s if use_rec else a.observed_s
        if not use_rec: self.time["unsupported_s"] += a.reconstructed_s          # strict mode: reconstructed time is unsupported, not zero heat
        if tmean is None or dur <= 0 or q_m3_s is None:
            self.history.append(dict(supported_s=0.0, note="no supported time: no heat allocated, contact/flow NOT inferred")); return None
        r = PE.predict_interval(self.sections, q_m3_s, tmean, dur, self.ice_temperature_C, self.method)
        floor = sum(s["bed_heat_J"] for s in r["sections"]); wall = sum(s["sidewall_heat_J"] for s in r["sections"])
        if abs(floor + wall - r["removed_water_heat_J"]) > 1e-9 * max(1.0, r["removed_water_heat_J"]): raise AssertionError("heat allocated more than once")
        self.ledger["incoming_water_heat_J"] += r["inlet_heat_J"]; self.ledger["floor_heat_J"] += floor; self.ledger["wall_heat_J"] += wall; self.ledger["outgoing_water_heat_J"] += r["outlet_heat_J"]
        self.ledger["latent_J"] += r["latent_heat_J"]; self.ledger["fresh_ice_sensible_J"] += r["fresh_ice_sensible_heat_J"]
        self.ledger["floor_ice_volume_m3"] += sum(s["bed_ice_volume_m3"] for s in r["sections"]); self.ledger["wall_ice_volume_m3"] += sum(s["sidewall_ice_volume_m3"] for s in r["sections"])
        self.history.append(dict(supported_s=dur, mean_T_C=tmean, Q=q_m3_s, energy_residual_J=r["energy_residual_J"], interface_residual_J=r["interface_energy_residual_J"])); return r
    def residuals(self):
        L = self.ledger
        return dict(water_side_J=L["incoming_water_heat_J"] - L["outgoing_water_heat_J"] - L["floor_heat_J"] - L["wall_heat_J"], interface_J=L["floor_heat_J"] + L["wall_heat_J"] - L["latent_J"] - L["fresh_ice_sensible_J"],
                    implemented=["finite water heat depletion", "floor/wall split", "latent", "fresh-ice sensible (algebraic)"], not_implemented=["atmospheric exchange", "frictional dissipation", "moving geometry", "cold-ice conduction PDE", "meltwater feedback on Q or heat capacity"])
    def state(self): return dict(config_fingerprint=self.config_fingerprint(), ledger=dict(self.ledger), time=dict(self.time), n_history=len(self.history), metadata=self.metadata)
    @classmethod
    def restore(cls, sections, method, state, ice_temperature_C=0.0, metadata=None):
        m = cls(sections, method, ice_temperature_C, metadata)
        if m.config_fingerprint() != state["config_fingerprint"]: raise ValueError("restart refused: configuration fingerprint differs")
        m.ledger = dict(state["ledger"]); m.time = dict(state["time"]); return m
