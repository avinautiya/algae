"""Field-forcing interface: elapsed-time temperature integration with SEPARATE observed / reconstructed / unsupported accounting.
Rules enforced here:
 * calibration offsets (P3 +0.035 C, P5 +0.040 C) are applied exactly once, inside the Codex source adapter; `extra_offset_C` is an explicit labelled sensitivity, never a correction;
 * integrals are over ELAPSED TIME of adjacent finite supported samples (irregular sampling handled by the piecewise-linear primitive);
 * strict native support (adjacent samples <=60 s) is reported alongside the <=24 h reconstructed version; reconstruction is never silent;
 * a missing/purged record yields UNSUPPORTED time; it is never converted to zero temperature, zero flow or lost contact.
"""
import hashlib, json, os
from dataclasses import dataclass, asdict
from . import _paths
from gap_reconstructed_temperature import GapReconstructedTemperatureSeries
from original_temperature_adapter import OFFSETS
@dataclass(frozen=True)
class SupportAccount:
    requested_s: float; observed_s: float; reconstructed_s: float; unsupported_s: float
    observed_integral_C_s: float; reconstructed_integral_C_s: float; method: str; extra_offset_C: float
    @property
    def supported_s(self): return self.observed_s + self.reconstructed_s
    @property
    def strict_native_only(self): return self.unsupported_s + self.reconstructed_s == 0.0
    def mean_temperature(self, include_reconstructed=True):
        """time-weighted mean over SUPPORTED time (None if none); unsupported time is not weighted as zero"""
        s = self.supported_s if include_reconstructed else self.observed_s; i = self.observed_integral_C_s + (self.reconstructed_integral_C_s if include_reconstructed else 0.0)
        return None if s <= 0 else i / s
def _sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()
class FieldForcing:
    """Wraps the Codex gap-reconstructed series for one station. `method`: linear_endpoint_bridge (primary, selected on a temperature-only mask experiment), constant_previous_endpoint, constant_next_endpoint."""
    def __init__(self, station, source_root=_paths.SWISS_SOURCE, method="linear_endpoint_bridge", max_reconstruction_gap_s=86400):
        self.station = int(str(station).removeprefix("P")); self.method = method; self.max_gap = max_reconstruction_gap_s; self.source_root = source_root
        self.series = GapReconstructedTemperatureSeries.from_archive(source_root, self.station, max_observed_gap_seconds=60, max_reconstruction_gap_seconds=max_reconstruction_gap_s)
        self.offset_applied_C = OFFSETS[self.station]
    def integrate(self, start, end, extra_offset_C=0.0):
        r = self.series.integrate_positive(start, end, extra_offset_C=extra_offset_C, method=self.method)
        return SupportAccount(requested_s=r["requested_seconds"], observed_s=r["observed_seconds"], reconstructed_s=r["reconstructed_seconds"], unsupported_s=r["requested_seconds"] - r["covered_seconds"],
                              observed_integral_C_s=r["observed_integral_C_s"], reconstructed_integral_C_s=r["reconstructed_integral_C_s"], method=self.method, extra_offset_C=extra_offset_C)
    def provenance(self):
        man = os.path.join(self.source_root, "SOURCE_MANIFEST.json")
        return dict(station=f"P{self.station}", source_manifest_sha256=_sha(man), calibration_offset_applied_once_C=self.offset_applied_C, method=self.method, max_observed_gap_s=60, max_reconstruction_gap_s=self.max_gap, clock="naive/unresolved (no timezone inferred)")
def account_from_saved_support(d, method="linear_endpoint_bridge"):
    """convert the support dict saved by the Codex V4/V5 runs (observed/reconstructed seconds and integrals) into a SupportAccount"""
    return SupportAccount(requested_s=d["requested_seconds"], observed_s=d["observed_seconds"], reconstructed_s=d["reconstructed_seconds"], unsupported_s=d["requested_seconds"] - d["covered_seconds"],
                          observed_integral_C_s=d["observed_integral_C_s"], reconstructed_integral_C_s=d["reconstructed_integral_C_s"], method=method, extra_offset_C=d.get("extra_offset_C", 0.0))
