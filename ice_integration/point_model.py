"""Point-control (station floor) predictions of the frozen V5 families expressed through the integrated interface. Heat transfer h = k Nu / D_h per hydraulic cell, melt-only positive temperature integral over SUPPORTED time (observed + labelled reconstruction), cost = rho_i L (temperate ice).
This reproduces the Codex V5 predictor exactly; it adds no model physics and no finite-heat depletion (a station sensor measures the local floor temperature; there is no supplied inlet-to-station path)."""
import math
K, RHO, LATENT = 0.56, 917.0, 334000.0
FAMILIES = ("constant_rate_fitted", "constant_h_fitted", "db04_unfitted", "db03_unfitted", "gnielinski_smooth_unfitted", "db04_scale_fitted")
NU_KEY = {"db04_unfitted": "db04", "db03_unfitted": "db03", "gnielinski_smooth_unfitted": "gn"}
def basis(cells, key):
    return math.fsum(K / c["hydraulic_diameter_m"] / (RHO * LATENT) * c["Nu"][key] * c["integral_C_s"] for c in cells)
def predict(interval, calibration):
    cells = interval["hydraulic_cells"]; b = {n: basis(cells, k) for n, k in NU_KEY.items()}; tint = math.fsum(c["integral_C_s"] for c in cells)
    out = {"constant_rate_fitted": calibration["constant_rate_m_s"] * interval["duration_s"], "constant_h_fitted": calibration["constant_h_W_m2_K"] * tint / (RHO * LATENT), "db04_scale_fitted": calibration["db04_scale"] * b["db04_unfitted"]}
    out.update(b); return out
def metrics(errors):
    n = len(errors); return dict(n=n, MAE_m=math.fsum(abs(e) for e in errors) / n, RMSE_m=math.sqrt(math.fsum(e * e for e in errors) / n), bias_m=math.fsum(errors) / n)
