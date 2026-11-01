"""Grain removal at the cliff toe: why a deep-bed transport law is not used, and what is implemented instead.
A conventional sediment-transport law (Meyer-Peter-Muller, Shields-based bedload, Hjulstrom) presumes a bed of the transported material much deeper than the grain and a flow depth >> grain size, with the bed shear stress acting on a mobile layer.
At a stream-fed cliff toe the bed is ICE (smooth, friction coefficient mu ~ 0.05, Petersen et al. 2024 Sect. 5.2), clasts are one layer thick or isolated, and clast diameter is comparable to or larger than the flow depth (measured depths 12-35 cm, velocities 0.6-1.13 m/s, Q 0.06-0.21 m3/s). The implemented treatment is therefore a
CLAST-SCALE FORCE BALANCE on a smooth ice bed (mobilisation threshold with partial submergence), plus an explicit, injectable clearing model that converts 'mobile' into a solid-volume export capacity. Both parts carry provenance flags; the clearing parameters are ASSUMPTIONS to be supplied, not calibrated.
Provenance: threshold = 'structural' (force balance, validated only against the published Petersen numbers: tests/test_transport.py); clearing = 'assumed' (no data)."""
import math
from .units import RHO_W, RHO_ROCK_DEFAULT, G, require
SHAPES = {"half_sphere": dict(Cw=0.42), "cube": dict(Cw=1.05)}
def critical_velocity_submerged(d, shape="half_sphere", mu=0.05, rho_s=RHO_ROCK_DEFAULT, rho_w=RHO_W, buoyancy=False, g=G):
    """fully submerged clast of characteristic diameter/edge d resting on ice: drag 1/2 rho_w Cw A v^2 = mu F_N.
    half_sphere (diameter d, flat face down): frontal area = semicircle pi d^2/8, V = pi d^3/12 ; cube (edge d): A = d^2, V = d^3. F_N = (rho_s - rho_w*buoyancy) V g  (Petersen et al. neglect buoyancy: buoyancy=False reproduces their Eqs. 2-3)."""
    d = require(d, "d", 0.0, lo_strict=True); Cw = SHAPES[shape]["Cw"]; A, V = (math.pi * d * d / 8, math.pi * d ** 3 / 12) if shape == "half_sphere" else (d * d, d ** 3)
    FN = (rho_s - (rho_w if buoyancy else 0.0)) * V * g; return math.sqrt(2.0 * mu * FN / (Cw * rho_w * A))
def max_mobile_diameter(v, shape="half_sphere", mu=0.05, rho_s=RHO_ROCK_DEFAULT, rho_w=RHO_W, buoyancy=False, g=G):
    """largest d with v >= v_c(d): v_c^2 is linear in d, so d_max = v^2 / (v_c(1)^2)."""
    v1 = critical_velocity_submerged(1.0, shape, mu, rho_s, rho_w, buoyancy, g); return v * v / (v1 * v1)
def partial_submergence_factor(d, h):
    """fraction of a sphere's vertical extent that is wetted when depth h < d (>=1 if h >= d): crude interface for the regime d ~ h; the force balance is NOT validated for h < d (flagged). Returns min(1, h/d)."""
    return min(1.0, require(h, "h", 0.0) / require(d, "d", 0.0, lo_strict=True))
class ClearingModel:
    """capacity [m3 solid/s] = mobile_fraction(v) * Vs * (zeta v / L_clear): the stock's mobile part leaves the toe zone of length L_clear at the clast speed zeta v. zeta, L_clear are ASSUMED inputs (no default)."""
    def __init__(self, d_clast, v_critical, zeta, L_clear, provenance="assumed"):
        self.d, self.vc, self.zeta, self.L = require(d_clast, "d", 0.0, lo_strict=True), require(v_critical, "v_c", 0.0, lo_strict=True), require(zeta, "zeta", 0.0, 1.0, lo_strict=True), require(L_clear, "L_clear", 0.0, lo_strict=True); self.provenance = provenance
    def capacity(self, Vs, v):
        if v < self.vc: return 0.0
        return Vs * self.zeta * v / self.L
