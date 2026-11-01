"""Exposure geometry of a planar cliff toe with a debris wedge (2-D cross-section per unit width, then x width W).
Cliff face rises at angle theta from the horizontal ground; the debris pile has free-surface repose angle beta < theta and leans against the face. The pile cross-section is a triangle with
  area A = 1/2 z_f^2 (cot beta - cot theta)   (z_f = pile height on the face; bulk area per unit width; solid area = (1 - porosity) A).
So z_f = sqrt(2 A / (cot beta - cot theta)). The face is geometrically buried when z_f >= H_cliff. This is a GEOMETRIC coverage statement only: water contact, export and rebuilding decide survival (see MODEL_FOUNDATIONS.md)."""
import math
from .units import require
def wedge_height(A_bulk, theta_deg, beta_deg):
    th, be = math.radians(require(theta_deg, "theta", 0.0, 90.0, lo_strict=True)), math.radians(require(beta_deg, "beta", 0.0, 90.0, lo_strict=True))
    if be >= th: raise ValueError("repose angle must be smaller than the cliff face angle for a wedge against the face")
    return math.sqrt(2.0 * require(A_bulk, "A_bulk", 0.0) / (1.0 / math.tan(be) - 1.0 / math.tan(th)))
def wedge_area(z_f, theta_deg, beta_deg):
    th, be = math.radians(theta_deg), math.radians(beta_deg); return 0.5 * z_f ** 2 * (1.0 / math.tan(be) - 1.0 / math.tan(th))
def exposure(Vs_solid, W, porosity, H_cliff, theta_deg, beta_deg, water_depth_at_toe=0.0):
    """returns exposed/covered/wetted face quantities. Vs_solid [m3 solid] in the toe wedge of width W [m]."""
    A_bulk = Vs_solid / (1.0 - porosity) / W; z_f = wedge_height(A_bulk, theta_deg, beta_deg); sin_t = math.sin(math.radians(theta_deg))
    z_cov = min(z_f, H_cliff); face_len = H_cliff / sin_t; cov_len = z_cov / sin_t; exp_len = face_len - cov_len
    wet_exposed = max(0.0, min(water_depth_at_toe, H_cliff) - z_cov) / sin_t         # water standing ABOVE the pile on the exposed face (structural assumption: water in the pile pores is not ice contact)
    return dict(pile_height_on_face=z_f, covered_fraction=cov_len / face_len, exposed_area=W * exp_len, covered_area=W * cov_len, wetted_exposed_area=W * wet_exposed, fully_buried=z_f >= H_cliff)
def critical_solid_volume(W, porosity, H_cliff, theta_deg, beta_deg):
    """solid volume at which the face is geometrically fully covered: (1 - phi) W 1/2 H^2 (cot beta - cot theta). A geometric covering volume, not a universal burial threshold."""
    return (1.0 - porosity) * W * wedge_area(H_cliff, theta_deg, beta_deg)
def _check_angles(theta_deg, beta_deg):
    for n, v in (("theta", theta_deg), ("beta", beta_deg)):
        if not (isinstance(v, (int, float)) and math.isfinite(v)): raise ValueError(f"{n} must be finite")
    if not (0.0 < beta_deg < theta_deg <= 90.0): raise ValueError("require 0 < beta < theta <= 90 degrees")
def debris_normal_thickness(z, z_f, theta_deg, beta_deg):
    """normal-to-face thickness of the static debris wedge at height z on the face (M2 CORRECTED).
    Geometry (open side to the right): face line through (-z cot th, z) (back-tilted, ice behind), ground z=0, unit normal into the open side n=(sin th, cos th); pile surface descends from the face point at height z_f at angle beta. The normal ray from the face point at height z hits the pile surface after
        d = (z_f - z) tan(theta - beta) / sin(theta);  zero for z >= z_f.  (theta = 90 deg: d = (z_f - z) tan(90deg - beta) = (z_f-z) cot(beta).)
    Requires 0 < beta < theta <= 90 deg. The M1 formula is kept as debris_normal_thickness_v1_projected_horizontal (provenance only)."""
    _check_angles(theta_deg, beta_deg)
    if not (math.isfinite(z) and math.isfinite(z_f) and z_f >= 0.0): raise ValueError("z, z_f must be finite, z_f >= 0")
    if z >= z_f: return 0.0
    th, be = math.radians(theta_deg), math.radians(beta_deg)
    return (z_f - z) * math.tan(th - be) / math.sin(th)
def debris_normal_thickness_v1_projected_horizontal(z, z_f, theta_deg, beta_deg):
    """DEPRECATED M1 (a680fae) formula, wrong: horizontal width times sin(theta). Kept so old results are reproducible; do not use."""
    th, be = math.radians(theta_deg), math.radians(beta_deg)
    return 0.0 if z >= z_f else (z_f - z) * (1.0 / math.tan(be) - 1.0 / math.tan(th)) * math.sin(th)
def shield_conductance(z, z_f, theta_deg, beta_deg, k_deb, d_min=1e-4):
    """(M2: uses the corrected normal thickness.) static-shield interface conductance h_eff = k_deb / d [W m-2 K-1]; k/d is a 1-D effective-resistance (massless, normal-conduction) approximation, not a resolved wedge-conduction result. for the 2-D heat/moving-ice solver (hShield list; 1e30 = bare). d floored at d_min to avoid unbounded h at the pile crest. Debris heat capacity and pore water are ignored (quasi-steady dry layer; k_deb is an assumed input)."""
    d = debris_normal_thickness(z, z_f, theta_deg, beta_deg)
    return 1e30 if d <= 0 else k_deb / max(d, d_min)
