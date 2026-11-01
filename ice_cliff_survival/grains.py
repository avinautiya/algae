"""Grain bookkeeping on SOLID volume (pore space excluded).
Conventions: n points from ice into the control volume (CV)/water; w = velocity of the CV boundary; flux of grains across a boundary element = c (u_s - w).n dA with the grain solid-volume concentration c (solid volume per mixture volume),
positive INTO the CV when (u_s - w).n > 0 for an inward-pointing n (the caller supplies the inward normal). Frozen englacial grains have u_s = u_ice.
Englacial concentration: c_mix = solid volume per MIXTURE (ice + grains) volume; c_ice = solid volume per PURE-ICE volume; c_mix = c_ice / (1 + c_ice)."""
from .units import require
def c_mix_from_c_ice(c_ice): c = require(c_ice, "c_ice", 0.0); return c / (1.0 + c)
def c_ice_from_c_mix(c_mix): c = require(c_mix, "c_mix", 0.0, 1.0, hi_strict=True); return c / (1.0 - c)
def solid_from_bulk(V_bulk, porosity): return require(V_bulk, "V_bulk", 0.0) * (1.0 - require(porosity, "porosity", 0.0, 1.0, hi_strict=True))
def bulk_from_solid(V_solid, porosity): return require(V_solid, "V_solid", 0.0) / (1.0 - require(porosity, "porosity", 0.0, 1.0, hi_strict=True))
def relative_flux(c, u_s_n, w_n, area):
    """grain solid-volume flux [m3/s] INTO the CV through a boundary element: c (u_s - w).n_in * area (u_s_n, w_n are components along the INWARD normal)."""
    return require(c, "c", 0.0, 1.0) * (u_s_n - w_n) * require(area, "area", 0.0)
def mantle_source(porosity, W, d_top, v_top, delivery_fraction):
    """S_mantle = (1 - porosity) W d_top v_top f_delivery  [m3 solid / s]; v_top = retreat speed of the UPPER lip (m/s), NOT the toe incision rate."""
    return (1.0 - require(porosity, "porosity", 0.0, 1.0, hi_strict=True)) * require(W, "W", 0.0) * require(d_top, "d_top", 0.0) * require(v_top, "v_top", 0.0) * require(delivery_fraction, "delivery_fraction", 0.0, 1.0)
def englacial_source(c_mix, face_area, retreat_speed_normal):
    """S_eng = c_mix A r_n [m3 solid / s]: the CV boundary (cliff face) moves into ice at geometric speed r_n and sweeps mixture volume A r_n; equivalent to relative_flux with u_s = u_ice = 0 and w.n_in = -r_n."""
    return require(c_mix, "c_mix", 0.0, 1.0) * require(face_area, "face_area", 0.0) * require(retreat_speed_normal, "retreat_speed_normal", 0.0)
def ice_mass_melt_rate(c_mix, face_area, retreat_speed_normal, rho_ice):
    """kg/s of pure ice melted when the mixture recedes at r_n: rho_ice (1 - c_mix) A r_n  (grains are not melted)."""
    return require(rho_ice, "rho_ice", 0.0, lo_strict=True) * (1.0 - require(c_mix, "c_mix", 0.0, 1.0)) * face_area * retreat_speed_normal
def mantle_split(porosity, W, d_top, v_top, delivery_fraction):
    """(delivered, undelivered) solid-volume rates [m3/s]; delivered + undelivered = total released solid (1-phi) W d v. Undelivered material is a recorded sink of the mantle stock, NOT available as external rockfall."""
    tot = mantle_source(porosity, W, d_top, v_top, 1.0); d = mantle_source(porosity, W, d_top, v_top, delivery_fraction); return d, tot - d
