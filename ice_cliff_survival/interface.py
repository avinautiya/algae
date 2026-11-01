"""Ice/water interface kinematics (n points from ice into water). m_n = (u_ice - v_I).n >= 0 is the ice mass flux per unit ice density = recession speed of the ice along -n.
Mass:   rho_i (u_ice - v_I).n = rho_w (u_w - v_I).n   =>  (u_w - v_I).n = (rho_i/rho_w) m_n
Energy: rho_i L_f m_n = k_w d_n T_w - k_i d_n T_i      (d_n = derivative along n; sensible heat of the ice-side temperature gradient enters through k_i d_n T_i)
For ice at rest (u_ice = 0): v_I.n = -m_n and u_w.n = (rho_i/rho_w - 1) m_n  (small negative = liquid drawn toward the interface because meltwater is denser than ice)."""
from .units import RHO_I, RHO_W, LF, K_W, require
def melt_speed(dTw_dn, dTi_dn, k_w=K_W, k_i=2.1, rho_i=RHO_I, Lf=LF):
    """m_n [m/s] from the Stefan condition; temperate ice (T_i = T_m everywhere): dTi_dn = 0. Negative value = freezing (not supported by the caller unless explicitly allowed)."""
    return (k_w * dTw_dn - k_i * dTi_dn) / (rho_i * Lf)
def water_normal_velocity(m_n, rho_i=RHO_I, rho_w=RHO_W, u_ice_n=0.0):
    """u_w.n for an ice boundary moving with v_I.n = u_ice_n - m_n :  u_w.n = v_I.n + (rho_i/rho_w) m_n."""
    return (u_ice_n - m_n) + (rho_i / rho_w) * m_n
def meltwater_volume_rate(m_n, area, rho_i=RHO_I, rho_w=RHO_W): return (rho_i / rho_w) * m_n * area
def covered_interface_melt(q_conductive_from_debris, rho_i=RHO_I, Lf=LF):
    """debris-covered/dry interface: no stream contact; melt driven by the heat conducted from the debris side q [W/m2] (positive toward the ice) and the ice-side gradient; supplied by the caller (not modelled here)."""
    return q_conductive_from_debris / (rho_i * Lf)
