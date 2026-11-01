"""Observation operators: what a model output can be compared with, and the geometry needed.
Stake record: D(t) = emerged HORIZONTAL stake length [cm] = horizontal distance between the face and the stake collar position. For a planar face of inclination theta from the horizontal retreating by normal distance M·dt,
  d(D)/dt = (normal recession speed) / sin(theta)   =>   normal recession speed = dD/dt * sin(theta).   (theta = 90 deg: vertical face, horizontal = normal.)
The conversion needs the face slope AT THE STAKE AT THAT TIME; the archive provides one slope per site-year (± unknown), so the normal rate carries that uncertainty. It is NOT exposed-ice melt energy: it includes undercut/collapse/shade effects and the stake-face geometry."""
import math
def normal_recession_from_horizontal(dD_dt, slope_deg):
    s = math.sin(math.radians(slope_deg))
    if not 0 < slope_deg <= 90: raise ValueError("face slope must be in (0, 90] degrees")
    return dD_dt * s
def horizontal_from_normal(M_normal, slope_deg): return M_normal / math.sin(math.radians(slope_deg))
def interval_rates(t_utc, D_cm, sigma_D=2.0, sigma_t_h=1.0):
    """per-interval horizontal rates [cm/day] with first-order uncertainty from ±sigma_D at both ends and ±sigma_t on the interval length (independent errors, quadrature). No interpolation; irregular intervals kept."""
    out = []
    for i in range(len(D_cm) - 1):
        dt_d = (t_utc[i + 1] - t_utc[i]).total_seconds() / 86400.0
        if dt_d <= 0: continue
        r = (D_cm[i + 1] - D_cm[i]) / dt_d; sig = math.sqrt((math.sqrt(2) * sigma_D / dt_d) ** 2 + (r * (math.sqrt(2) * sigma_t_h / 24.0) / dt_d) ** 2)
        out.append(dict(t0=t_utc[i], t1=t_utc[i + 1], days=dt_d, rate_cm_day=r, sigma_cm_day=sig))
    return out
