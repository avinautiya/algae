"""Independent budget reconstruction (does not call the stepping code): from a recorded history of source/export rates, recompute the stock change by trapezoid/rectangle integration and compare."""
def grain_budget(Vs0, t, sources_series, Q_out_series):
    """t: list of step end times starting at t0 = t[0] (so len = n+1); sources_series: list (n) of dict(S_up,...) mid-step rates; Q_out_series: n rates. Piecewise-constant rates => exact sum."""
    Vs = Vs0
    for i in range(len(t) - 1):
        dt = t[i + 1] - t[i]; Vs += dt * (sum(sources_series[i].values()) - Q_out_series[i])
    return Vs
def energy_budget_water(records):
    """sum of per-step residuals of sensible + dissipation - export - heat_to_ice - storage change; returns (sum, sum of |terms|) for a relative denominator."""
    res = sum(r.water_energy_residual() for r in records); den = sum(abs(r.sensible_in_J) + abs(r.sensible_out_J) + abs(r.q_wall_W * (r.t1 - r.t0)) for r in records); return res, den
