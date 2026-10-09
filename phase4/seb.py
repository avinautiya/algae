#!/usr/bin/env python
"""
Offline point surface energy balance (SEB) for bare glacier ice, driven by measured PROMICE KAN_M hourly
meteorology, used to turn an algal albedo reduction into ACTUAL melt (not "potential melt").

Energy balance at the surface (W m^-2, positive towards the surface):
    Q_M = SW_down (1 - alpha) + LW_down - eps sigma T_s^4 + H + LE          (G = 0, see assumptions)
  * if Q_M(T_s = 0 C) > 0 the surface melts at 0 C and the melt rate is Q_M / (L_f rho_w);
  * otherwise T_s < 0 C is solved from Q_M(T_s) = 0 and there is no melt.
Turbulent fluxes: bulk aerodynamic, H = rho_a c_p C U (T_a - T_s), LE = rho_a L_s C U (q_a - q_s(T_s)),
C = kappa^2 / ln(z/z0)^2 x stability factor from the bulk Richardson number (stable: (1 - 5 Ri)^2 for
Ri < 0.2, 0 above; unstable: (1 - 16 Ri)^0.5). z0 = 1 mm for ice (sensitivity x10, /10).

Assumptions (stated, with sensitivity where measurable):
  * no subsurface conduction (G = 0) and no cold-content memory: night-time cooling of the ice is not
    carried into the next day, so melt onset in the morning is slightly early (checked against measured
    surface temperature from the upward longwave);
  * eps = 0.98 for ice; measured LW_down, T_a, q_a, U, p, SW_down;
  * bare-ice periods only (PROMICE snow_height < 0.02 m).

Validation against independent observations: hourly modelled melt with the MEASURED albedo, summed over
bare-ice periods, vs the measured ice-surface lowering (z_ice_surf) x rho_ice / rho_w.

Algal coupling without double counting: the measured station albedo already contains whatever algae
were there. Paired runs therefore use a MODEL albedo for the same surface with and without algae:
    off: alpha_clean(t)                 on: alpha_clean(t) - d_alpha_algae(t)
with identical meteorology; algal melt = melt_on - melt_off. The surface-temperature, longwave and
turbulent responses and the 'no melt when cold' threshold are thus included, which the potential-melt
estimate (RF x 86400 / L_f) ignores.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "phase2"))

import empirical_data as ED  # noqa: E402

SIGMA_SB = 5.670374419e-8
EPS_ICE = 0.98
L_F = 3.34e5            # J kg^-1
L_S = 2.834e6           # J kg^-1 (sublimation; deposition)
L_V = 2.501e6
CP = 1005.0
KAPPA = 0.4
RHO_W = 1000.0
T0 = 273.15
FORCING = "promice_KAN_M_hour_JJA_2016_2019_seb.csv"


def q_sat(T_k, p_hpa):
    """Saturation specific humidity over ice (kg/kg), Goff-Gratch-like Magnus fit."""
    Tc = T_k - T0
    es = 6.112 * np.exp(22.46 * Tc / (272.62 + Tc))          # hPa, over ice
    return 0.622 * es / (p_hpa - 0.378 * es)


def turbulent(T_s, T_a, q_a, U, p_hpa, z=2.5, z0=1e-3):
    """Sensible and latent heat flux (W m^-2, positive towards the surface)."""
    rho_a = p_hpa * 100.0 / (287.05 * T_a)
    U = np.maximum(U, 0.1)
    C = KAPPA ** 2 / np.log(z / z0) ** 2
    Ri = 9.81 * z * (T_a - T_s) / (0.5 * (T_a + T_s) * U ** 2)
    f = np.where(Ri >= 0, np.where(Ri < 0.2, (1 - 5 * Ri) ** 2, 0.0), np.sqrt(1 - 16 * np.minimum(Ri, 0)))
    H = rho_a * CP * C * f * U * (T_a - T_s)
    L = np.where(T_s >= T0, L_V, L_S)
    LE = rho_a * L * C * f * U * (q_a - q_sat(T_s, p_hpa))
    return H, LE


def balance(T_s, sw_net, lw_in, T_a, q_a, U, p, z, z0):
    H, LE = turbulent(T_s, T_a, q_a, U, p, z, z0)
    return sw_net + lw_in - EPS_ICE * SIGMA_SB * T_s ** 4 - (1 - EPS_ICE) * 0 + H + LE, H, LE


def solve_hour(sw_net, lw_in, T_a, q_a, U, p, z=2.5, z0=1e-3):
    """Vectorised over hours. Returns (melt_flux W m^-2, T_s K, H, LE)."""
    sw_net, lw_in, T_a, q_a, U, p, z = (np.asarray(v, float) for v in (sw_net, lw_in, T_a, q_a, U, p, z))
    q0, H0, LE0 = balance(np.full_like(T_a, T0), sw_net, lw_in, T_a, q_a, U, p, z, z0)
    melt = np.maximum(q0, 0.0)
    T_s = np.full_like(T_a, T0)
    cold = q0 < 0
    if cold.any():                                        # bisection for T_s < 0 C with zero net flux
        lo = np.full(cold.sum(), 200.0)
        hi = np.full(cold.sum(), T0)
        args = [v[cold] for v in (sw_net, lw_in, T_a, q_a, U, p, z)]
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            f, _, _ = balance(mid, args[0], args[1], args[2], args[3], args[4], args[5], args[6], z0)
            hi = np.where(f < 0, mid, hi)
            lo = np.where(f < 0, lo, mid)
        T_s[cold] = 0.5 * (lo + hi)
    _, H, LE = balance(T_s, sw_net, lw_in, T_a, q_a, U, p, z, z0)
    return melt, T_s, H, LE


C_ICE = 2097.0           # J kg^-1 K^-1
RHO_ICE = 900.0


def solve_series(sw_net, lw_in, T_a, q_a, U, p, z, z0=1e-3, slab_m=0.1, dt=3600.0):
    """Hourly time stepping with a surface ice slab of thickness slab_m carrying cold content: the slab
    temperature is the surface temperature; net flux warms/cools it; melt only when it is at 0 C and the
    flux is positive. slab_m = 0 reduces to the instantaneous balance (solve_hour)."""
    if slab_m <= 0:
        return solve_hour(sw_net, lw_in, T_a, q_a, U, p, z, z0)
    n = len(T_a)
    cap = C_ICE * RHO_ICE * slab_m
    Ts = T0
    melt, T_s, Hs, LEs = (np.zeros(n) for _ in range(4))
    for i in range(n):
        q, H, LE = balance(np.array([Ts]), sw_net[i], lw_in[i], T_a[i], q_a[i], U[i], p[i], z[i], z0)
        q = float(q[0])
        Tn = Ts + q * dt / cap
        if Tn > T0:
            melt[i] = (Tn - T0) * cap / dt
            Tn = T0
        Ts = max(Tn, 200.0)
        T_s[i], Hs[i], LEs[i] = Ts, float(H[0]), float(LE[0])
    return melt, T_s, Hs, LEs


def load_forcing(year=None):
    d = pd.read_csv(ED.path(FORCING), parse_dates=["time"])
    if year is not None:
        d = d[d.time.dt.year == year]
    d = d.copy()
    # shortwave: tilt-corrected where available (2019 has only the uncorrected dsr)
    d["sw_down"] = d.dsr_cor.where(d.dsr_cor.notna(), d.dsr)
    d["sw_source"] = np.where(d.dsr_cor.notna(), "dsr_cor", "dsr")
    d["T_a"] = d.t_u + T0
    d["q_a"] = d.qh_u / 1000.0
    d["z_meas"] = d.z_boom_u.clip(1.0, 4.0).fillna(2.5)
    return d.reset_index(drop=True)


def run(forcing, albedo, z0=1e-3, slab_m=0.1):
    """Hourly SEB with a given albedo series (array or scalar). Returns a DataFrame. Missing forcing
    hours are filled by interpolation (at most 6 h) so the slab can be stepped continuously."""
    f = forcing
    alb = np.broadcast_to(np.asarray(albedo, float), (len(f),))
    cols = {k: f[k].interpolate(limit=6, limit_direction="both").to_numpy(float)
            for k in ("sw_down", "dlr", "T_a", "q_a", "wspd_u", "p_u", "z_meas")}
    cols["sw_down"] = np.nan_to_num(cols["sw_down"], nan=0.0)
    melt, T_s, H, LE = solve_series(cols["sw_down"] * (1 - alb), cols["dlr"], cols["T_a"], cols["q_a"], cols["wspd_u"],
                                    cols["p_u"], cols["z_meas"], z0, slab_m)
    return pd.DataFrame(dict(time=f.time, melt_w=melt, melt_mwe=melt * 3600.0 / (L_F * RHO_W), T_s=T_s, H=H, LE=LE,
                             sw_net=cols["sw_down"] * (1 - alb)))


def observed_ablation(forcing, rho_ice=900.0):
    """Hourly surface lowering (m w.e.) from the PROMICE ice-surface height on bare-ice hours."""
    z = forcing.z_ice_surf.interpolate(limit=6)
    dz = -z.diff().fillna(0.0)
    return dz * rho_ice / RHO_W


def validate(year, z0=1e-3, rho_ice=900.0, slab_m=0.1):
    """Modelled melt with the MEASURED albedo vs measured ablation, over bare-ice days with complete data
    (daily sums; days with snow, missing albedo or missing height excluded)."""
    f = load_forcing(year)
    ok = f[["sw_down", "dlr", "T_a", "q_a", "wspd_u", "p_u"]].notna().all(axis=1)
    # PROMICE albedo is reported only for high sun; each hour gets that day's mean measured albedo
    day = f.time.dt.floor("D")
    alb_day = f.albedo.groupby(day).transform("mean")
    ok &= alb_day.notna()
    m = run(f, alb_day.fillna(0.6).clip(0.05, 0.95), z0, slab_m)
    obs = observed_ablation(f, rho_ice)
    g = pd.DataFrame(dict(day=day, model=m.melt_mwe, obs=obs, ok=ok, snow=f.snow_height.fillna(0) > 0.02,
                          T_s_mod=m.T_s, T_s_meas=f.t_surf + T0, has_z=f.z_ice_surf.notna()))
    dd = g.groupby("day").agg(model=("model", "sum"), obs=("obs", "sum"), ok=("ok", "mean"), snow=("snow", "max"),
                              has_z=("has_z", "mean"))
    keep = (dd.ok > 0.95) & (~dd.snow) & (dd.has_z > 0.9)
    dd = dd[keep]
    ts = g.dropna(subset=["T_s_meas"])
    return dict(year=year, n_days=int(len(dd)), model_total_mwe=float(dd.model.sum()), obs_total_mwe=float(dd.obs.sum()),
                daily_rmse_mwe=float(np.sqrt(np.mean((dd.model - dd.obs) ** 2))),
                daily_r=float(np.corrcoef(dd.model, dd.obs)[0, 1]) if len(dd) > 2 else np.nan,
                T_s_bias_K=float((ts.T_s_mod - ts.T_s_meas).mean()),
                T_s_rmse_K=float(np.sqrt(np.mean((ts.T_s_mod - ts.T_s_meas) ** 2))),
                sw_source=str(f.sw_source.iloc[0])), dd


def paired_algae(year, alpha_clean, dalpha, z0=1e-3, slab_m=0.1):
    """Paired SEB runs with and without algae on a model surface: alpha_clean (scalar or hourly) and the
    algal albedo reduction d_alpha (scalar or hourly). Returns totals (m w.e.) and the potential-melt
    estimate from the same albedo change for comparison."""
    f = load_forcing(year)
    ok = f[["sw_down", "dlr", "T_a", "q_a", "wspd_u", "p_u"]].notna().all(axis=1).to_numpy()
    f = f[ok].reset_index(drop=True)
    ac = np.broadcast_to(np.asarray(alpha_clean, float), (len(f),)) if np.ndim(alpha_clean) == 0 else np.asarray(alpha_clean)[ok]
    da = np.broadcast_to(np.asarray(dalpha, float), (len(f),)) if np.ndim(dalpha) == 0 else np.asarray(dalpha)[ok]
    off = run(f, ac, z0, slab_m)
    on = run(f, ac - da, z0, slab_m)
    potential = float(np.sum(f.sw_down * da) * 3600.0 / (L_F * RHO_W))
    return dict(year=year, hours=int(len(f)), melt_off_mwe=float(off.melt_mwe.sum()), melt_on_mwe=float(on.melt_mwe.sum()),
                algal_melt_mwe=float(on.melt_mwe.sum() - off.melt_mwe.sum()), potential_algal_melt_mwe=potential)


def algae_2019(summary_csv, slab_m=0.1, z0=1e-3, method="rf_ours"):
    """Paired SEB runs for the 2019 sampled dates, KAN_M meteorology, S6 retrieved surface:
    alpha_on(t) = retrieved mean bare-ice BBA (algae included), alpha_off = alpha_on + d_alpha (algae
    removed), both linearly interpolated in time between the sampled dates (no extrapolation: the run is
    restricted to the first..last sampled date)."""
    t = pd.read_csv(summary_csv).sort_values("date")
    f = load_forcing(2019)
    d0, d1 = pd.Timestamp(t.date.iloc[0]), pd.Timestamp(t.date.iloc[-1]) + pd.Timedelta(days=1)
    f = f[(f.time >= d0) & (f.time < d1)].reset_index(drop=True)
    x = (f.time - d0).dt.total_seconds().to_numpy() / 86400.0
    xd = ((pd.to_datetime(t.date) - d0).dt.total_seconds() / 86400.0 + 0.5).to_numpy()
    a_on = np.interp(x, xd, t.bba_mean.to_numpy())
    da = np.interp(x, xd, t[f"dalpha_{method}"].to_numpy())
    on = run(f, a_on, z0, slab_m)
    off = run(f, a_on + da, z0, slab_m)
    potential = float(np.sum(np.nan_to_num(f.sw_down.to_numpy()) * da) * 3600.0 / (L_F * RHO_W))
    return dict(method=method, start=str(d0.date()), end=str((d1 - pd.Timedelta(days=1)).date()), days=float(len(f) / 24),
                slab_m=slab_m, z0=z0, melt_on_mwe=float(on.melt_mwe.sum()), melt_off_mwe=float(off.melt_mwe.sum()),
                algal_melt_mwe=float(on.melt_mwe.sum() - off.melt_mwe.sum()), potential_algal_melt_mwe=potential,
                actual_over_potential=float((on.melt_mwe.sum() - off.melt_mwe.sum()) / potential))
