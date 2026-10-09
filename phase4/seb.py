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
    """Net surface flux (W m^-2). Longwave: the surface absorbs eps LW_down (the rest, (1 - eps) LW_down,
    is reflected - Kirchhoff) and emits eps sigma T_s^4, so LW_net = eps (LW_down - sigma T_s^4)."""
    H, LE = turbulent(T_s, T_a, q_a, U, p, z, z0)
    return sw_net + EPS_ICE * (lw_in - SIGMA_SB * T_s ** 4) + H + LE, H, LE


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


def solve_series(sw_net, lw_in, T_a, q_a, U, p, z, z0=1e-3, slab_m=0.1, dt=3600.0, n_sub=1, restart=None,
                 return_budget=False):
    """Time stepping with a surface ice slab of thickness slab_m carrying cold content (explicit Euler,
    n_sub sub-steps per forcing interval, forcing held constant over the interval - PROMICE values are
    interval means). The slab temperature is the surface temperature; net flux warms/cools it; melt only
    when it is at 0 C and the flux is positive. dt: scalar or per-step interval lengths (s). restart:
    boolean per step - re-initialise the slab (after a data gap) at min(T_a, 0 C). NaN forcing gives NaN
    output for that step (no filling). slab_m = 0 reduces to the instantaneous balance (solve_hour).
    With return_budget, also returns the energy bookkeeping (J m^-2): sum of net flux x dt, melt energy,
    change of slab heat content and energy discarded at restarts."""
    if slab_m <= 0:
        out = solve_hour(sw_net, lw_in, T_a, q_a, U, p, z, z0)
        return (*out, None) if return_budget else out
    n = len(T_a)
    dt = np.broadcast_to(np.asarray(dt, float), (n,))
    restart = np.zeros(n, bool) if restart is None else np.asarray(restart, bool)
    cap = C_ICE * RHO_ICE * slab_m
    Ts = T0
    melt, T_s, Hs, LEs = (np.full(n, np.nan) for _ in range(4))
    e_in = e_melt = e_restart = 0.0
    T_start = Ts
    for i in range(n):
        args = (sw_net[i], lw_in[i], T_a[i], q_a[i], U[i], p[i], z[i])
        if not np.all(np.isfinite(args)):
            continue
        if restart[i]:
            Tn = min(T_a[i], T0)
            e_restart += cap * (Tn - Ts)
            Ts = Tn
        h = dt[i] / n_sub
        m_e = 0.0
        Hm = LEm = 0.0
        for _ in range(n_sub):
            q, H, LE = balance(np.array([Ts]), *args[:6], args[6], z0)
            q = float(q[0])
            e_in += q * h
            Tn = Ts + q * h / cap
            if Tn > T0:
                m_e += (Tn - T0) * cap
                Tn = T0
            if Tn < 150.0:
                raise FloatingPointError(f"slab temperature {Tn:.1f} K at step {i}: unphysical forcing")
            Ts = Tn
            Hm += float(H[0]) / n_sub
            LEm += float(LE[0]) / n_sub
        e_melt += m_e
        melt[i], T_s[i], Hs[i], LEs[i] = m_e / dt[i], Ts, Hm, LEm
    if not return_budget:
        return melt, T_s, Hs, LEs
    budget = dict(energy_in=e_in, energy_melt=e_melt, energy_slab_change=cap * (Ts - T_start) - e_restart,
                  energy_restarts=e_restart)
    budget["closure_residual"] = budget["energy_in"] - budget["energy_melt"] - budget["energy_slab_change"]
    return melt, T_s, Hs, LEs, budget


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


REQUIRED = ("sw_down", "dlr", "T_a", "q_a", "wspd_u", "p_u", "z_meas")
SHORT_GAP_H = 2          # gaps up to this length are interpolated for the non-shortwave forcing only


def prepare(forcing, interp_short_gaps=True):
    """Put the forcing on its own complete hourly time axis (the actual timestamps; missing hours become
    rows of NaN rather than being skipped), interpolate gaps <= SHORT_GAP_H h for LW, T, q, U, p only
    (shortwave is never filled), and mark steps that restart after a longer gap."""
    f = forcing.set_index("time").sort_index()
    if f.index.has_duplicates:
        raise ValueError("duplicate timestamps in forcing")
    full = pd.date_range(f.index[0], f.index[-1], freq="h")
    f = f.reindex(full)
    f.index.name = "time"
    filled = {}
    for k in REQUIRED:
        if k != "sw_down" and interp_short_gaps:
            before = f[k].isna()
            run_id = (before != before.shift()).cumsum()
            run_len = before.groupby(run_id).transform("sum")
            short = before & (run_len <= SHORT_GAP_H)          # only gaps of <= SHORT_GAP_H hours, whole
            interp = f[k].interpolate(limit_area="inside")
            f[k] = f[k].where(~short, interp)
            filled[k] = int((before & f[k].notna()).sum())
    ok = f[list(REQUIRED)].notna().all(axis=1).to_numpy()
    prev_ok = np.r_[True, ok[:-1]]
    restart = ok & ~prev_ok
    f = f.reset_index()
    f.attrs.update(filled_hours=filled, missing_hours=int((~ok).sum()), restarts=int(restart.sum()))
    return f, ok, restart


def check_albedo(alb, ok):
    a = np.asarray(alb, float)
    bad = ok & ~(np.isfinite(a) & (a >= 0.0) & (a <= 1.0))
    if bad.any():
        raise ValueError(f"albedo missing or outside [0, 1] on {int(bad.sum())} hours with forcing")


def run(forcing, albedo, z0=1e-3, slab_m=0.1, n_sub=1, return_budget=False, sw_subsurface_frac=0.0):
    """SEB on the forcing's actual hourly time axis. albedo: scalar, or an array aligned with `forcing`
    rows (it is carried through `prepare` by timestamp). Hours with missing forcing (and the shortwave is
    never filled) give NaN; the slab restarts after such gaps. Albedo must lie in [0, 1] wherever forcing
    exists (raises otherwise)."""
    f0 = forcing.copy()
    f0["_alb"] = np.broadcast_to(np.asarray(albedo, float), (len(f0),))
    f, ok, restart = prepare(f0)
    alb = f["_alb"].to_numpy(float)
    check_albedo(alb, ok)
    cols = {k: f[k].to_numpy(float) for k in REQUIRED}
    sw_net = cols["sw_down"] * (1 - alb)
    # a fraction of the absorbed shortwave deposited BELOW the surface (penetration into bubbly ice / the
    # weathering crust) is withheld from the surface balance and reported as subsurface energy; it melts
    # ice internally without (immediately) lowering the surface - a tested hypothesis, default 0
    if not 0.0 <= sw_subsurface_frac < 1.0:
        raise ValueError("sw_subsurface_frac must be in [0, 1)")
    sw_sub = sw_net * sw_subsurface_frac
    sw_net = sw_net - sw_sub
    res = solve_series(sw_net, cols["dlr"], cols["T_a"], cols["q_a"], cols["wspd_u"], cols["p_u"], cols["z_meas"],
                       z0, slab_m, 3600.0, n_sub, restart, return_budget=return_budget)
    melt, T_s, H, LE = res[:4]
    out = pd.DataFrame(dict(time=f.time, valid=ok, restart=restart, melt_w=melt, melt_mwe=melt * 3600.0 / (L_F * RHO_W),
                            T_s=T_s, H=H, LE=LE, sw_net=sw_net, sw_subsurface=sw_sub,
                            subsurface_melt_mwe=sw_sub * 3600.0 / (L_F * RHO_W),
                            lw_net=EPS_ICE * (cols["dlr"] - SIGMA_SB * T_s ** 4)))
    out.attrs.update(f.attrs)
    if return_budget:
        out.attrs["budget"] = res[4]
    return out


def observed_ablation(forcing, rho_ice=900.0, record="z_ice_surf"):
    """Hourly surface lowering (m w.e.) from a PROMICE height record (default the combined ice-surface
    height; z_pt_cor, z_stake_cor or z_boom_u for sensor cross-checks; the boom/stake distances increase
    as the surface lowers, so their sign is flipped). Gaps up to 6 h are interpolated."""
    z = forcing[record].interpolate(limit=6, limit_area="inside")
    if record in ("z_boom_u", "z_stake_cor", "z_stake"):
        z = -z
    dz = -z.diff().fillna(0.0)
    return dz * rho_ice / RHO_W


def validate(year, z0=1e-3, rho_ice=900.0, slab_m=0.1, n_sub=1, ablation="z_ice_surf", sw_subsurface_frac=0.0):
    """Modelled melt with the MEASURED albedo vs measured ablation, daily sums over bare-ice days with
    complete forcing, albedo and height (days with snow or any missing hour excluded). PROMICE albedo is
    reported only for high sun: each hour gets that day's mean measured albedo; days without any albedo
    are excluded explicitly (their shortwave is set missing, never filled). `ablation` selects the height
    record (z_ice_surf, or the separate sensors z_pt_cor / z_stake_cor / z_boom_u for cross-checks)."""
    f = load_forcing(year)
    day = f.time.dt.floor("D")
    alb_day = f.albedo.groupby(day).transform("mean")
    f.loc[alb_day.isna(), "sw_down"] = np.nan
    m = run(f, alb_day.fillna(0.5).to_numpy(), z0, slab_m, n_sub,       # fill value never used: SW is missing
            sw_subsurface_frac=sw_subsurface_frac)
    f2 = f.set_index("time").reindex(m.time).reset_index()
    obs = observed_ablation(f2, rho_ice, ablation)
    dayf = m.time.dt.floor("D")
    g = pd.DataFrame(dict(day=dayf, model=m.melt_mwe, obs=obs, ok=m.valid, snow=f2.snow_height.fillna(0) > 0.02,
                          T_s_mod=m.T_s, T_s_meas=f2.t_surf + T0, has_z=f2[ablation].notna()))
    dd = g.groupby("day").agg(model=("model", "sum"), obs=("obs", "sum"), ok=("ok", "mean"), snow=("snow", "max"),
                              has_z=("has_z", "mean"))
    dd = dd[(dd.ok == 1.0) & (~dd.snow) & (dd.has_z > 0.9)]
    ts = g[g.ok].dropna(subset=["T_s_meas", "T_s_mod"])
    return dict(year=year, ablation_record=ablation, n_days=int(len(dd)), model_total_mwe=float(dd.model.sum()),
                obs_total_mwe=float(dd.obs.sum()), ratio=float(dd.model.sum() / dd.obs.sum()) if dd.obs.sum() else np.nan,
                daily_rmse_mwe=float(np.sqrt(np.mean((dd.model - dd.obs) ** 2))) if len(dd) else np.nan,
                daily_r=float(np.corrcoef(dd.model, dd.obs)[0, 1]) if len(dd) > 2 else np.nan,
                T_s_bias_K=float((ts.T_s_mod - ts.T_s_meas).mean()),
                T_s_rmse_K=float(np.sqrt(np.mean((ts.T_s_mod - ts.T_s_meas) ** 2))),
                sw_source=str(f.sw_source.iloc[0]), missing_hours=m.attrs["missing_hours"],
                restarts=m.attrs["restarts"]), dd


def paired_algae(year, alpha_on, alpha_off, z0=1e-3, slab_m=0.1, n_sub=1, forcing=None, sw_subsurface_frac=0.0):
    """Paired SEB runs on identical forcing (actual hourly axis, gaps kept): surface WITH algae (alpha_on)
    and WITHOUT (alpha_off), scalars or arrays aligned with the forcing rows. Returns totals (m w.e.) over
    hours valid in both runs and the 'potential' estimate sum(SW (alpha_off - alpha_on)) / L_f for
    comparison. The difference is a MODELLED melt increment, conditional on the SEB, the forcing and the
    albedo inputs."""
    f = load_forcing(year) if forcing is None else forcing
    on = run(f, alpha_on, z0, slab_m, n_sub, sw_subsurface_frac=sw_subsurface_frac)
    off = run(f, alpha_off, z0, slab_m, n_sub, sw_subsurface_frac=sw_subsurface_frac)
    ok = on.valid & off.valid
    f2 = f.assign(_a_on=np.broadcast_to(np.asarray(alpha_on, float), (len(f),)),
                  _a_off=np.broadcast_to(np.asarray(alpha_off, float), (len(f),)))
    f2 = f2.set_index("time").reindex(on.time)
    potential = float(np.nansum((f2.sw_down * (f2._a_off - f2._a_on)).to_numpy()[ok.to_numpy()]) * 3600.0 / (L_F * RHO_W))
    m_on, m_off = float(on.melt_mwe[ok].sum()), float(off.melt_mwe[ok].sum())
    return dict(year=year, valid_hours=int(ok.sum()), missing_hours=int((~ok).sum()), melt_on_mwe=m_on,
                melt_off_mwe=m_off, modelled_algal_melt_increment_mwe=m_on - m_off,
                potential_algal_melt_mwe=potential,
                increment_over_potential=(m_on - m_off) / potential if potential else np.nan)


def algae_2019(summary_csv, slab_m=0.1, z0=1e-3, method="rf_ours", n_sub=1, sw_subsurface_frac=0.0):
    """Paired SEB runs for the 2019 sampled dates (KAN_M meteorology, retrieved S6 surface), each optical
    model with ITS OWN albedos: alpha_on(t) = that model's retrieved mean bare-ice BBA (algae included),
    alpha_off = alpha_on + that model's algal albedo reduction; both linearly interpolated in time between
    the sampled dates, run restricted to first..last sampled date (no extrapolation)."""
    t = pd.read_csv(summary_csv).sort_values("date")
    lab = "ours" if method == "rf_ours" else "tierA"
    bcol = "bba_mean" if lab == "ours" else "bba_tierA_mean"
    if bcol not in t:
        raise ValueError(f"{summary_csv} has no {bcol}: re-aggregate with multi_scene.py (own-albedo runs need it)")
    f = load_forcing(2019)
    d0, d1 = pd.Timestamp(t.date.iloc[0]), pd.Timestamp(t.date.iloc[-1]) + pd.Timedelta(days=1)
    f = f[(f.time >= d0) & (f.time < d1)].reset_index(drop=True)
    x = (f.time - d0).dt.total_seconds().to_numpy() / 86400.0
    xd = ((pd.to_datetime(t.date) - d0).dt.total_seconds() / 86400.0 + 0.5).to_numpy()
    a_on = np.interp(x, xd, t[bcol].to_numpy())
    da = np.interp(x, xd, t[f"dalpha_{method}"].to_numpy())
    r = paired_algae(2019, a_on, a_on + da, z0, slab_m, n_sub, forcing=f, sw_subsurface_frac=sw_subsurface_frac)
    r.update(sw_subsurface_frac=sw_subsurface_frac, method=method, albedo_source=bcol, start=str(d0.date()), end=str((d1 - pd.Timedelta(days=1)).date()),
             slab_m=slab_m, z0=z0, n_sub=n_sub)
    return r
