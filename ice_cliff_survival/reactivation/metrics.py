"""Metrics (definitions fixed in the freeze). contact: toe-route discharge fraction f_toe = Q_toe/Q_in >= theta (theta in THETAS; primary 1e-3).
loss time = first interpolated crossing below theta; temporary rewetting = f_toe back above theta after loss for < SUSTAIN_D days; post-removal reconnection = f_toe >= theta continuously for >= SUSTAIN_D days after removal following a loss before removal.
Local toe melt = ice volume lost in the TOE-LANE REACH (excludes entrance ramp); not whole-cliff or glacier mass loss."""
import numpy as np
THETAS = (1e-4, 1e-3, 1e-2, 5e-2); SUSTAIN_D = 1.0; PRIMARY = 1e-3
def arrays(run):
    r = run["records"]; g = lambda k: np.array([x.get(k, np.nan) for x in r], float); return dict(day=g("day"), f=g("f_toe"), Qt=g("Q_toe"), Qb=g("Q_byp"), melt_toe=g("cum_melt_toe_reach"), melt_byp=g("cum_melt_byp_reach"), crest_toe=g("crest_toe"), crest_byp=g("crest_byp"), wet_toe=g("wet_ice_area_toe"))
def events(day, f, theta, removal=None):
    lost = None; rewet = []; back = None; t_sus = None
    above = f >= theta
    for i in range(1, len(f)):
        if lost is None and above[i - 1] and not above[i]: lost = day[i - 1] + (day[i] - day[i - 1]) * (f[i - 1] - theta) / (f[i - 1] - f[i])
    if lost is None: return dict(loss=None, rewet_episodes=0, sustained_reconnection=None)
    # episodes of f>=theta after the loss
    i = int(np.searchsorted(day, lost)); eps = []; start = None
    for k in range(i, len(f)):
        if above[k] and start is None: start = day[k - 1] + (day[k] - day[k - 1]) * (theta - f[k - 1]) / (f[k] - f[k - 1]) if k > 0 and f[k] != f[k - 1] else day[k]
        if (not above[k] or k == len(f) - 1) and start is not None:
            end = day[k] if (above[k] and k == len(f) - 1) else day[k - 1]; eps.append((float(start), float(end), bool(above[k] and k == len(f) - 1))); start = None
    sust = [e for e in eps if (e[1] - e[0]) >= SUSTAIN_D or e[2]]; brief = [e for e in eps if e not in sust]
    rc = (max(sust[0][0], removal) if (sust and removal is not None) else (sust[0][0] if sust else None))      # reconnection cannot precede the removal (the bed update interval smears the jump)
    return dict(loss=float(lost), rewet_episodes=len(brief), sustained_reconnection=rc)
def summarize(run, theta=PRIMARY):
    a = arrays(run); ev = events(a["day"], a["f"], theta, run.get("removal_day")); Qin = 0.07; dt = np.diff(a["day"]) * 86400
    water_toe = float(np.sum(0.5 * (a["Qt"][1:] + a["Qt"][:-1]) * dt)); rem = run.get("removal_day"); post = float(np.mean(a["f"][a["day"] >= rem + 1.0])) if rem is not None and a["day"][-1] >= rem + 1 else None
    return dict(name=run.get("name"), layout=run["layout"], dx=run["dx"], dtm_h=run["dtm"] / 3600, removal_day=rem, theta=theta, **ev, f0=float(a["f"][0]), f_end=float(a["f"][-1]), delivered_toe_water_m3=water_toe, delivered_fraction=water_toe / (Qin * (a["day"][-1] * 86400)),
                toe_melt_all_m3=float(run["records"][-1]["cum_melt_toe"]), toe_melt_reach_m3=float(a["melt_toe"][-1]), byp_melt_reach_m3=float(a["melt_byp"][-1]), crest_toe_end=float(a["crest_toe"][-1]), crest_byp_end=float(a["crest_byp"][-1]), post_removal_mean_f=post,
                max_water_resid=float(np.nanmax(np.abs([x["water_resid"] for x in run["records"]]))), max_energy_resid=float(np.nanmax(np.abs([x["energy_resid_rel"] for x in run["records"]]))), wall_s=run["wall_s"], maxrss_MB=run["maxrss_MB"], cells=run["cells"])
def f_interp(run, days):
    a = arrays(run); return np.interp(days, a["day"], a["f"])

# ---- corrected recovery classification (v2) --------------------------------------------------------------------------------
RECOVERY_FULL, RECOVERY_PARTIAL = 0.5, 0.05
def classify(day, f, removal, theta=PRIMARY, f_ref=None):
    """access class relative to the bare-route share f_ref (ordinary access). maintained: never below theta before removal (a removal then cannot 'reconnect' anything); lost-recovered-full: R>=0.5; lost-recovered-partial: 0.05<=R<0.5; lost-not-recovered: R<0.05, where R = mean f over [removal+1 d, end] / f_ref. decline flag: f at end < 0.5 f(removal+1 d).
    For runs without removal: maintained / lost. (v1 'sustained reconnection' counted any f>=theta as reconnection and labelled runs that never lost contact as reconnected.)"""
    day = np.asarray(day); f = np.asarray(f); pre = f[day <= (removal if removal is not None else day[-1])]; lost_pre = bool(np.any(pre < theta))
    if removal is None: return dict(cls="lost" if np.any(f < theta) else "maintained", R=None, decline=False)
    if not lost_pre: cls = "maintained"; w = day >= removal + 1.0
    else: w = day >= removal + 1.0
    R = float(np.mean(f[w]) / f_ref) if (np.any(w) and f_ref) else None; f1 = float(np.interp(removal + 1.0, day, f)); decline = bool(f[-1] < 0.5 * f1) if f1 > theta else False
    if lost_pre:
        cls = "lost-recovered-full" if R >= RECOVERY_FULL else ("lost-recovered-partial" if R >= RECOVERY_PARTIAL else "lost-not-recovered")
    return dict(cls=cls, R=R, decline=decline, f_post1=f1)
