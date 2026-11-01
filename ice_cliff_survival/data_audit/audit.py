"""Reproducible audit of the two pinned Zenodo archives (raw data untouched; read via zipfile). python audit.py -> AUDIT.json, AUDIT.md.
Reports actual dates, site ids, coordinates, missingness, units, uncertainties, and product defects found. Never interpolates or fabricates forcing."""
import csv, datetime as dt, io, json, math, os, re, struct, sys, zipfile
HERE = os.path.dirname(os.path.abspath(__file__)); RAW = os.path.join(HERE, "raw")
K = zipfile.ZipFile(os.path.join(RAW, "Petersen_etal_2025_Kennicott_Glacier_Data.zip")); S = zipfile.ZipFile(os.path.join(RAW, "Petersen_etal_Streams_IceCliffs_DataRepo.zip"))
def names(z, key, ext=None): return [n for n in z.namelist() if "__MACOSX" not in n and key in n and not n.endswith("/") and (ext is None or n.endswith(ext))]
def rd(z, n): return z.read(n).decode("utf-8-sig", "replace")
def parse_dt(s): return dt.datetime.strptime(s.strip(), "%Y-%m-%d %H:%M:%S")
def backwasting():
    slope = {r["Stake"]: r for r in csv.DictReader(io.StringIO(rd(K, names(K, "Ice_Cliff_Backwasting", "ice_cliff_slope.csv")[0])))}
    aspect = {r["Stake"]: r for r in csv.DictReader(io.StringIO(rd(K, names(K, "Ice_Cliff_Backwasting", "ice_cliff_aspect.csv")[0])))}
    out = []
    for n in sorted(names(K, "Backwasting_Data", ".csv")):
        year = n.split("Backwasting_Data/")[1].split("/")[0]; site = n.rsplit("/", 1)[1][:-4]; rows = list(csv.DictReader(io.StringIO(rd(K, n)))); tot = len(rows); obs = []
        for r in rows:
            if r["Datetime"].strip() == "" or r["D"].strip() == "": continue
            obs.append((parse_dt(r["Datetime"]), parse_dt(r["Datetime_ADT"]), float(r["D"])))
        tz_ok = all(abs((u - a).total_seconds() / 3600 - 8) < 1e-9 for u, a, _ in obs); t = [o[0] for o in obs]; D = [o[2] for o in obs]
        gaps_h = [(b - a).total_seconds() / 3600 for a, b in zip(t, t[1:])]; dD = [b - a for a, b in zip(D, D[1:])]
        mono = all(x >= 0 for x in dD); sl = slope.get(site, {}).get(year); asp = aspect.get(site, {}).get(year)
        rates = [d / (g / 24) for d, g in zip(dD, gaps_h) if g > 0]
        out.append(dict(file=n.split("Petersen_etal_2025_Kennicott_Glacier_Data/")[1], year=int(year), site=site, rows_in_file=tot, observations=len(obs), blank_rows=tot - len(obs), first=str(t[0]) if t else None, last=str(t[-1]) if t else None, span_days=(t[-1] - t[0]).total_seconds() / 86400 if len(t) > 1 else None,
                        D_first_cm=D[0] if D else None, D_last_cm=D[-1] if D else None, D_monotone_nondecreasing=mono, utc_minus_adt_is_8h=tz_ok, gap_h_min=min(gaps_h) if gaps_h else None, gap_h_median=sorted(gaps_h)[len(gaps_h) // 2] if gaps_h else None, gap_h_max=max(gaps_h) if gaps_h else None,
                        slope_deg=float(sl) if sl not in (None, "") else None, aspect_deg=float(asp) if asp not in (None, "") else None, rate_cm_per_day_median=sorted(rates)[len(rates) // 2] if rates else None, rate_cm_per_day_max=max(rates) if rates else None,
                        regular_daily=bool(gaps_h) and max(gaps_h) < 40 and min(gaps_h) > 12))
    return out
def gauges():
    summ = list(csv.DictReader(io.StringIO(rd(S, names(S, "gauges_summary.csv")[0])))); res = []
    for r in summ:
        gid = r["Gauge #"]; date = r["Date"]; folders = sorted({n.split("/")[-2] for n in names(S, "Stream_Discharge_Measurements/gauge_" + gid + "_")}); day = date.split("-")[0]; fold = [f for f in folders if f.endswith("_07_" + day)]
        if not fold: res.append(dict(gauge=gid, exp=r["Experiment #"], date=date, status="no data folder")); continue
        base = [n for n in names(S, "Stream_Discharge_Measurements/" + fold[0] + "/")]; tr = list(csv.DictReader(io.StringIO(rd(S, [n for n in base if n.endswith("transect.csv")][0])))); L = float(r["Length (m)"])
        areas = []
        for n in sorted(n for n in base if re.search(r"prof\d+\.csv$", n)):
            pts = [tuple(map(float, l.split(","))) for l in rd(S, n).strip().splitlines() if l.strip()]; pts.sort(); a = sum(0.5 * (pts[i][1] + pts[i + 1][1]) * (pts[i + 1][0] - pts[i][0]) for i in range(len(pts) - 1)); areas.append(a)   # cm2 (distance cm x depth cm)
        tm = list(csv.DictReader(io.StringIO(rd(S, [n for n in base if n.endswith("time.csv")][0])))); use = [x for x in tm if x["Use"].strip() == "1"]; T = [float(x["Time"]) if x["Time"].strip() else float(x["Lap1"]) + float(x["Lap2"]) for x in use]
        Tt = [x for x in T if x > 0]; v = [L / x for x in Tt]; A = sum(areas) / len(areas) * 1e-4
        res.append(dict(gauge=gid, exp=r["Experiment #"], date=date, time_of_day=r["Time of Day"], conditions=r["Conditions"], lat=float(r["Y"]), lon=float(r["X"]), elev_m=float(r["Z"]), transect_length_m=L, profiles=len(areas), area_cm2=[round(a) for a in areas], mean_area_m2=A, floats_used=len(Tt), floats_total=len(tm), declared_good=r["Good_Datapoints"], float_velocity_m_s_mean=sum(v) / len(v), float_velocity_m_s_sd=(sum((x - sum(v) / len(v)) ** 2 for x in v) / max(len(v) - 1, 1)) ** 0.5, Q_surface_float_m3_s=A * sum(v) / len(v), status="ok"))
    return res
def dbf_count(z, n): b = z.read(n); return struct.unpack("<I", b[4:8])[0], struct.unpack("<H", b[8:10])[0]
def stats_csv():
    n = names(S, "ice_cliffs_stats.csv")[0]; rows = list(csv.DictReader(io.StringIO(rd(S, n)))); sw = sum(1 for r in rows if float(r["slope_min"]) > float(r["slope_max"]))
    return dict(rows=len(rows), columns=list(rows[0].keys()), rows_with_slope_min_gt_slope_max=sw, fraction_swapped=sw / len(rows), planform_area_total_m2=sum(float(r["area_planform"]) for r in rows), area_sloped_total_m2=sum(float(r["area_sloped"]) for r in rows),
                note="slope_min/slope_max appear to have swapped labels (min > max in most rows); do not use either as written without checking against slope_mean/median")
def profiles():
    p = names(S, "NPS_DEM_cliff_ramp_profile_samples", ".csv"); kinds = {}
    for n in p: k = n.rsplit("/", 1)[1].rsplit("_", 1)[0]; kinds[k] = kinds.get(k, 0) + 1
    ad = rd(S, names(S, "Ice_cliff_backwasting_ArcticDEM.csv")[0]).splitlines()
    return dict(ramp_profile_files=len(p), kinds=kinds, arcticdem_columns=ad[0], arcticdem_rows=len(ad) - 1, arcticdem_years=[2013, 2014, 2017, 2018])
def main():
    out = dict(kennicott_backwasting=backwasting(), gauges=gauges(), ice_cliffs_stats=stats_csv(), profiles=profiles(), dbf_records={n.rsplit("/", 1)[1]: dbf_count(S, n)[0] for n in names(S, ".dbf")})
    # other Kennicott products inventory (counts only)
    inv = {}
    for n in K.namelist():
        if "__MACOSX" in n or n.endswith("/") or n.endswith(".DS_Store"): continue
        grp = n.split("Petersen_etal_2025_Kennicott_Glacier_Data/")[1].split("/")[0]; inv.setdefault(grp, 0); inv[grp] += 1
    out["kennicott_inventory_files_by_folder"] = inv; bw = out["kennicott_backwasting"]; out["summary"] = dict(csv_files=len(bw), observations=sum(b["observations"] for b in bw), sites=sorted({b["site"] for b in bw}), site_years=len(bw), years=sorted({b["year"] for b in bw}),
        independent_cliffs_upper_bound=len({b["site"] for b in bw}), all_utc_adt_8h=all(b["utc_minus_adt_is_8h"] for b in bw), all_monotone=all(b["D_monotone_nondecreasing"] for b in bw), regular_daily_series=sum(b["regular_daily"] for b in bw))
    json.dump(out, open(os.path.join(HERE, "AUDIT.json"), "w"), indent=1, default=str); print(json.dumps(out["summary"], indent=1)); return out
if __name__ == "__main__": main()
