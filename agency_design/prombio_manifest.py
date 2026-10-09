"""PROMBIO matching manifest: one row per biological observation, with candidate satellite and station
radiation matches. No coordinates, timestamps or labels are invented.

Roles are frozen in docs/new_data_inventory.md (2021/2023 development; 2024 final test, withheld).
The manifest records MATCHING information only: radiation and imagery availability, never model
predictions or scores.

    python common/run_budgeted.py --background --mem-mb 600 --name prombio_manifest -- \
        python3 -m agency_design.prombio_manifest --outdir records/prombio_matching
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from agency_design.prombio import load_2024, observation

DATA = Path(__file__).parent / "data"
RAW = Path(__file__).resolve().parents[1] / "data" / "promice_raw"
SAT_WINDOW_DAYS = 3
STRATA_2123 = {"radiometer (a)": "under radiometer", "representative (b)": "average ice in the area",
               "dark (c)": "dark ice"}
PROVIDER_QC_EXCLUDED = {("2023", "KAN_L"), ("2023", "KAN_M")}   # Anesio 2024 TR314 §2: KAN samples not of good quality
MISSING_SAMPLE_NOTES = ("no sample", "no liquid", "missing")


def load_2021_2023():
    """Rows x strata of the 2021-2023 workbook; min/max/average summary columns are NOT observations."""
    import openpyxl
    path = DATA / "prombio_2021_2023_original.xlsx"
    meta = json.loads((DATA / "prombio_2021_2023_metadata.json").read_text())
    ver = meta["data"]["latestVersion"]
    src = next(f["dataFile"] for f in ver["files"] if f["dataFile"]["filename"].endswith(".xlsx"))
    if hashlib.md5(path.read_bytes()).hexdigest() != src["md5"]:
        raise ValueError("provider checksum mismatch")
    ws = openpyxl.load_workbook(path, data_only=True).worksheets[0]
    rows = list(ws.iter_rows(values_only=True))
    top, head = rows[0], rows[1]
    groups = {}
    cur = None
    for j, v in enumerate(top):
        if v in STRATA_2123:
            cur = STRATA_2123[v]
        elif v is not None:
            cur = None
        if cur and head[j] in ("cells/ml", "approx. abiotic load**", "notes"):
            groups.setdefault(cur, {})[head[j]] = j
    i_lat, i_lon, i_samp = head.index("latitude"), head.index("longitude"), head.index("Sampling")
    out = []
    for ln, r in enumerate(rows[2:], start=3):
        if r[1] is None or not str(r[1])[:3].isupper() or r[i_samp] is None:
            continue                                            # footnotes
        station = str(r[1]).strip().replace("-", "_")
        date = r[0].date().isoformat() if r[0] is not None else None
        year = str(r[i_samp]).split()[-1]
        for stratum, cols in groups.items():
            raw = r[cols["cells/ml"]]
            ob = observation("" if raw is None else str(raw))
            out.append(dict(obs_id=f"PROMBIO{year}_{station}_{stratum.split()[0]}", source="PROMBIO 2021-2023",
                            source_line=ln, campaign=f"PROMBIO {year}", year=year, station=station, date=date,
                            sampling_stratum=stratum, sample_tool="ice screw (report §2, 2023)" if year == "2023" else "not stated",
                            obs_state=ob["kind"], ice_cells_ml=ob["value"], upper_bound=ob["upper_bound"], raw_value=ob["raw"],
                            abiotic_load=r[cols.get("approx. abiotic load**")], note=r[cols.get("notes")],
                            station_lat=r[i_lat], station_lon=r[i_lon],
                            doi=meta["data"]["persistentUrl"], version=f"{ver['versionNumber']}.{ver['versionMinorNumber']}",
                            sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    return out


def load_2024_rows():
    recs, prov = load_2024()
    out = []
    for r in recs:
        ob = r["ice_cells_ml"]
        note = (r["raw"].get("other observations") or "").strip()
        state = ob["kind"]
        if state == "unresolved_annotation" and any(k in note.lower() for k in MISSING_SAMPLE_NOTES):
            state = "missing_sample"
        out.append(dict(obs_id=r["sample_id"], source="PROMBIO 2024", source_line=r["source_line"], campaign="PROMBIO 2024",
                        year="2024", station=r["station_area"], date=r["date"], sampling_stratum=r["sampling_stratum"].lower(),
                        sample_tool={"C": "chisel", "IS": "ice screw"}.get(r["raw"]["Chisel (C) vs ice screw (IS)"].strip(), "not stated"),
                        obs_state=state, ice_cells_ml=ob["value"], upper_bound=ob["upper_bound"], raw_value=ob["raw"],
                        abiotic_load=r["raw"]["mineral/cell ratio (qualitative measurement)"], note=note,
                        station_lat=None, station_lon=None, doi=prov["doi"], version=prov["version"], sha256=prov["sha256"]))
    return out


def station_day(station, date):
    """Station position (L3 lat/lon median of the day) and radiation/snow availability on that day."""
    f = RAW / f"{station}_hour.csv"
    if not f.exists() or date is None:
        return dict(radiation_status="no station file" if date else "no date")
    d = pd.read_csv(f, usecols=["time", "dsr", "dsr_cor", "usr_cor", "albedo", "snow_height", "lat", "lon"], parse_dates=["time"])
    day = d[d.time.dt.date == pd.Timestamp(date).date()]
    if day.empty:
        return dict(radiation_status="no station data on the day")
    return dict(gps_lat=float(day.lat.median()), gps_lon=float(day.lon.median()),
                hours_dsr_cor=int(day.dsr_cor.notna().sum()), hours_usr_cor=int(day.usr_cor.notna().sum()),
                hours_albedo=int(day.albedo.notna().sum()), snow_height_max=float(day.snow_height.max()),
                daily_albedo=float(day.usr_cor.sum() / day.dsr_cor.sum()) if day.dsr_cor.notna().any() else np.nan,
                radiation_status="available" if day.dsr_cor.notna().sum() >= 6 else "insufficient tilt-corrected SW")


def satellite_candidates(lat, lon, date):
    if lat is None or date is None:
        return dict(sat_status="no position or date")
    from pystac_client import Client
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "phase4"))
    import s2_io as io
    t = pd.Timestamp(date)
    s = Client.open(io.STAC_URL).search(collections=["sentinel-2-l2a"], intersects=dict(type="Point", coordinates=[lon, lat]),
                                        datetime=f"{(t - pd.Timedelta(days=SAT_WINDOW_DAYS)).date()}/{(t + pd.Timedelta(days=SAT_WINDOW_DAYS + 1)).date()}",
                                        max_items=100)
    items = [it.to_dict() for it in s.items()]
    if not items:
        return dict(sat_status="no L2A acquisition within ±3 d", n_candidates=0)
    rows = sorted(((abs((pd.Timestamp(i["properties"]["datetime"]).tz_localize(None) - t).total_seconds()) / 86400,
                    i["properties"].get("eo:cloud_cover", np.nan), i["id"]) for i in items))
    clear = [r for r in rows if r[1] < 50]
    best = clear[0] if clear else rows[0]
    return dict(sat_status="candidate" if clear else "only cloudy tiles (>= 50 %)", n_candidates=len(items),
                best_scene=best[2], best_scene_dt_days=round(best[0], 2), best_scene_tile_cloud=best[1],
                candidates=json.dumps([r[2] for r in rows]))


def qc(row):
    if (row["year"], row["station"]) in PROVIDER_QC_EXCLUDED:
        return "excluded: provider QC (TR314 §2, KAN samples not of good quality)"
    if row["date"] is None:
        return "excluded: no date"
    if row["obs_state"] in ("missing", "missing_sample", "unresolved_annotation"):
        return f"excluded: observation state {row['obs_state']}"
    return ""


def build(outdir=None, with_satellite=True):
    rows = load_2021_2023() + load_2024_rows()
    cache = {}
    for r in rows:
        key = (r["station"], r["date"])
        if key not in cache:
            sd = station_day(*key)
            lat = r["station_lat"] if r["station_lat"] is not None else sd.get("gps_lat")
            lon = r["station_lon"] if r["station_lon"] is not None else sd.get("gps_lon")
            sat = satellite_candidates(lat, lon, r["date"]) if with_satellite else {}
            cache[key] = dict(sd, **sat, match_lat=lat, match_lon=lon)
        r.update(cache[key])
        r["position_basis"] = "station position (workbook/GPS), NOT the sample position"
        r["spatial_support"] = {"under radiometer": "sample under the station radiometer (~cm² scrape, top 2 cm); radiometer footprint ~10-100 m²",
                                "randomized grid": "grid around the station; grid extent not documented",
                                "average ice in the area": "one sample chosen as representative (subjective)",
                                "dark ice": "targeted dark ice (biased high; not an area mean)"}.get(r["sampling_stratum"], "unknown")
        r["counting_protocol"] = "haemocytometer, 2 µl counted, LOD 1000 cells/mL, PFA-preserved"
        r["target"] = "ice algae (snow algae separate)"
        r["role"] = "final_test" if r["year"] == "2024" else "development"
        r["qc_exclusion"] = qc(r)
    df = pd.DataFrame(rows)
    if outdir:
        out = Path(outdir)
        out.mkdir(parents=True, exist_ok=True)
        df.to_csv(out / "manifest.csv", index=False)
        sd = (df.groupby(["role", "campaign", "station", "date"], dropna=False)
                .agg(n_obs=("obs_id", "size"), n_usable=("qc_exclusion", lambda x: int((x == "").sum())),
                     strata=("sampling_stratum", lambda x: ",".join(sorted(set(x)))),
                     radiation=("radiation_status", "first"), satellite=("sat_status", "first"),
                     sat_dt_days=("best_scene_dt_days", "first")).reset_index())
        sd.to_csv(out / "station_days.csv", index=False)
    return df


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--outdir", default="records/prombio_matching")
    p.add_argument("--no-satellite", action="store_true")
    a = p.parse_args()
    df = build(a.outdir, not a.no_satellite)
    pd.set_option("display.width", 250)
    print(pd.read_csv(Path(a.outdir) / "station_days.csv").to_string(index=False))
