#!/usr/bin/env python
"""
Sentinel-2 station-pixel extraction for H2 and H4 (docs/glacier_model_validation_protocol.md §3, §5).

For each PROMICE station (KAN_L, KAN_M) and each Sentinel-2 L2A acquisition, June-August 2016-2023 (protocol amendment A1):
  * the station position on that day comes from the station's own hourly lat/lon (the stations move
    with the ice);
  * a 3 x 3 window of 20 m pixels is read, snapped to the tile grid, with the station in the centre
    pixel. Reflectance is area-averaged and SCL is taken by nearest;
  * station targets come from the hours starting at floor(t) - 1 h, floor(t) and floor(t) + 1 h
    (PROMICE hour-start convention), t = acquisition time:
      - absorbed SW = mean(dsr_cor - usr_cor);
      - albedo = sum(usr_cor) / sum(dsr_cor);
      - forcing for predictions: mean(dsr_cor).
    The station's SW-up and albedo are targets only; they are never inputs.

Population rules (protocol): all 9 SCL values = 11; snow_height < 0.02 m in all 3 hours; dsr_cor and
usr_cor finite in all 3 hours; one acquisition per station-day (the earliest that passes). Every
rejected acquisition is kept with its reason.

Sources are read remotely and never modified. Station data: data/promice_raw/<station>_hour.csv
(checksums in SHA256SUMS).

    python phase4/station_pixels.py --outdir phase4/results/station_pixels
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "phase2"))

RAW = os.path.join(HERE, "..", "data", "promice_raw")
STATIONS = ("KAN_L", "KAN_M")
YEARS = tuple(range(2016, 2024))            # protocol amendment A1
RES = 20.0
COLS = ["time", "dsr_cor", "usr_cor", "dsr", "usr", "snow_height", "lat", "lon", "tilt_x", "tilt_y", "cc"]


def station_hourly(station, years=YEARS):
    d = pd.read_csv(os.path.join(RAW, f"{station}_hour.csv"), usecols=COLS, parse_dates=["time"])
    return d[d.time.dt.year.isin(years) & d.time.dt.month.isin((6, 7, 8))].set_index("time").sort_index()


def search(lat, lon, year):
    from pystac_client import Client
    import s2_io as io
    cat = Client.open(io.STAC_URL)
    s = cat.search(collections=["sentinel-2-l2a"], intersects=dict(type="Point", coordinates=[lon, lat]),
                   datetime=f"{year}-06-01/{year}-08-31T23:59:59Z", max_items=1000)
    return sorted((it.to_dict() for it in s.items()), key=lambda d: d["properties"]["datetime"])


def station_targets(h, t):
    """Targets over the hours starting at floor(t) - 1 h ... floor(t) + 1 h; None and a reason if invalid."""
    t0 = pd.Timestamp(t).tz_localize(None).floor("h")
    hrs = [t0 + pd.Timedelta(hours=k) for k in (-1, 0, 1)]
    sub = h.reindex(hrs)
    if sub[["dsr_cor", "usr_cor"]].isna().any().any():
        return None, "station: dsr_cor/usr_cor missing in the 3 overpass hours"
    if (sub.snow_height.isna() | (sub.snow_height >= 0.02)).any():
        return None, "station: snow_height >= 0.02 m or missing"
    return dict(hours=",".join(x.strftime("%H") for x in hrs), sw_down=float(sub.dsr_cor.mean()),
                absorbed_obs=float((sub.dsr_cor - sub.usr_cor).mean()),
                albedo_obs=float(sub.usr_cor.sum() / sub.dsr_cor.sum()), cloud_cover=float(sub.cc.mean()),
                lat=float(sub.lat.mean()), lon=float(sub.lon.mean())), None


def window(item, lat, lon):
    """3 x 3 pixels at RES m snapped to the tile grid, the station in the centre pixel."""
    import s2_io as io
    from pyproj import Transformer
    from rasterio.crs import CRS
    p = item["properties"]
    crs = CRS.from_epsg(int(p.get("proj:epsg") or str(p.get("proj:code", "")).split(":")[-1]))
    x, y = Transformer.from_crs(4326, crs, always_xy=True).transform(lon, lat)
    tr = item["assets"]["blue"].get("proj:transform")
    x0, y0 = tr[2], tr[5]
    cx = x0 + (np.floor((x - x0) / RES) + 0.5) * RES                 # centre of the pixel containing x
    cy = y0 + (np.floor((y - y0) / RES) + 0.5) * RES
    b = (cx - 1.5 * RES, cy - 1.5 * RES, cx + 1.5 * RES, cy + 1.5 * RES)
    return io.read_scene(item, b, RES, keep_scl=tuple(range(1, 12)), dst_crs=crs, aoi_lonlat=(lon, lat)), (x, y)


def run(outdir, stations=STATIONS, years=YEARS, limit=None):
    rows = []
    for st in stations:
        h = station_hourly(st, years)
        for yr in years:
            hy = h[h.index.year == yr]
            if hy.empty:
                continue
            items = search(float(hy.lat.median()), float(hy.lon.median()), yr)
            print(f"[{st} {yr}] {len(items)} acquisitions", flush=True)
            done_days = set()
            for it in items[:limit]:
                t = pd.Timestamp(it["properties"]["datetime"])
                row = dict(station=st, year=yr, scene_id=it["id"], datetime=str(t), day=str(t.date()),
                           processing_baseline=it["properties"].get("s2:processing_baseline"),
                           tile_cloud=it["properties"].get("eo:cloud_cover"))
                tg, why = station_targets(h, t)
                if tg is None:
                    rows.append(dict(row, status="excluded", reason=why))
                    continue
                row.update(tg)
                if row["day"] in done_days:
                    rows.append(dict(row, status="excluded", reason="station-day already has an acquisition"))
                    continue
                try:
                    sc, xy = window(it, tg["lat"], tg["lon"])
                except Exception as e:  # noqa: BLE001 - recorded, not hidden
                    rows.append(dict(row, status="excluded", reason=f"read failed: {type(e).__name__}: {e}"[:200]))
                    continue
                R = sc.stack()
                row.update(sza=sc.sza, scl=json.dumps(sc.scl.ravel().tolist()),
                           **{f"{b}": float(R[1, 1, i]) for i, b in enumerate(("B2", "B3", "B4", "B8"))},
                           **{f"{b}_sd3x3": float(np.nanstd(R[..., i])) for i, b in enumerate(("B2", "B3", "B4", "B8"))},
                           scaling=json.dumps(sc.item["_read"]["scaling"]))
                if not np.all(sc.scl == 11):
                    rows.append(dict(row, status="excluded", reason="SCL != 11 within 30 m"))
                    continue
                if not np.all(np.isfinite(R)):
                    rows.append(dict(row, status="excluded", reason="nodata in window"))
                    continue
                done_days.add(row["day"])
                rows.append(dict(row, status="valid", reason=""))
    df = pd.DataFrame(rows)
    if outdir:
        import provenance as PV
        os.makedirs(outdir, exist_ok=True)
        PV.atomic_to_csv(df, os.path.join(outdir, "station_pixels.csv"), index=False, float_format="%.6g")
        PV.atomic_write_text(os.path.join(outdir, "settings.json"), json.dumps(dict(
            stations=list(stations), years=list(years), resolution_m=RES,
            raw_sha256=open(os.path.join(RAW, "SHA256SUMS")).read(),
            raw_downloaded_utc=open(os.path.join(RAW, "DOWNLOADED_UTC")).read().strip(),
            environment=PV.environment()), indent=1))
    return df


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--outdir", default=os.path.join(HERE, "results", "station_pixels"))
    p.add_argument("--stations", nargs="+", default=list(STATIONS))
    p.add_argument("--years", nargs="+", type=int, default=list(YEARS))
    p.add_argument("--limit", type=int, default=None, help="first N acquisitions per station-year (smoke test)")
    a = p.parse_args(argv)
    df = run(a.outdir, a.stations, a.years, a.limit)
    print(df.groupby(["station", "year", "status"]).size().to_string())
    print(df[df.status == "excluded"].reason.value_counts().to_string())


if __name__ == "__main__":
    main()
