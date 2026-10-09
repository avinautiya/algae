"""
Sentinel-2 L2A ingestion: scene search, windowed band reads, scaling, masking,
DEM / slope on the same grid, GeoTIFF export.

Surface reflectance
-------------------
We ingest Level-2A (Sen2Cor bottom-of-atmosphere reflectance, R_sur) Cloud-Optimised
GeoTIFFs from the public AWS archive (Element 84 earth-search STAC, bucket
sentinel-cogs). Reflectance = DN * scale + offset; for processing baseline >= 04.00
(data from 25 Jan 2022) ESA introduced BOA_ADD_OFFSET = -1000 DN, i.e. offset -0.1.
The offset rules are in `_scale_offset` (an explicit recorded value is never silently replaced;
contradictory metadata raises).

If you start from Level-1C (top of atmosphere), run ESA Sen2Cor
(https://step.esa.int/main/snap-supported-plugins/sen2cor/) or download L2A instead;
over bright ice, also consider that Sen2Cor's aerosol retrieval is weakly constrained
(no dark targets) - the inversion's multiplicative nuisance k and per-band sigma absorb
residual errors of a few percent.

The Scene Classification Layer (SCL) masks everything but class 11 (snow/ice) by
default (removes cloud, cloud shadow, water, rock).
"""

from __future__ import annotations

import json
import os
import urllib.request
from dataclasses import dataclass

import numpy as np

BANDS = {"B2": "blue", "B3": "green", "B4": "red", "B8": "nir"}
S3_BASE = "https://sentinel-cogs.s3.us-west-2.amazonaws.com/sentinel-s2-l2a-cogs"
STAC_URL = "https://earth-search.aws.element84.com/v1"
GDAL_ENV = dict(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR", CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif",
                GDAL_HTTP_MULTIRANGE="YES", GDAL_HTTP_MERGE_CONSECUTIVE_RANGES="YES",
                GDAL_HTTP_MAX_RETRY="8", GDAL_HTTP_RETRY_DELAY="1",
                # the bucket (through the session proxy) intermittently answers 404 for existing files
                GDAL_HTTP_RETRY_CODES="404,429,500,502,503,504")

# SW Greenland "Dark Zone" around PROMICE S6 (well-documented glacier-algal blooms)
DEFAULT_AOI = dict(lat=67.08, lon=-49.35, size_km=6.0)
DEFAULT_SCENE = "S2A_22WEV_20190723_0_L2A"      # 2.6 % cloud, full tile coverage


@dataclass
class Scene:
    item: dict
    bands: dict            # name -> (H, W) reflectance float32 (NaN = masked)
    mask: np.ndarray       # True = valid snow/ice pixel
    transform: object
    crs: object
    sza: float
    date: str
    scl: np.ndarray | None = None      # raw Sen2Cor scene classification on the same grid (0 = nodata)

    @property
    def shape(self):
        return self.mask.shape

    def stack(self):
        return np.stack([self.bands[b] for b in BANDS], axis=-1)


# --------------------------------------------------------------------------- #
# Scene discovery                                                               #
# --------------------------------------------------------------------------- #
def search_scenes(lat, lon, start, end, max_cloud=10.0, limit=20):
    """STAC search (pystac-client) for L2A scenes covering the point, least cloudy first."""
    from pystac_client import Client
    cat = Client.open(STAC_URL)
    s = cat.search(collections=["sentinel-2-l2a"], intersects=dict(type="Point", coordinates=[lon, lat]),
                   datetime=f"{start}/{end}", query={"eo:cloud_cover": {"lt": max_cloud}}, max_items=limit)
    items = [it.to_dict() for it in s.items()]
    items.sort(key=lambda d: (d["properties"].get("s2:nodata_pixel_percentage", 0) > 5,
                              d["properties"].get("eo:cloud_cover", 100)))
    return items


def item_from_s3(scene_id: str, retries: int = 3):
    """Fetch a STAC item JSON directly from the public bucket, e.g. 'S2A_22WEV_20190723_0_L2A'.
    The bucket intermittently answers 404 for existing items, so each source is retried, and the
    earth-search STAC API (/collections/sentinel-2-l2a/items/<id>) is the fallback."""
    import time
    _, tile, date, _, _ = scene_id.split("_")
    urls = [f"{S3_BASE}/{tile[:2]}/{tile[2]}/{tile[3:]}/{date[:4]}/{int(date[4:6])}/{scene_id}/{scene_id}.json",
            f"{STAC_URL}/collections/sentinel-2-l2a/items/{scene_id}"]
    last = None
    for url in urls:
        for k in range(retries):
            try:
                with urllib.request.urlopen(url, timeout=60) as r:
                    return json.loads(r.read())
            except Exception as e:  # noqa: BLE001 - HTTP 404/5xx or network
                last = e
                time.sleep(2 * (k + 1))
    raise RuntimeError(f"STAC item {scene_id} not available from S3 or the STAC API: {last}")


# --------------------------------------------------------------------------- #
# Reading                                                                       #
# --------------------------------------------------------------------------- #
def aoi_bounds(lat, lon, size_km, crs):
    from pyproj import Transformer
    x, y = Transformer.from_crs(4326, crs, always_xy=True).transform(lon, lat)
    h = size_km * 500.0
    return (x - h, y - h, x + h, y + h)


BOA_OFFSET_PB4 = -0.1                 # ESA BOA_ADD_OFFSET = -1000 DN for processing baseline >= 04.00


def _scale_offset(item, asset):
    """(scale, offset, source) for DN -> reflectance.

    Rules (an explicit value is never silently replaced):
      * processing baseline < 04.00: the recorded offset (expected 0).
      * baseline >= 04.00 and the item says the offset is applied in its metadata
        (earthsearch:boa_offset_applied true): the recorded offset, which must then be -0.1; a recorded
        0 contradicts the flag and raises (the DNs or the metadata are not what the reader assumes).
      * baseline >= 04.00 without that flag: Element 84 items made before the flag existed record
        offset 0 although the DNs carry ESA's -1000 DN offset, so -0.1 is applied once; a recorded
        nonzero offset is used as given.
      * baseline unknown (e.g. local files without metadata): raises - pass it explicitly."""
    rb = (asset.get("raster:bands") or [{}])[0]
    scale = float(rb.get("scale", 1e-4))
    rec = rb.get("offset")
    p = item["properties"]
    pb_raw = p.get("s2:processing_baseline")
    if pb_raw in (None, "", "unknown"):
        raise ValueError("processing baseline unknown: cannot decide the BOA offset")
    pb = float(pb_raw)
    flag = p.get("earthsearch:boa_offset_applied")
    if pb < 4.0:
        return scale, float(rec or 0.0), "recorded (baseline < 04.00)"
    if flag is True:
        if rec is None or not np.isclose(float(rec), BOA_OFFSET_PB4):
            raise ValueError(f"boa_offset_applied but recorded offset {rec!r} != {BOA_OFFSET_PB4}")
        return scale, float(rec), "recorded (boa_offset_applied)"
    if rec is None or float(rec) == 0.0:
        return scale, BOA_OFFSET_PB4, "ESA BOA_ADD_OFFSET applied (metadata predates the offset field)"
    return scale, float(rec), "recorded"


def read_scene(item: dict, bounds, resolution: float = 20.0, keep_scl=(11,), retries: int = 4,
               dst_crs=None, aoi_lonlat=None, max_sun_mismatch_deg: float = 1.0) -> Scene:
    """Windowed read with retries: the public COG bucket intermittently refuses existing files."""
    import time
    for k in range(retries):
        try:
            return _read_scene(item, bounds, resolution, keep_scl, dst_crs, aoi_lonlat, max_sun_mismatch_deg)
        except (ValueError, KeyError):
            raise                                             # metadata errors are not transient
        except Exception as e:  # noqa: BLE001 - rasterio open/read errors
            if k == retries - 1:
                raise
            print(f"read_scene: {type(e).__name__}; retry {k + 1}/{retries - 1}", flush=True)
            time.sleep(10 * (k + 1))


def solar_zenith_deg(lat: float, lon: float, when_utc) -> float:
    """Solar zenith angle (NOAA / Spencer-series approximation, about 0.1 deg; no refraction)."""
    import pandas as pd
    t = pd.Timestamp(when_utc)
    t = t.tz_convert("UTC") if t.tzinfo else t.tz_localize("UTC")
    g = 2 * np.pi / 365.0 * (t.dayofyear - 1 + (t.hour - 12 + t.minute / 60 + t.second / 3600) / 24)
    eqt = 229.18 * (0.000075 + 0.001868 * np.cos(g) - 0.032077 * np.sin(g) - 0.014615 * np.cos(2 * g)
                    - 0.040849 * np.sin(2 * g))
    dec = (0.006918 - 0.399912 * np.cos(g) + 0.070257 * np.sin(g) - 0.006758 * np.cos(2 * g)
           + 0.000907 * np.sin(2 * g) - 0.002697 * np.cos(3 * g) + 0.00148 * np.sin(3 * g))
    tst = t.hour * 60 + t.minute + t.second / 60 + eqt + 4 * lon
    ha = np.radians(tst / 4 - 180)
    la = np.radians(lat)
    cz = np.sin(la) * np.sin(dec) + np.cos(la) * np.cos(dec) * np.cos(ha)
    return float(np.degrees(np.arccos(np.clip(cz, -1, 1))))


def _dst_grid(bounds, resolution):
    """Output grid: bounds must be whole multiples of the resolution (no silent half-pixel shifts)."""
    import rasterio
    W = (bounds[2] - bounds[0]) / resolution
    H = (bounds[3] - bounds[1]) / resolution
    if not (np.isclose(W, round(W)) and np.isclose(H, round(H))):
        raise ValueError(f"bounds {bounds} are not a whole number of {resolution} m pixels")
    W, H = int(round(W)), int(round(H))
    return rasterio.transform.from_bounds(*bounds, W, H), H, W


def aligned_bounds(bounds, origin_xy, resolution):
    """Snap bounds outward onto the source pixel grid (origin_xy = tile corner, e.g. (499980, 7500000))
    so that native pixels are read without resampling shifts."""
    x0, y0 = origin_xy
    lo_x = x0 + np.floor((bounds[0] - x0) / resolution) * resolution
    hi_x = x0 + np.ceil((bounds[2] - x0) / resolution) * resolution
    lo_y = y0 + np.floor((bounds[1] - y0) / resolution) * resolution
    hi_y = y0 + np.ceil((bounds[3] - y0) / resolution) * resolution
    return (float(lo_x), float(lo_y), float(hi_x), float(hi_y))


def _read_scene(item: dict, bounds, resolution: float = 20.0, keep_scl=(11,), dst_crs=None, aoi_lonlat=None,
                max_sun_mismatch_deg: float = 1.0) -> Scene:
    """Bands are warped (GDAL WarpedVRT) onto ONE output grid defined by `bounds` in `dst_crs` (default:
    the item's CRS): area-average for reflectance (source nodata 0 excluded), nearest for SCL. The
    returned transform is that grid's, so every band and the mask share it exactly; a different
    dst_crs is reprojected rather than assumed."""
    import rasterio
    from rasterio.crs import CRS
    from rasterio.enums import Resampling
    from rasterio.vrt import WarpedVRT

    p = item["properties"]
    if dst_crs is None:
        code = p.get("proj:epsg") or str(p.get("proj:code", "")).split(":")[-1]
        dst_crs = CRS.from_epsg(int(code)) if code else None
    transform, H, W = _dst_grid(bounds, resolution)
    out, provenance = {}, {}

    def warp(href, resampling, dtype):
        with rasterio.open(href) as src:
            crs = dst_crs or src.crs
            with WarpedVRT(src, crs=crs, transform=transform, width=W, height=H, resampling=resampling,
                           src_nodata=0, nodata=0) as vrt:
                return vrt.read(1).astype(dtype), crs

    with rasterio.Env(**GDAL_ENV):
        for b, key in BANDS.items():
            a = item["assets"][key]
            dn, crs = warp(a["href"], Resampling.average, np.float32)
            scale, offset, src = _scale_offset(item, a)
            provenance[b] = dict(scale=scale, offset=offset, source=src)
            r = dn * scale + offset
            r[dn == 0] = np.nan                           # nodata
            out[b] = r
        scl, _ = warp(item["assets"]["scl"]["href"], Resampling.nearest, np.int16)
    mask = np.isin(scl, keep_scl) & np.all([np.isfinite(v) & (v > 0) for v in out.values()], axis=0)
    for b in out:
        out[b] = np.where(mask, out[b], np.nan).astype(np.float32)
    # illumination: the item's (tile-mean) sun elevation is required and checked against the solar
    # position at the AOI centre; the AOI value is used (the tile spans ~110 km)
    if "view:sun_elevation" not in p:
        raise ValueError("item has no view:sun_elevation; pass the solar zenith explicitly")
    sza_tile = 90.0 - float(p["view:sun_elevation"])
    sza = sza_tile
    sun = dict(sza_tile_mean=sza_tile)
    if aoi_lonlat is None and crs is not None:
        from pyproj import Transformer
        cx, cy = 0.5 * (bounds[0] + bounds[2]), 0.5 * (bounds[1] + bounds[3])
        aoi_lonlat = Transformer.from_crs(crs, 4326, always_xy=True).transform(cx, cy)
    dt = p.get("datetime", "")
    if aoi_lonlat is not None and dt and dt != "local":
        sza_aoi = solar_zenith_deg(aoi_lonlat[1], aoi_lonlat[0], dt)
        sun.update(sza_aoi=sza_aoi, aoi_lonlat=list(map(float, aoi_lonlat)), datetime=dt)
        if abs(sza_aoi - sza_tile) > max_sun_mismatch_deg:
            raise ValueError(f"solar zenith at the AOI ({sza_aoi:.2f}) differs from the item's tile mean "
                             f"({sza_tile:.2f}) by more than {max_sun_mismatch_deg} deg - check time/location")
        sza = sza_aoi
    item = dict(item)
    item["_read"] = dict(scaling=provenance, sun=sun, grid=dict(bounds=list(bounds), resolution=resolution,
                                                                crs=str(crs), shape=[H, W]),
                         valid_fraction=float(mask.mean()))
    return Scene(item, out, mask, transform, crs, sza, dt[:10], scl)


# --------------------------------------------------------------------------- #
# DEM (Copernicus GLO-30) and slope                                             #
# --------------------------------------------------------------------------- #
def _cop30_url(lat_floor, lon_floor):
    ns = "N" if lat_floor >= 0 else "S"
    ew = "E" if lon_floor >= 0 else "W"
    name = f"Copernicus_DSM_COG_10_{ns}{abs(lat_floor):02d}_00_{ew}{abs(lon_floor):03d}_00_DEM"
    return f"https://copernicus-dem-30m.s3.amazonaws.com/{name}/{name}.tif"


def read_dem(bounds, crs, shape, transform):
    """Copernicus GLO-30 elevation reprojected (bilinear) onto the scene grid, plus slope (deg)."""
    import rasterio
    from pyproj import Transformer
    from rasterio.warp import reproject, Resampling

    tr = Transformer.from_crs(crs, 4326, always_xy=True)
    xs = [bounds[0], bounds[2], bounds[0], bounds[2]]
    ys = [bounds[1], bounds[1], bounds[3], bounds[3]]
    lon, lat = tr.transform(xs, ys)
    dem = np.full(shape, np.nan, dtype=np.float32)
    with rasterio.Env(**GDAL_ENV):
        for la in range(int(np.floor(min(lat))), int(np.floor(max(lat))) + 1):
            for lo in range(int(np.floor(min(lon))), int(np.floor(max(lon))) + 1):
                tmp = np.full(shape, np.nan, dtype=np.float32)
                try:
                    with rasterio.open(_cop30_url(la, lo)) as src:
                        reproject(rasterio.band(src, 1), tmp, dst_transform=transform, dst_crs=crs,
                                  dst_nodata=np.nan, resampling=Resampling.bilinear)
                except Exception as e:                      # tile absent (ocean) or offline
                    print(f"DEM tile {la},{lo}: {e}")
                    continue
                dem = np.where(np.isfinite(tmp), tmp, dem)
    return dem, slope_deg(dem, abs(transform.a))


def slope_deg(dem, res_m):
    gy, gx = np.gradient(dem.astype(float), res_m)
    return np.degrees(np.arctan(np.hypot(gx, gy)))


def write_geotiff(path, arrays: dict, transform, crs, nodata=np.nan):
    """Multi-band float32 GeoTIFF with band descriptions."""
    import rasterio
    names = list(arrays)
    H, W = arrays[names[0]].shape
    with rasterio.open(path, "w", driver="GTiff", height=H, width=W, count=len(names), dtype="float32",
                       crs=crs, transform=transform, nodata=nodata, compress="deflate") as dst:
        for i, n in enumerate(names, 1):
            dst.write(arrays[n].astype(np.float32), i)
            dst.set_band_description(i, n)
