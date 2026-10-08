"""
Sentinel-2 L2A ingestion: scene search, windowed band reads, scaling, masking,
DEM / slope on the same grid, GeoTIFF export.

Surface reflectance
-------------------
We ingest Level-2A (Sen2Cor bottom-of-atmosphere reflectance, R_sur) Cloud-Optimised
GeoTIFFs from the public AWS archive (Element 84 earth-search STAC, bucket
sentinel-cogs). Reflectance = DN * scale + offset; for processing baseline >= 04.00
(data from 25 Jan 2022) ESA introduced BOA_ADD_OFFSET = -1000 DN, i.e. offset -0.1.
We use the scale/offset recorded in the STAC item and enforce the -0.1 offset for
baselines >= 04.00 if the item does not carry it.

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
                GDAL_HTTP_MAX_RETRY="5", GDAL_HTTP_RETRY_DELAY="3")

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


def _scale_offset(item, asset):
    rb = (asset.get("raster:bands") or [{}])[0]
    scale, offset = rb.get("scale", 1e-4), rb.get("offset", 0.0)
    pb = float(item["properties"].get("s2:processing_baseline", "0") or 0)
    if pb >= 4.0 and offset == 0.0:
        offset = -0.1                                 # BOA_ADD_OFFSET = -1000 DN
    return scale, offset


def read_scene(item: dict, bounds, resolution: float = 20.0, keep_scl=(11,), retries: int = 4) -> Scene:
    """Windowed read with retries: the public COG bucket intermittently refuses existing files."""
    import time
    for k in range(retries):
        try:
            return _read_scene(item, bounds, resolution, keep_scl)
        except Exception as e:  # noqa: BLE001 - rasterio open/read errors
            if k == retries - 1:
                raise
            print(f"read_scene: {type(e).__name__}; retry {k + 1}/{retries - 1}", flush=True)
            time.sleep(10 * (k + 1))


def _read_scene(item: dict, bounds, resolution: float = 20.0, keep_scl=(11,)) -> Scene:
    import rasterio
    from rasterio.enums import Resampling
    from rasterio.windows import from_bounds

    out, transform, crs = {}, None, None
    with rasterio.Env(**GDAL_ENV):
        for b, key in BANDS.items():
            a = item["assets"][key]
            with rasterio.open(a["href"]) as src:
                w = from_bounds(*bounds, transform=src.transform).round_offsets().round_lengths()
                H = int(round((bounds[3] - bounds[1]) / resolution))
                W = int(round((bounds[2] - bounds[0]) / resolution))
                dn = src.read(1, window=w, out_shape=(H, W), resampling=Resampling.average).astype(np.float32)
                transform = rasterio.transform.from_bounds(*bounds, W, H)
                crs = src.crs
            scale, offset = _scale_offset(item, a)
            r = dn * scale + offset
            r[dn == 0] = np.nan                           # nodata
            out[b] = r
        with rasterio.open(item["assets"]["scl"]["href"]) as src:
            w = from_bounds(*bounds, transform=src.transform).round_offsets().round_lengths()
            scl = src.read(1, window=w, out_shape=out["B2"].shape, resampling=Resampling.nearest)
    mask = np.isin(scl, keep_scl) & np.all([np.isfinite(v) & (v > 0) for v in out.values()], axis=0)
    for b in out:
        out[b] = np.where(mask, out[b], np.nan).astype(np.float32)
    p = item["properties"]
    sza = 90.0 - float(p.get("view:sun_elevation", 45.0))
    return Scene(item, out, mask, transform, crs, sza, p.get("datetime", "")[:10])


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
