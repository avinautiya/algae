#!/usr/bin/env python
"""
Research-grade Sentinel-2 scene interpretation product (Phase 5 items 4, 5 and 13).

Question answered per pixel: is the observed darkening consistent with glacier algae, and how
large is the modelled algal albedo and radiative effect? Dark ice is NOT classified as algae
because it is dark. The interpretation classes are (decision order):

  0 nodata
  1 cloud_or_shadow          SCL 1 (saturated), 3, 8, 9, 10
  2 water_or_pond            SCL 6, or NDWI_ice = (B2-B4)/(B2+B4) > 0.25 (Yang & Smith 2013)
  3 snow_or_bright_surface   B3 > 0.80 and B8 > 0.70: the bare-ice model does not apply
  9 non_ice_surface          SCL 4 (vegetation), 5 (bare soil)
  8 forward_model_inconsistent   posterior predictive p < 0.01
  7 prior_dominated          posterior SD(log B) > 0.8 x prior SD, and mean within 0.5 prior SD
                             of the prior mean
  5 biologically_consistent_darkening_supported
                             P(dalpha_algae >= 0.01) >= 0.9 AND Bayes factor (full model vs
                             no-algae) >= 10
  4 bare_ice_no_meaningful_algal_darkening   P(dalpha_algae >= 0.01) <= 0.1
  6 ambiguous_algae_vs_dust_or_ice           everything else

  SCL 2 (dark area) and 7 (unclassified) are NOT excluded: Sen2Cor often labels dark bare ice
  that way. They are interpreted like ice and flagged (bit 4).

Thresholds (fixed before any scene output was examined):
  * 0.01 albedo: the minimum meaningful improvement of the frozen protocol
    (docs/glacier_model_validation_protocol.md);
  * Bayes factor 10: "strong" on the Jeffreys scale;
  * 0.9 / 0.1: posterior probabilities.

The physics:
  * Frozen held-out emulators (records/heldout_v2), posterior-mean optics (deviation D1).
  * Primary-fold field calibration, trained on S6 2017: sigma, radius prior, tau.
  * Exact grid posterior (inversion.GridPosterior), run three times:
      - full prior;
      - no algae (log B at its lowest node, 10 cells/mL);
      - Tier A optics, as a structural alternative.
  * Species fractions are NOT reported: f_n is prior-dominated, and the shrinkage is in the
    report.

Uncertainty:
  * posterior SD;
  * structural tau (dex), added in quadrature to the abundance interval;
  * the optics alternative (Tier A minus TD-DFT tier D);
  * calibration draws: a single shared draw id "posterior_mean" (D1). Joint draws are a recorded
    follow-up.

Source imagery is read remotely and never modified. All outputs share one explicit grid:
UTM of the tile, `resolution` m, bounds snapped to the tile's pixel grid, area-average
resampling for reflectance and nearest for SCL.

    python common/run_budgeted.py --name scene_product --mem-mb 2500 --threads 1 -- \
        python3 phase4/scene_product.py --scene-id S2A_22WEV_20190723_0_L2A --size-km 3 \
        --outdir phase4/results/scene_product
"""

from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "phase2"))

SCHEMA_VERSION = "1.0.0"
SCHEMA_ID = f"algae-scene-product/{SCHEMA_VERSION}"
EMU_DIR = os.path.join(HERE, "results", "heldout_v2", "cache")
SETTINGS = os.path.join(HERE, "results", "heldout_v2", "heldout_settings.json")
PRIMARY, ALTERNATIVE = "tddft_D", "tierA_empirical"
CAL_FOLD = "primary"                     # trained on S6 2017
DALPHA_MIN = 0.01
LOG_BF_MIN = float(np.log(10.0))
P_SUPPORT, P_NONE = 0.9, 0.1
PRIOR_DOMINATED_SD_RATIO, PRIOR_DOMINATED_SHIFT = 0.8, 0.5
PPP_MIN = 0.01
NDWI_ICE_WATER = 0.25
SNOW_B3, SNOW_B8 = 0.80, 0.70
DRAW_ID = "posterior_mean"

CLASSES = {0: "nodata", 1: "cloud_or_shadow", 2: "water_or_pond", 3: "snow_or_bright_surface",
           4: "bare_ice_no_meaningful_algal_darkening", 5: "biologically_consistent_darkening_supported",
           6: "ambiguous_algae_vs_dust_or_ice", 7: "prior_dominated", 8: "forward_model_inconsistent",
           9: "non_ice_surface"}
FLAGS = {0: "nodata", 1: "scl_cloud_or_shadow", 2: "water_test", 3: "snow_test", 4: "scl_dark_or_unclassified",
         5: "ppp_below_0.01", 6: "prior_dominated", 7: "algae_dust_ambiguous", 8: "posterior_at_grid_edge",
         9: "scl_non_ice"}
SCL_CLOUD = (1, 3, 8, 9, 10)
SCL_WATER = (6,)
SCL_NON_ICE = (4, 5)
SCL_DARK = (2, 7)

# variable -> (units, long name, dtype, valid range)
VARIABLES = {
    "interpretation_class": ("1", "interpretation class (see flag_meanings)", "uint8", (0, 9)),
    "quality_flags": ("1", "quality bit flags (see flag_meanings)", "uint16", (0, 1023)),
    "log10_cells_mean": ("log10(cells mL-1 meltwater)", "posterior mean algal abundance", "float32", (0.5, 6.5)),
    "log10_cells_sd_posterior": ("log10(cells mL-1 meltwater)", "posterior SD of abundance", "float32", (0, 3)),
    "log10_cells_q025_total": ("log10(cells mL-1 meltwater)", "2.5 % (posterior SD (+) tau)", "float32", (-5, 8)),
    "log10_cells_q975_total": ("log10(cells mL-1 meltwater)", "97.5 % (posterior SD (+) tau)", "float32", (-5, 12)),
    "dalpha_algae_mean": ("1", "algal broadband albedo reduction (signed; >0 darkening)", "float32", (-0.2, 0.6)),
    "dalpha_algae_sd": ("1", "posterior SD of the algal albedo reduction", "float32", (0, 0.6)),
    "p_dalpha_ge_min": ("1", "posterior probability that the algal albedo reduction >= 0.01", "float32", (0, 1)),
    "log_bayes_factor_algae": ("1", "ln evidence(full) - ln evidence(no algae)", "float32", (-1e4, 1e4)),
    "rf_algae_mean": ("W m-2", "instantaneous clear-sky algal SW forcing at overpass SZA (modelled SW)", "float32", (-200, 600)),
    "bba_mean": ("1", "modelled broadband (300-2500 nm) albedo", "float32", (0, 1)),
    "bba_sd": ("1", "posterior SD of modelled broadband albedo", "float32", (0, 1)),
    "ppp": ("1", "posterior predictive p-value (k-marginal, 4 bands)", "float32", (0, 1)),
    "prior_sd_ratio_log10_cells": ("1", "posterior SD / prior SD of abundance", "float32", (0, 5)),
    "dalpha_algae_alt_optics": ("1", "algal albedo reduction with Tier A optics (structural alternative)", "float32", (-0.2, 0.6)),
}


# --------------------------------------------------------------------------- emulators and calibration
def _sha(path, n=None):
    h = hashlib.sha256(open(path, "rb").read()).hexdigest()
    return h[:n] if n else h


def load_physics(sza):
    """Frozen emulators at the scene's SZA node and the primary-fold field calibration.
    Refuses (no interpolation, no inferred illumination) when the node does not exist."""
    import emulator as E
    node = int(round(sza))
    out = {}
    for opt in (PRIMARY, ALTERNATIVE):
        p = os.path.join(EMU_DIR, f"heldout_{opt}_sza{node}.npz")
        if not os.path.isfile(p):
            avail = sorted(int(f.split("sza")[1][:-4]) for f in os.listdir(EMU_DIR) if f.startswith(f"heldout_{opt}_sza"))
            raise ValueError(f"scene SZA {sza:.2f} deg -> node {node} has no frozen emulator for {opt} "
                             f"(available {avail}); out of domain - refused")
        out[opt] = dict(em=E.Emulator.load(p).refine_log_b(0.05), path=p, sha256=_sha(p))
    s = json.load(open(SETTINGS))["settings"]
    for opt in out:
        out[opt]["cal"] = s[f"{CAL_FOLD}/{opt}"]
    return out


def add_indicator(em):
    """Derived node quantities: signed algal albedo reduction and the indicator dalpha >= DALPHA_MIN."""
    d = em.data["rf_algae"] / float(em.meta["sw_down"])
    em.data["dalpha_algae"] = d
    em.data["ind_dalpha"] = (d >= DALPHA_MIN).astype(float)
    return em


def prior_config(cal):
    import dataclasses
    from priors import PriorConfig
    return dataclasses.replace(PriorConfig.for_density(690.0), mu_lnr=float(np.log(cal["radius_median_um"])),
                               sd_lnr=float(cal["radius_ln_sd"]))


def posterior(em, R, cal, restrict_no_algae=False):
    import inversion as INV
    from priors import prior_logpdfs
    pc = prior_config(cal)
    lp, sk, mk, _ = prior_logpdfs(em.axes, R.shape[0], pc)
    if restrict_no_algae:
        lb = np.full_like(lp["log_b"], -np.inf)
        lb[:, 0] = 0.0
        lp = dict(lp, log_b=lb)
    derived = ("bba", "rf_algae", "dalpha_algae", "ind_dalpha")
    gp = INV.GridPosterior(em, np.full(4, float(cal["sigma"])), derived=derived)
    r = gp.run(R, lp, sk, mk=mk, marginals=("log_b", "r_um", "dust_ppb") if not restrict_no_algae else ())
    return r, pc


# --------------------------------------------------------------------------- interpretation
def screen(R, scl):
    """Exclusion classes and flags from SCL and reflectance tests (no physics)."""
    P = R.shape[0]
    cls = np.full(P, 255, np.uint8)                       # 255 = to be interpreted
    flags = np.zeros(P, np.uint16)
    finite = np.all(np.isfinite(R), axis=1) & (scl != 0)
    b2, b3, b4, b8 = (R[:, i] for i in range(4))
    with np.errstate(invalid="ignore", divide="ignore"):
        ndwi = (b2 - b4) / (b2 + b4)
    cloud = np.isin(scl, SCL_CLOUD)
    water = np.isin(scl, SCL_WATER) | (ndwi > NDWI_ICE_WATER)
    snow = (b3 > SNOW_B3) & (b8 > SNOW_B8)
    nonice = np.isin(scl, SCL_NON_ICE)
    dark = np.isin(scl, SCL_DARK)
    flags |= (~finite).astype(np.uint16) << 0
    flags |= cloud.astype(np.uint16) << 1
    flags |= (water & finite).astype(np.uint16) << 2
    flags |= (snow & finite).astype(np.uint16) << 3
    flags |= dark.astype(np.uint16) << 4
    flags |= nonice.astype(np.uint16) << 9
    for c, m in ((0, ~finite), (1, cloud), (2, water), (3, snow), (9, nonice)):
        cls[(cls == 255) & m] = c
    return cls, flags, ndwi


def interpret(R, scl, phys):
    """Per-pixel interpretation. R: (P, 4) BOA reflectance B2, B3, B4, B8; scl: (P,) SCL codes."""
    cls, flags, ndwi = screen(R, scl)
    todo = np.flatnonzero(cls == 255)
    P = R.shape[0]
    out = {k: np.full(P, np.nan, np.float32) for k in VARIABLES if k not in ("interpretation_class", "quality_flags")}
    diag = {}
    if todo.size:
        em = add_indicator(phys[PRIMARY]["em"])
        cal = phys[PRIMARY]["cal"]
        full, pc = posterior(em, R[todo], cal)
        none, _ = posterior(em, R[todo], cal, restrict_no_algae=True)
        ema = add_indicator(phys[ALTERNATIVE]["em"])
        alt, _ = posterior(ema, R[todo], phys[ALTERNATIVE]["cal"])
        tau = float(cal["tau_dex"])
        m, sd = full["log_b_mean"], full["log_b_sd"]
        tot = np.hypot(sd, tau)
        lbf = full["log_evidence"] - none["log_evidence"]
        p = full["ind_dalpha_mean"]
        ratio = sd / pc.sd_b
        prior_dom = (ratio > PRIOR_DOMINATED_SD_RATIO) & (np.abs(m - pc.mu_b) < PRIOR_DOMINATED_SHIFT * pc.sd_b)
        edge = np.zeros(todo.size, bool)
        for ax in ("log_b", "r_um", "dust_ppb"):
            mg = full[f"marg_{ax}"]
            edge |= (mg[:, 0] > 0.05) & (ax != "log_b")          # low log B is a physical state, not an edge
            edge |= mg[:, -1] > 0.05
        c = np.full(todo.size, 6, np.uint8)
        c[p <= P_NONE] = 4
        c[(p >= P_SUPPORT) & (lbf >= LOG_BF_MIN)] = 5
        c[prior_dom] = 7
        c[full["ppp"] < PPP_MIN] = 8
        cls[todo] = c
        f = flags[todo]
        f |= (full["ppp"] < PPP_MIN).astype(np.uint16) << 5
        f |= prior_dom.astype(np.uint16) << 6
        f |= (c == 6).astype(np.uint16) << 7
        f |= edge.astype(np.uint16) << 8
        flags[todo] = f
        vals = dict(log10_cells_mean=m, log10_cells_sd_posterior=sd, log10_cells_q025_total=m - 1.96 * tot,
                    log10_cells_q975_total=m + 1.96 * tot, dalpha_algae_mean=full["dalpha_algae_mean"],
                    dalpha_algae_sd=full["dalpha_algae_sd"], p_dalpha_ge_min=p, log_bayes_factor_algae=lbf,
                    rf_algae_mean=full["rf_algae_mean"], bba_mean=full["bba_mean"], bba_sd=full["bba_sd"],
                    ppp=full["ppp"], prior_sd_ratio_log10_cells=ratio, dalpha_algae_alt_optics=alt["dalpha_algae_mean"])
        for k, v in vals.items():
            out[k][todo] = v
        diag = dict(tau_dex=tau, prior_mu_log_b=float(pc.mu_b), prior_sd_log_b=float(pc.sd_b),
                    f_n_posterior_prior_sd_ratio_median=float(np.nanmedian(full["f_n_sd"]) / _beta_sd(pc)),
                    sw_down_model=float(em.meta["sw_down"]))
    out["interpretation_class"] = cls
    out["quality_flags"] = flags
    return out, diag


def _beta_sd(pc):
    a, b = pc.f_alpha, pc.f_beta
    return float(np.sqrt(a * b / ((a + b) ** 2 * (a + b + 1))))


# --------------------------------------------------------------------------- export
def _grid_coords(transform, H, W):
    x = transform.c + transform.a * (np.arange(W) + 0.5)
    y = transform.f + transform.e * (np.arange(H) + 0.5)
    return x, y


def export(maps, transform, crs, outdir, stem, provenance):
    """COG (continuous float32 + class/flag uint16), CF NetCDF4, provenance sidecar. Returns paths."""
    import rasterio
    os.makedirs(outdir, exist_ok=True)
    H, W = maps["interpretation_class"].shape
    cont = [k for k in VARIABLES if VARIABLES[k][2] == "float32"]
    paths = dict(cog_continuous=os.path.join(outdir, f"{stem}_continuous.tif"),
                 cog_class=os.path.join(outdir, f"{stem}_class_flags.tif"),
                 netcdf=os.path.join(outdir, f"{stem}.nc"), sidecar=os.path.join(outdir, f"{stem}_provenance.json"))
    for p in paths.values():
        if os.path.exists(p):
            raise FileExistsError(f"{p} exists: outputs are immutable, choose a new outdir")
    common = dict(driver="COG", height=H, width=W, crs=crs, transform=transform, compress="DEFLATE",
                  BLOCKSIZE=256, OVERVIEWS="NONE")
    tmp = paths["cog_continuous"] + ".tmp.tif"
    with rasterio.open(tmp, "w", count=len(cont), dtype="float32", nodata=np.nan, **common) as dst:
        for i, k in enumerate(cont, 1):
            dst.write(maps[k].astype(np.float32), i)
            dst.set_band_description(i, k)
            dst.update_tags(i, units=VARIABLES[k][0], long_name=VARIABLES[k][1])
        dst.update_tags(schema=SCHEMA_ID, calibration_draw_id=DRAW_ID)
    os.replace(tmp, paths["cog_continuous"])
    tmp = paths["cog_class"] + ".tmp.tif"
    with rasterio.open(tmp, "w", count=2, dtype="uint16", nodata=None, **common) as dst:
        dst.write(maps["interpretation_class"].astype(np.uint16), 1)
        dst.write(maps["quality_flags"].astype(np.uint16), 2)
        dst.set_band_description(1, "interpretation_class")
        dst.set_band_description(2, "quality_flags")
        dst.update_tags(1, flag_values=json.dumps(list(CLASSES)), flag_meanings=json.dumps(CLASSES))
        dst.update_tags(2, flag_masks=json.dumps([1 << b for b in FLAGS]), flag_meanings=json.dumps(FLAGS))
        dst.update_tags(schema=SCHEMA_ID)
    os.replace(tmp, paths["cog_class"])
    _write_netcdf(paths["netcdf"], maps, transform, crs, provenance)
    prov = dict(provenance, schema=SCHEMA_ID, calibration_draw_ids=[DRAW_ID],
                grid=dict(crs_wkt=crs.to_wkt(), transform=list(transform)[:6], shape=[H, W]),
                outputs={k: dict(path=os.path.basename(v), sha256=_sha(v)) for k, v in paths.items() if k != "sidecar"})
    tmp = paths["sidecar"] + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(prov, fh, indent=1, allow_nan=False, default=str)
    os.replace(tmp, paths["sidecar"])
    return paths


def _write_netcdf(path, maps, transform, crs, provenance):
    import netCDF4
    H, W = maps["interpretation_class"].shape
    x, y = _grid_coords(transform, H, W)
    tmp = path + ".tmp"
    with netCDF4.Dataset(tmp, "w", format="NETCDF4") as nc:
        nc.Conventions = "CF-1.8"
        nc.title = "Glacier algae scene interpretation (research-grade)"
        nc.schema = SCHEMA_ID
        nc.history = f"{_dt.datetime.now(_dt.timezone.utc).isoformat()} scene_product.py"
        nc.source_scene = str(provenance.get("scene_id"))
        nc.validation_status = provenance.get("validation_status", "")
        nc.createDimension("y", H)
        nc.createDimension("x", W)
        nc.createDimension("draw", 1)
        vx = nc.createVariable("x", "f8", ("x",))
        vx[:] = x
        vx.units, vx.standard_name = "m", "projection_x_coordinate"
        vy = nc.createVariable("y", "f8", ("y",))
        vy[:] = y
        vy.units, vy.standard_name = "m", "projection_y_coordinate"
        vd = nc.createVariable("calibration_draw_id", str, ("draw",))
        vd[0] = DRAW_ID
        g = nc.createVariable("crs", "i4")
        g.crs_wkt = crs.to_wkt()
        g.spatial_ref = crs.to_wkt()
        g.GeoTransform = " ".join(str(v) for v in (transform.c, transform.a, transform.b, transform.f, transform.d, transform.e))
        for k, (units, ln, dtype, rng) in VARIABLES.items():
            if dtype == "float32":
                v = nc.createVariable(k, "f4", ("y", "x"), zlib=True, fill_value=np.float32(np.nan))
                v[:] = maps[k].astype(np.float32)
            else:
                v = nc.createVariable(k, "u2" if dtype == "uint16" else "u1", ("y", "x"), zlib=True)
                v[:] = maps[k].astype(np.uint16 if dtype == "uint16" else np.uint8)
            v.units, v.long_name, v.grid_mapping = units, ln, "crs"
            cast = {"float32": np.float32, "uint16": np.uint16, "uint8": np.uint8}[dtype]
            v.valid_min, v.valid_max = cast(rng[0]), cast(rng[1])
            if k == "interpretation_class":
                v.flag_values = np.array(list(CLASSES), np.uint8)
                v.flag_meanings = " ".join(CLASSES.values())
            if k == "quality_flags":
                v.flag_masks = np.array([1 << b for b in FLAGS], np.uint16)
                v.flag_meanings = " ".join(FLAGS.values())
    os.replace(tmp, path)


def validate_outputs(paths, maps, transform, crs):
    """Round trip: values (NaN-aware), dtypes, CRS/transform, coordinates, masks, units, checksums,
    COG layout. Returns a dict of checks; raises AssertionError on the first failure."""
    import netCDF4
    import rasterio
    checks = {}
    side = json.load(open(paths["sidecar"]))
    for k, meta in side["outputs"].items():
        assert _sha(paths[k]) == meta["sha256"], f"checksum {k}"
    checks["checksums"] = True
    cont = [k for k in VARIABLES if VARIABLES[k][2] == "float32"]
    with rasterio.open(paths["cog_continuous"]) as src:
        assert src.tags(ns="IMAGE_STRUCTURE").get("LAYOUT") == "COG", "not a COG layout"
        assert src.crs == crs and src.transform.almost_equals(transform), "COG georeference"
        assert list(src.descriptions) == cont, "band order"
        for i, k in enumerate(cont, 1):
            np.testing.assert_array_equal(src.read(i), maps[k].astype(np.float32), err_msg=k)
            assert src.tags(i)["units"] == VARIABLES[k][0]
    with rasterio.open(paths["cog_class"]) as src:
        assert src.tags(ns="IMAGE_STRUCTURE").get("LAYOUT") == "COG"
        np.testing.assert_array_equal(src.read(1), maps["interpretation_class"])
        np.testing.assert_array_equal(src.read(2), maps["quality_flags"])
    checks["cog_roundtrip"] = True
    H, W = maps["interpretation_class"].shape
    x, y = _grid_coords(transform, H, W)
    with netCDF4.Dataset(paths["netcdf"]) as nc:
        np.testing.assert_allclose(nc["x"][:], x)
        np.testing.assert_allclose(nc["y"][:], y)
        for k, (units, *_r) in VARIABLES.items():
            v = nc[k]
            assert v.units == units, k
            a = np.ma.filled(v[:], np.nan) if VARIABLES[k][2] == "float32" else np.asarray(v[:])
            np.testing.assert_array_equal(a, maps[k].astype(a.dtype), err_msg=k)
        assert nc["calibration_draw_id"][0] == DRAW_ID
    checks["netcdf_roundtrip"] = True
    # masks: excluded pixels carry no physical values; interpreted pixels have them
    cls = maps["interpretation_class"]
    excluded = np.isin(cls, (0, 1, 2, 3, 9))
    for k in cont:
        assert np.all(np.isnan(maps[k][excluded])), f"{k} has values on excluded pixels"
    assert np.all(np.isfinite(maps["log10_cells_mean"][~excluded])), "interpreted pixel without abundance"
    for k, (_u, _l, dtype, (lo, hi)) in VARIABLES.items():
        v = maps[k][np.isfinite(maps[k])] if dtype == "float32" else maps[k]
        assert v.size == 0 or (v.min() >= lo and v.max() <= hi), f"{k} outside its valid range"
    checks["masks_and_ranges"] = True
    return checks


# --------------------------------------------------------------------------- report
def report(maps, diag, provenance, res_m):
    cls = maps["interpretation_class"]
    px_km2 = (res_m / 1000.0) ** 2
    n = cls.size
    lines = [f"# Scene interpretation report: {provenance.get('scene_id')}", "",
             f"Schema `{SCHEMA_ID}`. Research-grade product: NOT operationally validated. See the "
             "validation status at the end.", "",
             "## Pixel classes", "", "| class | pixels | % | km² |", "|---|---|---|---|"]
    for c, name in CLASSES.items():
        k = int((cls == c).sum())
        lines.append(f"| {c} {name} | {k} | {100 * k / n:.1f} | {k * px_km2:.2f} |")
    sup = cls == 5
    amb = cls == 6
    interp = np.isin(cls, (4, 5, 6, 7, 8))

    def q(v, m):
        x = v[m & np.isfinite(v)]
        return "n/a" if x.size == 0 else f"{np.median(x):.3f} (IQR {np.quantile(x, .25):.3f}–{np.quantile(x, .75):.3f})"
    lines += ["", "## 1. Where is biologically consistent darkening supported?",
              f"- Class 5 covers {sup.sum()} pixels ({100 * sup.sum() / max(interp.sum(), 1):.1f} % of the interpreted pixels).",
              f"- Criterion: P(Δα ≥ {DALPHA_MIN}) ≥ {P_SUPPORT} and Bayes factor ≥ 10 against a no-algae explanation.",
              "## 2. Where are algae vs dust/ice/water explanations ambiguous?",
              f"- Class 6 covers {amb.sum()} pixels.",
              f"- Class 7 (prior-dominated) covers {(cls == 7).sum()} pixels.",
              f"- Class 8 (model inconsistent) covers {(cls == 8).sum()} pixels.",
              f"- Pixels with dark/unclassified SCL that were still interpreted: {int(((maps['quality_flags'] >> 4) & 1)[interp].sum())}.",
              "## 3. How large is the modelled radiative effect?",
              f"- Δα_algae, median (IQR), class 5: {q(maps['dalpha_algae_mean'], sup)}; all interpreted pixels: {q(maps['dalpha_algae_mean'], interp)}.",
              f"- Instantaneous clear-sky forcing at the overpass SZA (modelled SW {diag.get('sw_down_model', float('nan')):.0f} W m⁻²), class 5: {q(maps['rf_algae_mean'], sup)} W m⁻².",
              "- This is not a daily or seasonal forcing, and not melt.",
              "## 4. What uncertainty dominates?"]
    if interp.any():
        psd = np.nanmedian(maps["dalpha_algae_sd"][interp])
        struct = np.nanmedian(np.abs(maps["dalpha_algae_alt_optics"] - maps["dalpha_algae_mean"])[interp])
        abd = np.nanmedian(maps["log10_cells_sd_posterior"][interp])
        lines += [f"- Δα: posterior SD (median) {psd:.4f} vs optics structural difference |Tier A − TD-DFT D| (median) {struct:.4f}.",
                  f"- Abundance: posterior SD (median) {abd:.2f} dex vs structural τ {diag['tau_dex']:.2f} dex.",
                  f"- **Dominant:** {'optics structure' if struct > psd else 'retrieval (posterior) spread'} for Δα; "
                  f"{'structural τ' if diag['tau_dex'] > abd else 'posterior spread'} for abundance.",
                  "- Calibration uncertainty is not sampled: single draw `posterior_mean`, deviation D1.",
                  f"- Species fraction: median posterior/prior SD of f_n = {diag['f_n_posterior_prior_sd_ratio_median']:.2f}, so it is not reported."]
    lines += ["## 5. Which pixels should be excluded?",
              "- Classes 0–3 and 9 (no physical values are written).",
              "- Classes 7 and 8 should not be used for abundance or forcing.",
              "- Quality bit 8 marks posterior mass at a grid edge (radius or dust): extrapolation risk.",
              "## 6. What is actually independently validated?",
              "- Abundance retrieval: held-out at plot scale (`records/heldout_v2`), weak evidence (P5-CORR-1).",
              "- Pixel-scale abundance, albedo and forcing: **not independently validated**.",
              "- H1–H4 status: see `docs/evidence_table.md`.",
              "- Thresholds, priors and calibration fold: `docs/scene_product.md`."]
    return "\n".join(lines) + "\n"


def schema():
    """Machine-readable output contract (docs/schema/scene_product_v1.json is generated from this)."""
    return dict(schema=SCHEMA_ID, grid="projected UTM of the source tile; cell-centre x/y coordinates; north-up",
                nodata=dict(float32="NaN", class_value=0),
                variables={k: dict(units=u, long_name=ln, dtype=dt, valid_range=list(r)) for k, (u, ln, dt, r) in VARIABLES.items()},
                classes={str(k): v for k, v in CLASSES.items()}, flags={str(1 << b): v for b, v in FLAGS.items()},
                validity_rules=[
                    "physical variables are NaN on classes 0-3 and 9",
                    "classes 7 and 8 carry values but must not be used for abundance or forcing",
                    "quality bit 8: posterior mass at a radius/dust grid edge (extrapolation risk)",
                    "dalpha and forcing are SIGNED (negative = brightening); never clipped",
                    "rf_algae_mean is instantaneous clear-sky at the overpass SZA with modelled SW, not daily/seasonal",
                    "no species fractions are reported (f_n prior-dominated)"],
                thresholds=dict(dalpha_min=DALPHA_MIN, bayes_factor_min=10.0, p_support=P_SUPPORT, p_none=P_NONE,
                                ppp_min=PPP_MIN, ndwi_ice_water=NDWI_ICE_WATER, snow_b3=SNOW_B3, snow_b8=SNOW_B8,
                                prior_dominated_sd_ratio=PRIOR_DOMINATED_SD_RATIO,
                                prior_dominated_shift_prior_sd=PRIOR_DOMINATED_SHIFT),
                uncertainty=dict(calibration_draw_ids=[DRAW_ID], structural_tau="added in quadrature to abundance interval",
                                 optics_alternative="dalpha_algae_alt_optics (Tier A)"),
                files=["<scene>_continuous.tif (COG float32)", "<scene>_class_flags.tif (COG uint16)",
                       "<scene>.nc (NetCDF4 CF-1.8)", "<scene>_provenance.json", "<scene>_report.md"])


# --------------------------------------------------------------------------- driver
def code_provenance():
    import subprocess
    mods = ["phase4/scene_product.py", "phase4/inversion.py", "phase4/s2_io.py", "phase4/priors.py", "phase4/emulator.py"]
    root = os.path.join(HERE, "..")
    try:
        commit = subprocess.run(["git", "-C", root, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
        dirty = bool(subprocess.run(["git", "-C", root, "status", "--porcelain", "--"] + mods,
                                    capture_output=True, text=True).stdout.strip())
    except OSError:
        commit, dirty = "unknown", True
    return dict(commit=commit, modules_dirty=dirty, module_sha256={m: _sha(os.path.join(root, m)) for m in mods})


def run(scene_id, lat, lon, size_km, resolution, outdir):
    import s2_io as io
    item = io.item_from_s3(scene_id)
    item_sha = hashlib.sha256(json.dumps(item, sort_keys=True).encode()).hexdigest()
    p = item["properties"]
    code = p.get("proj:epsg") or str(p.get("proj:code", "")).split(":")[-1]
    from rasterio.crs import CRS
    crs = CRS.from_epsg(int(code))
    b = io.aoi_bounds(lat, lon, size_km, crs)
    origin = item["assets"]["blue"].get("proj:transform")
    b = io.aligned_bounds(b, (origin[2], origin[5]), resolution) if origin else b
    scene = io.read_scene(item, b, resolution, keep_scl=tuple(range(1, 12)), dst_crs=crs, aoi_lonlat=(lon, lat))
    phys = load_physics(scene.sza)
    H, W = scene.shape
    R = scene.stack().reshape(-1, 4).astype(float)
    scl = scene.scl.reshape(-1)
    vals, diag = interpret(R, scl, phys)
    maps = {k: v.reshape(H, W) for k, v in vals.items()}
    prov = dict(scene_id=scene_id, item_sha256=item_sha, processing_baseline=p.get("s2:processing_baseline"),
                datetime=p.get("datetime"), assets={k: item["assets"][k]["href"] for k in ("blue", "green", "red", "nir", "scl")},
                read=scene.item["_read"], sza_used=scene.sza, sza_source="solar position at AOI centre (checked against item)",
                emulators={k: dict(path=os.path.relpath(v["path"], os.path.join(HERE, "..")), sha256=v["sha256"],
                                   calibration=v["cal"]) for k, v in phys.items()},
                calibration_fold=CAL_FOLD, thresholds=dict(dalpha_min=DALPHA_MIN, log_bf_min=LOG_BF_MIN, p_support=P_SUPPORT,
                                                           p_none=P_NONE, ppp_min=PPP_MIN, ndwi_ice_water=NDWI_ICE_WATER,
                                                           snow_b3=SNOW_B3, snow_b8=SNOW_B8,
                                                           prior_dominated=[PRIOR_DOMINATED_SD_RATIO, PRIOR_DOMINATED_SHIFT]),
                diagnostics=diag, code=code_provenance(), aoi=dict(lat=lat, lon=lon, size_km=size_km, resolution_m=resolution),
                deviations=["D1 posterior-mean optics: calibration uncertainty only via tau",
                            "D3 k prior includes S6 2017 plot anisotropy"],
                validation_status="research-grade; pixel-scale products not independently validated")
    stem = scene_id
    paths = export(maps, scene.transform, scene.crs, outdir, stem, prov)
    checks = validate_outputs(paths, maps, scene.transform, scene.crs)
    rep = report(maps, diag, prov, resolution) + f"\n## Output validation\n\n```\n{json.dumps(checks)}\n```\n"
    with open(os.path.join(outdir, f"{stem}_report.md"), "w") as fh:
        fh.write(rep)
    return maps, paths, checks, rep


def main(argv=None):
    import s2_io as io
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--scene-id", default=io.DEFAULT_SCENE)
    p.add_argument("--lat", type=float, default=io.DEFAULT_AOI["lat"])
    p.add_argument("--lon", type=float, default=io.DEFAULT_AOI["lon"])
    p.add_argument("--size-km", type=float, default=3.0)
    p.add_argument("--resolution", type=float, default=20.0)
    p.add_argument("--outdir", default=os.path.join(HERE, "results", "scene_product"))
    p.add_argument("--write-schema", default=None, help="write the JSON output contract to this path and exit")
    a = p.parse_args(argv)
    if a.write_schema:
        with open(a.write_schema, "w") as fh:
            json.dump(schema(), fh, indent=1)
        return
    maps, paths, checks, rep = run(a.scene_id, a.lat, a.lon, a.size_km, a.resolution, a.outdir)
    print(rep)


if __name__ == "__main__":
    main()
