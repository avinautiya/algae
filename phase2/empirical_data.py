"""
Empirical inputs for every phase, loaded from the raw source files in data/empirical/
(see data/empirical/SOURCES.md for citations, licences and verbatim quotes).

Nothing here is tuned by hand: each value is either read directly from a published
dataset or computed from it by the documented operation. Where no empirical value
exists, the function says so in its docstring and the caller must treat it as an
explicit assumption.
"""

from __future__ import annotations

import os
from functools import lru_cache

import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "data", "empirical"))
WVL_480_NM = np.round(np.arange(0.205, 4.999, 0.01), 3) * 1000.0

CITATIONS = {
    "williamson2020": "Williamson C.J. et al. (2020) Algal photophysiology drives darkening and melt of the "
                      "Greenland Ice Sheet. PNAS 117:5694-5705, doi:10.1073/pnas.1918412117; data: UK Polar Data "
                      "Centre GB/NERC/BAS/PDC/01248 (OGL v3)",
    "williamson2018": "Williamson C.J. et al. (2018) Ice algal bloom development on the surface of the Greenland Ice "
                      "Sheet. FEMS Microbiol. Ecol. 94:fiy025, doi:10.1093/femsec/fiy025",
    "chevrollier2022": "Chevrollier L.-A. et al. (2022) Pigment signatures of algal communities and their "
                       "implications for glacier surface darkening. Sci. Rep. 12:18101, doi:10.1038/s41598-022-22271-4",
    "prochazkova2021": "Prochazkova L. et al. (2021) Unicellular versus filamentous: the glacial alga Ancylonema "
                       "alaskana comb. et stat. nov. ... Microorganisms 9:1103, doi:10.3390/microorganisms9051103",
    "halbach2025": "Halbach L. et al. (2025) Single-cell imaging reveals efficient nutrient uptake and growth of "
                   "microalgae darkening the Greenland Ice Sheet. Nat. Commun. 16, doi:10.1038/s41467-025-56664-6",
    "cook2020": "Cook J.M. et al. (2020) Glacier algae accelerate melt rates on the south-western Greenland Ice "
                "Sheet. The Cryosphere 14:309-330, doi:10.5194/tc-14-309-2020; data distributed with "
                "biosnicar-py (data/additional_data, MIT licence)",
    "cooper2018": "Cooper M.G. et al. (2018) Meltwater storage in low-density near-surface bare ice in the "
                  "Greenland ice sheet ablation zone. The Cryosphere 12:955-970, doi:10.5194/tc-12-955-2018",
    "cooper2021": "Cooper M.G. et al. (2021) Spectral attenuation coefficients from measurements of light "
                  "transmission in bare ice on the Greenland Ice Sheet. The Cryosphere 15:1931-1953, "
                  "doi:10.5194/tc-15-1931-2021",
    "esa_srf": "ESA (2024) Sentinel-2 Spectral Response Functions, COPE-GSEG-EOPG-TN-15-0007 v4.0 "
               "(sentiwiki.copernicus.eu)",
    "promice": "Fausto R.S. et al. (2021) PROMICE automatic weather station data, GEUS Dataverse "
               "doi:10.22008/FK2/IW73UU (CC-BY 4.0); stations KAN_L, KAN_M, KAN_U",
}


def path(name):
    return os.path.join(ROOT, name)


# --------------------------------------------------------------------------- #
# Sentinel-2 spectral response (ESA TN-15-0007 v4.0)                          #
# --------------------------------------------------------------------------- #
@lru_cache(maxsize=None)
def s2_srf_480(spacecraft: str = "S2A", bands=("B2", "B3", "B4", "B8")):
    """Official ESA spectral responses averaged into BioSNICAR's 10-nm bins (centre +/- 5 nm).
    Returns dict band -> (480,) weights (not normalised)."""
    d = pd.read_csv(path("esa_s2_srf_TN-15-0007_v4.0.csv")).set_index("wavelength_nm")
    out = {}
    for b in bands:
        s = d[f"{spacecraft}_{b}"].astype(float)
        w = np.zeros(WVL_480_NM.size)
        for i, c in enumerate(WVL_480_NM):
            seg = s.loc[(s.index >= c - 5) & (s.index < c + 5)]
            w[i] = seg.mean() if len(seg) else 0.0
        out[b] = np.nan_to_num(w)
    return out


def spacecraft_from_scene(scene_id: str) -> str:
    sid = (scene_id or "").upper()
    for s in ("S2A", "S2B", "S2C"):
        if sid.startswith(s):
            return s
    return "S2A"


# --------------------------------------------------------------------------- #
# Pigments (Williamson et al. 2020, S6, 2016)                                  #
# --------------------------------------------------------------------------- #
CAROTENOIDS = ("antheraxanthin", "carotene", "lutein", "zeaxanthin", "violaxanthin", "neoxanthin")


@lru_cache(maxsize=None)
def pigments_per_cell():
    """Mean (and SD) pigment mass per glacier-algal cell, ng cell^-1, from 53 surface-ice samples
    (file average_pigment_concentrations.csv). Carotenoids are summed, as in Williamson et al.
    (2020), who treated all carotenoids as photoprotective."""
    p = pd.read_csv(path("williamson2020_pigments_per_cell.csv")).set_index("pigment")
    return {
        "phenolics": (p.loc["phenolics", "ng.per.cell"], p.loc["phenolics", "ng.per.cell.sd"]),
        "chla": (p.loc["chla", "ng.per.cell"], p.loc["chla", "ng.per.cell.sd"]),
        "chlb": (p.loc["chlb", "ng.per.cell"], p.loc["chlb", "ng.per.cell.sd"]),
        "carotenoids": (p.loc[list(CAROTENOIDS), "ng.per.cell"].sum(),
                        float(np.sqrt((p.loc[list(CAROTENOIDS), "ng.per.cell.sd"] ** 2).sum()))),
    }


@lru_cache(maxsize=None)
def pigment_macs_480():
    """In vivo mass absorption coefficients (Napierian, m^2 kg^-1) on the 480-band grid.

    chla, chlb, carotenoids (photoprotective form, 'ppc') and phenolics, as used by
    Williamson et al. (2020) (phenolics derived there; chl/carotenoid spectra after
    Dauchet et al. 2015). Defined 250-750 nm; set to zero outside (no data)."""
    m = pd.read_csv(path("williamson2020_mass_absorption_coefficients.csv"))
    wl = m.wave.to_numpy(float)
    out = {}
    for key, col in (("chla", "chla.m2.mg"), ("chlb", "chlb.m2.mg"), ("carotenoids", "ppc.m2.mg"),
                     ("phenolics_williamson2020", "phenols.m2.mg")):
        out[key] = np.interp(WVL_480_NM, wl, m[col].to_numpy(float), left=0.0, right=0.0) * 1e6
    return out


# --------------------------------------------------------------------------- #
# Cells                                                                         #
# --------------------------------------------------------------------------- #
@lru_cache(maxsize=None)
def s6_biovolume_um3():
    """Pooled mean glacier-algal biovolume per cell at S6 (2016): sum(um^3 mL^-1)/sum(cells mL^-1)
    over natural surface-ice samples, the same assemblage whose pigments per cell were measured."""
    c = _s6_surface_counts()
    c = c[c["overall.cells.per.ml"] > 0]
    pooled = c["total.biomass.um3.ml"].sum() / c["overall.cells.per.ml"].sum()
    per = c["total.biomass.um3.ml"] / c["overall.cells.per.ml"]
    return float(pooled), float(per.mean()), float(per.std(ddof=1)), int(len(per))


# Procházková et al. (2021) Table 2: width x length means per population (um)
_PRO_ALASKANA = [(8.7, 12.1), (7.6, 12.1), (8.5, 12.8), (8.4, 11.4)]
_PRO_NORDENSKIOELDII = [(12.9, 29.4), (13.6, 37.5), (13.7, 28.3)]
# Chevrollier et al. (2022): Greenland (Mittivakkat) mean cell volumes, um^3
_VOL_GREENLAND = {"nordenskioeldii": (2307.0, 421.0), "alaskanum": (822.0, 133.0)}


def species_geometry():
    """Cylinder dimensions per species: Greenland mean cell volume (Chevrollier et al. 2022) with the
    mean length:width ratio of measured populations (Prochazkova et al. 2021, Table 2)."""
    out = {}
    for name, dims in (("nordenskioeldii", _PRO_NORDENSKIOELDII), ("alaskanum", _PRO_ALASKANA)):
        ar = float(np.mean([L / w for w, L in dims]))
        V, Vsd = _VOL_GREENLAND[name]
        d = (4.0 * V / (np.pi * ar)) ** (1.0 / 3.0)
        out[name] = dict(diameter_um=d, length_um=ar * d, volume_um3=V, volume_sd_um3=Vsd, aspect=ar)
    return out


def intracellular_concentration_kg_m3(pigment: str = "phenolics"):
    """Pigment mass per cell / pooled S6 biovolume per cell (both Williamson et al. 2020, S6 2016).
    ASSUMPTION (stated): the intracellular concentration is the same in both species, i.e. pigment
    per cell scales with cell volume (supported in direction, not significance, by Chevrollier et al.
    2022: purpurogallin per cell vs biovolume, R^2 = 0.49, p = 0.074)."""
    ng = pigments_per_cell()[pigment][0]
    V = s6_biovolume_um3()[0]
    return ng * 1e-12 / (V * 1e-18)                  # ng -> kg ; um^3 -> m^3


# --------------------------------------------------------------------------- #
# Abundance and community priors                                               #
# --------------------------------------------------------------------------- #
_SURFACE_HABITATS = ("h", "m", "l", "clean", "dark", "dirty", "snow", "h20log")


def _s6_surface_counts():
    c = pd.read_csv(path("williamson2020_count_biomass_S6_2016.csv"))
    return c[c.habitat.isin(_SURFACE_HABITATS)]     # excludes incubations, water, cryoconite, biofilm


def abundance_prior():
    """Normal fit to log10(cells mL^-1) of natural surface-ice samples with cells > 0 at S6 in 2016
    (Williamson et al. 2020 deposit). Independent of the 2017 field-validation samples."""
    c = _s6_surface_counts()
    x = np.log10(c.loc[c["overall.cells.per.ml"] > 0, "overall.cells.per.ml"])
    return float(x.mean()), float(x.std(ddof=1)), int(len(x)), float((c["overall.cells.per.ml"] == 0).mean())


COMMUNITY_SURVEYS = [  # fraction of A. nordenskioeldii among Ancylonema cells
    (0.65, "williamson2018", "SW Greenland K-transect 2016: 'Ancylonema ... (~65%) ... followed by Mesotaenium (~35%)'"),
    (0.66, "halbach2025", "SW tip of GrIS (QAS_M) 2020: '~66% were filamentous Ancylonema cf. nordenskioldii'"),
    (0.50, "chevrollier2022", "Mittivakkat, E Greenland: 'almost equal proportions of both species'"),
]


def community_prior():
    """Beta distribution matched (method of moments) to the three published Greenland surveys."""
    f = np.array([v for v, *_ in COMMUNITY_SURVEYS])
    mu, sd = f.mean(), f.std(ddof=1)
    k = mu * (1 - mu) / sd ** 2 - 1
    return float(mu * k), float((1 - mu) * k), float(mu), float(sd)


# --------------------------------------------------------------------------- #
# Ice and illumination                                                          #
# --------------------------------------------------------------------------- #
ICE_DENSITY = {  # Cooper et al. (2018), 67.049 N 49.022 W, 1215 m, July 2016 (g cm^-3 -> kg m^-3)
    "weathering_crust": dict(mean=450.0, lo=330.0, hi=560.0,
                             quote="thin (<0.5 m), even lower density (0.33-0.56 g cm-3, mu = 0.45 g cm-3) "
                                   "unsaturated weathering crust"),
    "near_surface": dict(mean=690.0, lo=430.0, hi=910.0,
                         quote="low-density (0.43-0.91 g cm-3, mu = 0.69 g cm-3) ice to at least 1.1 m depth"),
}

ICE_RADIUS_BOUNDS_UM = (300.0, 20000.0)
"""Bubbly-ice optical radius range. Lower bound: the brightest counted field spectrum (B8 = 0.88)
needs r ~ 300 um at the measured densities; upper bound: darkest NIR (B8 ~ 0.11) needs ~20 mm.
Cooper et al. (2021) measured optical effective radii of ~9.3 and ~10.6 mm in this ice."""


@lru_cache(maxsize=None)
def anisotropy_prior(bands=("B2", "B3", "B4", "B8"), spacecraft="S2A"):
    """Mean and SD of the field anisotropic reflectance factor (HCRF / albedo) band-averaged with the
    ESA SRFs over all spectra in biosnicar-py's ARF_master.csv (Cook et al. field campaigns). This is
    the empirical prior for the multiplicative factor k linking directional reflectance to albedo."""
    a = pd.read_csv(path("biosnicar_field_ARF.csv"))
    wl = a.wavelength_nm.to_numpy()
    srf = pd.read_csv(path("esa_s2_srf_TN-15-0007_v4.0.csv")).set_index("wavelength_nm")
    vals = []
    for col in a.columns[1:]:
        sp = a[col].to_numpy(float)
        per = []
        for b in bands:
            s = srf[f"{spacecraft}_{b}"].reindex(wl).fillna(0).to_numpy()
            ok = np.isfinite(sp)
            per.append(np.sum(sp[ok] * s[ok]) / s[ok].sum())
        vals.append(np.mean(per))
    v = np.array(vals)
    return float(v.mean()), float(v.std(ddof=1)), int(v.size)


# --------------------------------------------------------------------------- #
# PROMICE air temperature (melt diagnostics)                                    #
# --------------------------------------------------------------------------- #
def promice_daily(station: str):
    return pd.read_csv(path(f"promice_{station}_day_2019.csv"), parse_dates=["time"])


def pdd_at_elevation(elev_m, start="2019-05-01", end="2019-07-23", stations=("KAN_L", "KAN_M", "KAN_U")):
    """Positive degree days at arbitrary elevations: each day's air temperature is linearly
    interpolated (or extrapolated) in elevation between the PROMICE stations' daily means."""
    tabs = [promice_daily(s) for s in stations]
    alts = np.array([t.gps_alt.median() for t in tabs])
    days = pd.date_range(start, end, freq="D")
    T = np.full((len(days), len(stations)), np.nan)
    for j, t in enumerate(tabs):
        s = t.set_index(t.time.dt.normalize()).t_u
        T[:, j] = s.reindex(days).to_numpy()
    z = np.asarray(elev_m, dtype=float)
    pdd = np.zeros_like(z, dtype=float)
    for i in range(len(days)):
        ok = np.isfinite(T[i])
        if ok.sum() < 2:
            continue
        p = np.polyfit(alts[ok], T[i, ok], 1)
        pdd += np.clip(np.polyval(p, z), 0, None)
    return pdd


def lapse_rate_k_per_km(start="2019-06-01", end="2019-07-31"):
    tabs = {s: promice_daily(s) for s in ("KAN_L", "KAN_M", "KAN_U")}
    alts, tm = [], []
    for s, t in tabs.items():
        sel = t[(t.time >= start) & (t.time <= end)]
        alts.append(t.gps_alt.median())
        tm.append(sel.t_u.mean())
    return float(-np.polyfit(alts, tm, 1)[0] * 1000.0)


# --------------------------------------------------------------------------- #
# Field validation set (Cook et al. 2020)                                       #
# --------------------------------------------------------------------------- #
FIELD_SITE = dict(lat=67.04, lon=-49.07, year=2017)  # Cook et al. (2020) Sect. 2.2
COUNT_RESOLUTION_CELLS_ML = 62.5                     # 24 of 25 non-zero field counts are multiples of 62.5


def field_samples():
    """31 co-located field samples: HCRF spectrum (350-2499 nm) + haemocytometer cell count."""
    m = pd.read_csv(path("cook2020_field_metadata.csv"))
    h = pd.read_csv(path("cook2020_field_hcrf_counted_samples.csv"))
    return m, h


def solar_zenith_noon(day_of_year: int, lat: float = FIELD_SITE["lat"]):
    """Solar zenith angle at local solar noon (Cooper 1969 declination; Cook et al. (2020) measured
    'at solar noon +/- 2 h' under clear skies)."""
    decl = 23.45 * np.sin(np.radians(360.0 * (284 + day_of_year) / 365.0))
    return float(abs(lat - decl))
