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
    "halbach2022": "Halbach L., Chevrollier L.-A. et al. (2022) Pigment signatures of algal communities and their "
                   "implications for glacier surface darkening. Sci. Rep. 12:17643, doi:10.1038/s41598-022-22271-4",
    "chevrollier2023": "Chevrollier L.-A. et al. (2023) Light absorption and albedo reduction by pigmented microalgae "
                       "on snow and ice. J. Glaciol. 69:333-341, doi:10.1017/jog.2022.64",
    "prochazkova2025": "Prochazkova L. et al. (2025) Phenolic iron complexes protect glacier ice algae "
                       "(Zygnematophyceae) against excessive UV and VIS irradiation. Environ. Microbiol. Rep. "
                       "17:e70149, doi:10.1111/1758-2229.70149 (CC BY 4.0)",
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
# Halbach et al. (2022): Greenland (Mittivakkat) mean cell volumes, um^3
_VOL_GREENLAND = {"nordenskioeldii": (2307.0, 421.0), "alaskanum": (822.0, 133.0)}


def species_geometry():
    """Cylinder dimensions per species: Greenland mean cell volume (Halbach et al. 2022) with the
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
    per cell scales with cell volume (supported in direction, not significance, by Halbach et al.
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
    (0.50, "halbach2022", "Mittivakkat, E Greenland: 'almost equal proportions of both species'"),
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
"""Support of the bubble-radius grid (not a prior): the brightest counted field spectrum (B8 = 0.88)
needs r ~ 300 um at the measured densities, the darkest (B8 ~ 0.11) ~20 mm. The prior inside it is
ice_ssa_prior() converted with bubble_radius_um()."""


def bubble_radius_um(ssa_m2_kg, rho_kg_m3):
    """BioSNICAR bubbly-ice radius with the same specific surface area (m^2 per kg of ice). BioSNICAR
    (column_OPs.get_layer_OPs) scatters with sigma_air(r) * V_air / rho per kg, V_air = 1 - rho/917,
    i.e. bubble surface 3 V_air / (r rho) per kg (Whicker et al. 2022, Eq. 5)."""
    return 3.0 * (1.0 - rho_kg_m3 / 917.0) / (rho_kg_m3 * np.asarray(ssa_m2_kg, float)) * 1e6


def ssa_from_bubble_radius(r_um, rho_kg_m3):
    return 3.0 * (1.0 - rho_kg_m3 / 917.0) / (rho_kg_m3 * np.asarray(r_um, float) * 1e-6)


@lru_cache(maxsize=None)
def ice_ssa_prior():
    """Log-normal prior for the specific surface area of bubbly bare ice, fitted to all measured
    values found (data/empirical/ice_ssa_measurements.csv):
      Cooper et al. (2021), W Greenland (840 m a.s.l.): 0.35 and 0.31 m^2/kg (in-ice light attenuation);
      Dadic et al. (2013), Antarctic blue/white ice, micro-CT, 17 subsamples with rho > 830 kg/m^3.
    The two studies get equal total weight (Greenland is otherwise outnumbered 17:2).
    Returns (mean ln SSA, SD ln SSA, n values). SSA, not bubble radius, carries over between
    densities: micro-CT bubble radii are nearly constant while SSA tracks the air fraction."""
    d = pd.read_csv(path("ice_ssa_measurements.csv"))
    x = np.log(d.ssa_m2_kg.to_numpy())
    w = d.source.map(1.0 / d.source.value_counts()).to_numpy()
    w = w / w.sum()
    m = float(np.sum(w * x))
    return m, float(np.sqrt(np.sum(w * (x - m) ** 2))), int(len(x))


def bubble_radius_prior(rho_kg_m3: float = None):
    """(mean, SD) of ln(bubble radius / um) at the given layer density (default: near-surface ice,
    Cooper et al. 2018), from ice_ssa_prior()."""
    rho = ICE_DENSITY["near_surface"]["mean"] if rho_kg_m3 is None else rho_kg_m3
    m, sd, _ = ice_ssa_prior()
    return float(np.log(bubble_radius_um(np.exp(m), rho))), sd


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
COUNT_RESOLUTION_CELLS_ML = 62.5   # Fuchs-Rosenthal, 80 large squares: cells/mL = cells counted x 5000 / 80


def field_samples(include_stibal: bool = True):
    """Co-located field samples, each with a nadir HCRF spectrum and a haemocytometer count of the
    algae in the ice under the spectrometer. Returns (table, {dataset: spectra DataFrame}).

    s6_2017   Cook et al. (2020) data archive (Zenodo 10.5281/zenodo.3564501): 47 samples, S6,
              13-24 July 2017, counts from the primary workbook (cells counted in 80 large squares ->
              Poisson counting error), HCRF from HCRF_master.csv, solar zenith at noon (measurements
              'at solar noon +/- 2 h'). Two samples present in biosnicar-py's metadata with 0 cells
              (22_7_SB6, 22_7_SB7) are not in the count workbook and are therefore not used.
    sgris_2021 Chevrollier et al. (2023) Table S4 (counts, solar zenith of each spectrum) with spectra
              from Zenodo 10.5281/zenodo.18826013: 18 counted samples, southern Greenland ice sheet,
              5-6 Aug 2021, surface scraped 1-6 cm; quintuplicate counts (number counted not given).
    s6_2014   Stibal et al. (2017), only if the files were added by hand (phase2/stibal2017.py)."""
    c = pd.read_csv(path("cook2020_archive_cell_counts.csv"))
    rows = [dict(dataset="s6_2017", sample=r.sample, cells=r.cells_per_ml, cells_counted=r.cells_counted,
                 sza=solar_zenith_noon(pd.Timestamp(year=2017, month=int(r.sample.split("_")[1]),
                                                    day=int(r.sample.split("_")[0])).dayofyear))
            for r in c.itertuples()]
    t = pd.read_csv(path("chevrollier2023_tableS4.csv"))
    t = t[t.measured & t.sample_id.str.match(r"^\d{6}-S")]
    rows += [dict(dataset="sgris_2021", sample=r.sample_id, cells=float(r.ia_cells_per_ml), cells_counted=np.nan,
                  sza=float(r.sza_deg)) for r in t.itertuples()]
    spectra = {"s6_2017": pd.read_csv(path("cook2020_archive_hcrf.csv")),
               "sgris_2021": pd.read_csv(path("chevrollier2023_hcrf.csv"))}
    tab = pd.DataFrame(rows)
    tab["quantity"] = "hcrf"
    tab.attrs["albedo_k_sd"] = np.nan
    if include_stibal:
        import stibal2017
        if stibal2017.available():
            st, sp, meta = stibal2017.load()          # raises StibalDataError with a precise message
            tab = pd.concat([tab, st], ignore_index=True)
            spectra["s6_2014"] = sp
            tab.attrs["albedo_k_sd"] = meta["albedo_k_sd"]
    return tab, spectra


def cook2020_published_inversion():
    """BioSNICAR-GO retrievals published with Cook et al. (2020) for 31 S6 samples (biosnicar-py
    Spectra_Metadata.csv, column algae_cells_inv_model_S2), as a literature baseline."""
    m = pd.read_csv(path("cook2020_field_metadata.csv"))
    return dict(zip(m.filename, m.algae_cells_inv_model_S2))


def solar_zenith_noon(day_of_year: int, lat: float = FIELD_SITE["lat"]):
    """Solar zenith angle at local solar noon (Cooper 1969 declination; Cook et al. (2020) measured
    'at solar noon +/- 2 h' under clear skies)."""
    decl = 23.45 * np.sin(np.radians(360.0 * (284 + day_of_year) / 365.0))
    return float(abs(lat - decl))


# --------------------------------------------------------------------------- #
# Phenolic chromophore spectra (TD-DFT calibration, tier D)                    #
# --------------------------------------------------------------------------- #
PHENOL_UNITS = ("phenol equivalents: Williamson et al. (2020) quantified phenolics with US EPA Method 420.1 "
                "(4-AAP), so their per-cell phenolic mass and extract MAC are both per mass of phenol "
                "equivalents")
HPLC_CHROMOPHORE_PEAKS = ("peak2", "peak3", "peak4")


@lru_cache(maxsize=None)
def chromophore_hplc_shape(lo=265.0, hi=600.0, peaks=HPLC_CHROMOPHORE_PEAKS):
    """Absorption spectra of the chromatographically isolated glacier-algal phenolics (HPLC diode
    array; Williamson et al. 2020 deposit, phenolics_hplc_abs_spec.csv). Peaks 2-4 share the
    purpurogallin-type spectrum (maximum ~305 nm, shoulder ~395 nm); peak 1 (a different compound,
    no 395-nm band) is excluded. Each peak is normalised to unit area over [lo, hi]; returns
    (wavelength_nm, mean shape, SD between peaks). The isolated compound is the uncomplexed form:
    Prochazkova et al. (2025) note it is only yellowish, while whole extracts are dark."""
    d = pd.read_csv(path("williamson2020_phenolics_hplc_abs_spec.csv"))
    d = d[(d.wave >= lo) & (d.wave <= hi)]
    wl = d.wave.to_numpy(float)
    S = np.array([d[p].to_numpy(float) / np.trapezoid(d[p].to_numpy(float), wl) for p in peaks])
    return wl, S.mean(axis=0), S.std(axis=0, ddof=1)


@lru_cache(maxsize=None)
def phenolic_extract_mac():
    """Mass absorption coefficient of whole phenolic extracts of S6 surface ice (n = 53 samples
    regressed per wavelength; Williamson et al. 2020, phenolics_mac.csv), Napierian m^2 kg^-1 per
    kg of phenol equivalents (see PHENOL_UNITS), with its regression standard error.
    Returns (wavelength_nm, mac, mac_err)."""
    d = pd.read_csv(path("williamson2020_phenolics_mac_with_error.csv"))
    return d.wave.to_numpy(float), d["m2.mg"].to_numpy(float) * 1e6, d["m2.mg.err"].to_numpy(float) * 1e6


@lru_cache(maxsize=None)
def fe_purpurogallin_absorbance():
    """Absorbance of purpurogallin (PG) and of its Fe complex (PG-Fe) in water, 255-800 nm,
    digitised from Prochazkova et al. (2025) Environ. Microbiol. Rep. 17:e70149, Fig. 4 (image file EMI4-17-e70149-g007.jpg; CC BY 4.0;
    digitisation: axis ticks located in the image, curve pixels by colour, median per column;
    reproduces the published curves to within the line width, ~0.01 absorbance).
    Returns (wavelength_nm, A_PG, A_PG_Fe)."""
    d = pd.read_csv(path("prochazkova2025_fig4_PG_PGFe_absorbance_digitized.csv"))
    return d.wavelength_nm.to_numpy(float), d.A_PG.to_numpy(float), d.A_PG_Fe.to_numpy(float)


FE_NORM_RANGE_NM = (265.0, 600.0)


def fe_increment(wl_nm):
    """Fe-induced change of purpurogallin absorption per unit integrated PG absorbance (nm^-1):
    D(l) = [A_PG-Fe(l) - A_PG(l)] / integral_{265}^{600} A_PG dl, from Prochazkova et al. (2025) Fig. 4.
    ASSUMPTION (stated): both solutions in Fig. 4 hold the same amount of purpurogallin. The paper gives
    no concentrations; the equal ~315-nm maxima (1.04 vs 1.05) are consistent with equal amounts but
    could also mean the curves were normalised. Either way D is the Fe-induced change at equal
    UV-band absorbance; the fitted complexed fraction phi absorbs any residual scale difference. Zero outside 255-800 nm
    (NaN below 255 nm: no data)."""
    w, a_pg, a_fe = fe_purpurogallin_absorbance()
    a_pg = np.clip(a_pg, 0.0, None)
    k = (w >= FE_NORM_RANGE_NM[0]) & (w <= FE_NORM_RANGE_NM[1])
    d = (a_fe - a_pg) / np.trapezoid(a_pg[k], w[k])
    return np.interp(np.asarray(wl_nm, float), w, d, left=np.nan, right=0.0)


def fe_complex_shape(wl_nm):
    """Unit-peak absorption shape of the Fe-purpurogallin complex (Prochazkova et al. 2025, Fig. 4),
    zero beyond the measured range. Used as the measured spectral signature of the complexed
    pigment (tier D)."""
    w, _, a = fe_purpurogallin_absorbance()
    a = np.clip(a, 0.0, None) / np.nanmax(a)
    return np.interp(np.asarray(wl_nm, float), w, a, left=np.nan, right=0.0)


# --------------------------------------------------------------------------- #
# Pigment concentration vs cell size (species-specific concentrations)          #
# --------------------------------------------------------------------------- #
@lru_cache(maxsize=None)
def phenolic_size_scaling():
    """Test of the equal-concentration assumption with 64 S6 samples that have both phenolics per
    cell and mean biovolume per cell (Williamson et al. 2020, phenolic_pigmentation_concentrations).

    Model: ln(pg/cell) = a + b ln(V/V_ref) + e, e ~ N(0, 1/N_i + s0^2), where N_i = counted-cell proxy
    (cells mL^-1 / 62.5) gives the Poisson counting error of the per-cell normalisation and s0 the
    extra scatter; (a, b, s0) by maximum likelihood, SEs from the observed information.
    Equal intracellular concentration <=> b = 1. Concentration scales as V^gamma, gamma = b - 1.
    Returns dict(gamma, gamma_se, b, a, s0, n, v_ref)."""
    from scipy import optimize
    d = pd.read_csv(path("williamson2020_phenolics_per_sample.csv"))
    d = d[["algal.cells.per.ml", "um3.per.cell", "ng.per.cell.conc2"]].dropna()
    v_ref = s6_biovolume_um3()[0]
    x_all = np.log(d["um3.per.cell"].to_numpy() / v_ref)
    y_all = np.log(d["ng.per.cell.conc2"].to_numpy())
    n_all = d["algal.cells.per.ml"].to_numpy() / COUNT_RESOLUTION_CELLS_ML

    def fit(keep):
        x, y, n_c = x_all[keep], y_all[keep], n_all[keep]

        def nll(t):
            var = 1.0 / n_c + np.exp(2 * t[2])
            return 0.5 * np.sum((y - t[0] - t[1] * x) ** 2 / var + np.log(var))
        return optimize.minimize(nll, [y.mean(), 1.0, np.log(0.3)], method="Nelder-Mead",
                                 options=dict(xatol=1e-9, fatol=1e-12, maxiter=20000)).x

    n = len(y_all)
    a, b, ls0 = fit(np.ones(n, bool))
    # Jackknife SE: the estimate is sensitive to individual low-abundance samples (leave-one-out range
    # about -1.3 to -0.6; restricting to >= 300 cells/mL flips it to about +0.4), which the
    # information-matrix SE (about 0.47) does not capture.
    jk = np.array([fit(np.arange(n) != i)[1] for i in range(n)])
    se = float(np.sqrt((n - 1) / n * np.sum((jk - jk.mean()) ** 2)))
    return dict(gamma=float(b - 1.0), gamma_se=se, b=float(b), a=float(a), s0=float(np.exp(ls0)), n=int(n),
                v_ref=float(v_ref))


def species_concentrations_kg_m3(pigment: str = "phenolics", gamma: float | None = None):
    """Intracellular concentration per species: c_s = c_pooled (V_s / V_pooled)^gamma, gamma from
    phenolic_size_scaling() (measured; applied to every pigment, as only phenolics have per-sample
    size data). With gamma = 0 this reduces to the equal-concentration case."""
    g = phenolic_size_scaling()["gamma"] if gamma is None else gamma
    c = intracellular_concentration_kg_m3(pigment)
    v_ref = s6_biovolume_um3()[0]
    return {s: c * (geo["volume_um3"] / v_ref) ** g for s, geo in species_geometry().items()}


# --------------------------------------------------------------------------- #
# Clear-sky transmissivity (forcing in W m^-2)                                  #
# --------------------------------------------------------------------------- #
SOLAR_CONSTANT = 1361.0          # W m^-2 (Kopp & Lean 2011)
KAN_M = dict(lat=67.067, lon=-48.83, elev_m=1270.0)


def _solar_mu(time_utc, lat, lon):
    """cos(SZA) with the NOAA (Spencer 1971) declination / equation-of-time series."""
    t = pd.to_datetime(time_utc)
    doy = t.dt.dayofyear.to_numpy()
    hr = (t.dt.hour + t.dt.minute / 60.0).to_numpy()
    g = 2 * np.pi / 365 * (doy - 1 + (hr - 12) / 24)
    decl = (0.006918 - 0.399912 * np.cos(g) + 0.070257 * np.sin(g) - 0.006758 * np.cos(2 * g)
            + 0.000907 * np.sin(2 * g) - 0.002697 * np.cos(3 * g) + 0.00148 * np.sin(3 * g))
    eot = 229.18 * (0.000075 + 0.001868 * np.cos(g) - 0.032077 * np.sin(g) - 0.014615 * np.cos(2 * g)
                    - 0.040849 * np.sin(2 * g))
    ha = np.radians((hr * 60 + eot + 4 * lon) / 4 - 180)
    la = np.radians(lat)
    return np.sin(la) * np.sin(decl) + np.cos(la) * np.cos(decl) * np.cos(ha), doy


@lru_cache(maxsize=None)
def clear_sky_transmissivity(cc_max: float = 0.1, sza_max: float = 75.0):
    """Bulk clear-sky transmissivity T in SW = S0 E0 cos(SZA) T^(1/cos SZA), fitted to PROMICE KAN_M
    hourly tilt-corrected downwelling shortwave (dsr_cor), June-August of every available year, in
    hours with PROMICE cloud cover cc <= cc_max (cc is derived from downwelling longwave) and SZA <=
    sza_max. Hourly values are averages over the hour, timestamped at its start (centre used).
    Returns (median T, SD of per-hour T, n hours)."""
    h = pd.read_csv(path("promice_KAN_M_hour_JJA_radiation.csv"), parse_dates=["time"])
    mu, doy = _solar_mu(h.time + pd.Timedelta(minutes=30), KAN_M["lat"], KAN_M["lon"])
    e0 = 1 + 0.033 * np.cos(2 * np.pi * doy / 365)
    ok = (h.cc.to_numpy() <= cc_max) & (mu > np.cos(np.radians(sza_max))) & (h.dsr_cor.to_numpy() > 0)
    lt = np.log(h.dsr_cor.to_numpy()[ok] / (SOLAR_CONSTANT * e0[ok] * mu[ok])) * mu[ok]
    return float(np.exp(np.median(lt))), float(np.std(np.exp(lt), ddof=1)), int(ok.sum())


# --------------------------------------------------------------------------- #
# Mineral dust at S6                                                            #
# --------------------------------------------------------------------------- #
DUST_S6 = dict(
    lap_ug_per_ml=394.0, lap_sd_ug_per_ml=194.0, inorganic_fraction=0.95, mean_ug_per_g=342.0, max_ug_per_g=519.0,
    quote="They measured 394 ± 194 µgLAP mLice−1, of which ∼ 95 % was inorganic, giving mean and maximum mineral "
          "dust loadings of 373 and 567 µgLAP mLice−1. Assuming 1 mL of ice to weigh 0.917 g, this gives mean and "
          "maximum mass mixing ratios of 342 and 519 µgdust gice−1.",
    source="Cook et al. (2020) Sect. 2.5, citing McCutcheon et al. (heavy-biomass samples, S6, 2017)")


def dust_prior():
    """Log-normal prior for the mineral-dust mass mixing ratio of the surface layer (ppb = ng g^-1),
    moment-matched to the measured S6 loading: mean 342 ug/g, relative SD 194/394 (Cook et al. 2020).
    Returns (mean ln ppb, SD ln ppb)."""
    m = DUST_S6["mean_ug_per_g"] * 1e3
    cv = DUST_S6["lap_sd_ug_per_ml"] / DUST_S6["lap_ug_per_ml"]
    s2 = np.log(1.0 + cv ** 2)
    return float(np.log(m) - 0.5 * s2), float(np.sqrt(s2))
