#!/usr/bin/env python
"""
Verified data manifest for the glacier-darkening inversion study (milestone 1).

Every automatic field below (sha256, shape, wavelength grid, NaN and range checks, sample identifiers,
sampling groups, cross-file overlaps) is computed from the file as it exists when this script runs. The
curated fields (spectral quantity, geometry, documented uncertainty, ancestry, endpoints) are written by
hand from the sources named in `source`; where this session did not read the original publication the field
says so (`verified_here: false`).

Outputs:
  records/inversion_m1/data_manifest.json   one entry per dataset
  records/inversion_m1/sample_groups.csv    one row per sample: dataset, id, site, campaign, sampling group,
                                            which measurements exist, proposed role
Roots (environment variables; each is a pinned checkout, never modified here):
  BIOSNICAR_PATH   biosnicar-py at the revision recorded in docs/data_manifest.md (fe74eeef)
  SNOWLAPS_PATH    openosmia/snowlaps checkout (commit recorded in the manifest)
  MAIN_CHECKOUT    the main repository checkout, only for untracked read-only data (data/promice_raw)

    python inversion_m1/build_manifest.py
"""

from __future__ import annotations

import csv
import datetime as dt
import hashlib
import json
import os
import re
import subprocess
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, ".."))
OUT = os.path.join(REPO, "records", "inversion_m1")
SCR = "/tmp/claude-0/-home-user-algae/32b7214c-bd8e-5e37-94c0-4ae97c3c08dd/scratchpad"
BIOSNICAR = os.environ.get("BIOSNICAR_PATH", os.path.join(SCR, "ext", "biosnicar-py"))
SNOWLAPS = os.environ.get("SNOWLAPS_PATH", os.path.join(SCR, "ext_snowlaps", "src"))
MAIN = os.environ.get("MAIN_CHECKOUT", os.path.abspath(os.path.join(REPO, "..", "..")))


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def git_head(path: str) -> str | None:
    try:
        return subprocess.run(["git", "-C", path, "rev-parse", "HEAD"], capture_output=True, text=True,
                              timeout=20).stdout.strip() or None
    except (OSError, subprocess.SubprocessError):
        return None


def file_entry(path: str, role: str) -> dict:
    ok = os.path.isfile(path)
    return dict(role=role, path=path.replace(REPO + os.sep, "").replace(BIOSNICAR, "<BIOSNICAR_PATH>")
                .replace(SNOWLAPS, "<SNOWLAPS_PATH>").replace(MAIN, "<MAIN_CHECKOUT>"),
                exists=ok, bytes=os.path.getsize(path) if ok else None, sha256=sha256(path) if ok else None)


# ---------------------------------------------------------------------------------------------- readers
def wide_spectra(df: pd.DataFrame, wl_col: str | None, ids: list[str], window=(400.0, 1300.0)) -> dict:
    """Checks for a wavelength x sample table: grid, NaNs, physically implausible values."""
    wl = df[wl_col].to_numpy(float) if wl_col else df.index.to_numpy(float)
    sub = df[ids].to_numpy(float)
    w = (wl >= window[0]) & (wl <= window[1])
    step = np.unique(np.round(np.diff(wl), 6))
    empty = [ids[i] for i in range(len(ids)) if not np.isfinite(sub[w, i]).any()]
    bad = [i for i, c in enumerate(ids) if np.isfinite(sub[w, i]).any()
           and (np.nanmax(sub[w, i]) > 1.2 or np.nanmin(sub[w, i]) < 0.0)]
    return dict(n_samples=len(ids), wavelength_nm=dict(min=float(wl.min()), max=float(wl.max()), n=int(wl.size),
                                                        steps_nm=[float(s) for s in step[:4]]),
                nan_fraction=float(np.isnan(sub).mean()),
                fraction_samples_outside_0_1p2_in_window=round(len(bad) / max(len(ids), 1), 3),
                window_nm=list(window), samples_outside_window_qc=[ids[i] for i in bad][:12],
                samples_with_no_finite_value_in_window=empty,
                overall_min=float(np.nanmin(sub)), overall_max=float(np.nanmax(sub)))


def s6_day(sample_id: str) -> str | None:
    m = re.match(r"^(\d+)_(\d+)_", sample_id)
    return f"2017-{int(m.group(2)):02d}-{int(m.group(1)):02d}" if m else None


def sgreenland_day(sample_id: str) -> str | None:
    m = re.match(r"^(\d{2})(\d{2})(\d{2})-", sample_id)         # ddmmyy-S#
    return f"20{m.group(3)}-{m.group(2)}-{m.group(1)}" if m else None


def build() -> tuple[dict, list[dict]]:
    R = lambda p: os.path.join(REPO, p)                                      # noqa: E731
    B = lambda p: os.path.join(BIOSNICAR, "data", "additional_data", p)      # noqa: E731
    Nw = lambda p: os.path.join(SNOWLAPS, "data", p)                         # noqa: E731
    datasets, rows = [], []

    # ---------------------------------------------------------------- S6 2017 (Cook et al. 2020 archive)
    alb = pd.read_csv(B("Albedo_master.csv"))
    hc = pd.read_csv(R("data/empirical/cook2020_archive_hcrf.csv"))
    ct = pd.read_csv(R("data/empirical/cook2020_archive_cell_counts.csv"))
    arf = pd.read_csv(R("data/empirical/biosnicar_field_ARF.csv"))
    meta = pd.read_csv(R("data/empirical/cook2020_field_metadata.csv"))
    a_ids = [c for c in alb.columns if re.match(r"^\d+_\d+_", c)]
    h_ids = [c for c in hc.columns if c != "wavelength_nm"]
    c_ids = list(ct["sample"])
    r_ids = [c for c in arf.columns if c != "wavelength_nm"]
    both = sorted(set(a_ids) & set(h_ids) & set(c_ids))
    s6_days_counts = sorted({s6_day(i) for i in c_ids})
    s6_days_albedo = sorted({s6_day(i) for i in a_ids})
    overlaps = dict(counts_and_albedo=len(set(c_ids) & set(a_ids)), counts_and_hcrf=len(set(c_ids) & set(h_ids)),
                    albedo_and_hcrf=len(set(a_ids) & set(h_ids)), arf_and_albedo=len(set(r_ids) & set(a_ids)),
                    arf_and_counts=len(set(r_ids) & set(c_ids)), published_retrievals_and_counts=len(
                        set(meta["filename"]) & set(c_ids)),
                    all_three_measurements=len(both))
    datasets.append(dict(
        id="S6_2017_albedo", title="S6 2017 hemispherical spectral albedo (Cook et al. 2020 archive)",
        source="Zenodo 10.5281/zenodo.3564501 via biosnicar-py data/additional_data/Albedo_master.csv (MIT)",
        version="biosnicar-py revision " + str(git_head(BIOSNICAR)), license="MIT (biosnicar-py)",
        files=[file_entry(B("Albedo_master.csv"), "spectra")],
        site="S6, SW Greenland ablation area (about 1000 m a.s.l.; one site)", period="2017-07-13 to 2017-07-25",
        spectral_quantity="hemispherical spectral albedo (upward/downward irradiance ratio, cosine receptor)",
        geometry_illumination="natural sky at about solar noon; sky condition and exact time per plot not recorded in the file",
        units="dimensionless (0-1)", wavelengths=wide_spectra(alb, "Wavelength", a_ids),
        documented_uncertainty="about +-0.01-0.02 (not stated per plot); source docs/data_manifest.md; the publication was NOT re-read here",
        verified_here=dict(structure=True, uncertainty=False),
        measured_components="none in this file; algal counts for 47 of the 87 plots (S6_2017_counts)",
        colocation="the 47 SB plots with counts also have an albedo and an HCRF spectrum (same plot IDs); 40 plots (tag 'S') have albedo only",
        temporal_match="same day and plot; counts are surface samples taken at the plot, albedo is a footprint average",
        id_overlaps=overlaps, n_samples_total=len(a_ids), n_with_counts=len(both),
        sampling_groups=dict(level="sampling day (plots differ per day)", n_days_albedo=len(s6_days_albedo),
                             n_days_counts=len(s6_days_counts), days_with_counts=s6_days_counts),
        independent_glaciers=1,
        ancestry=["denominators of the ARF/k prior (50 of 87 plots); H1 albedo target; see docs/data_use_register.md",
                  "the S6 mean dust loading and the radius prior are campaign-level information from the same campaign"],
        leakage_notes="any S6 evaluation of a model whose k prior or dust prior came from S6 is not independent of S6",
        endpoints_available=["spectral albedo", "broadband albedo (integrated)"],
        non_plot_columns=[c for c in alb.columns if c != "Wavelength" and c not in a_ids]))
    for pid in a_ids:
        rows.append(dict(dataset="S6_2017_albedo", sample_id=pid, site="S6", campaign="Cook2020_2017",
                         sampling_group=f"S6:{s6_day(pid)}", has_albedo=True, has_hcrf=pid in h_ids,
                         has_counts=pid in c_ids, role="development" if pid in c_ids else "development_unlabelled",
                         note="" if pid in c_ids else "no algal count"))

    datasets.append(dict(
        id="S6_2017_hcrf", title="S6 2017 nadir HCRF (same plots)", source="Zenodo 10.5281/zenodo.3564501 (HCRF_master.csv)",
        version="repo copy", license="CC BY (Zenodo)", files=[file_entry(R("data/empirical/cook2020_archive_hcrf.csv"), "spectra")],
        site="S6", period="2017-07-13 to 2017-07-24",
        spectral_quantity="HCRF, nadir (ASD FieldSpec, bare fibre; directional, NOT albedo)",
        geometry_illumination="nadir view; solar zenith not recorded in the file", units="dimensionless reflectance factor",
        wavelengths=wide_spectra(hc, "wavelength_nm", h_ids),
        documented_uncertainty="not stated in the file", verified_here=dict(structure=True, uncertainty=False),
        measured_components="see S6_2017_counts", colocation="same plot IDs as the counts (47/47)", temporal_match="same day and plot",
        n_samples_total=len(h_ids), independent_glaciers=1,
        ancestry=["numerators of the ARF/k prior (29 of 51 ARF columns)", "sigma and radius prior of earlier retrievals (docs/data_use_register.md)"],
        leakage_notes="HCRF is directional; converting to albedo needs k = HCRF/albedo, which was itself estimated from these plots",
        endpoints_available=["directional spectral reflectance"], empty_columns=[c for c in h_ids if hc[c].isna().all()]))

    datasets.append(dict(
        id="S6_2017_counts", title="S6 2017 glacier-algae cell counts (47 plots)", source="Zenodo 10.5281/zenodo.3564501 (Cell_Counts)",
        version="repo copy", license="CC BY (Zenodo)", files=[file_entry(R("data/empirical/cook2020_archive_cell_counts.csv"), "table")],
        site="S6", period="2017-07-13 to 2017-07-24", spectral_quantity="n/a (component measurement)",
        geometry_illumination="surface sample, top about 2 cm, about 1 m2 plot", units="cells per mL of meltwater",
        wavelengths=None, documented_uncertainty="Poisson on cells counted (1/sqrt N); 5 zero counts; counted volume of zero counts not recorded",
        verified_here=dict(structure=True, uncertainty=False),
        measured_components="glacier-algal cells only. NO per-plot mineral dust, black carbon, cryoconite, pigment or ice-structure measurement",
        colocation="plot-level, matching albedo and HCRF IDs", temporal_match="same day",
        n_samples_total=len(c_ids), n_zero_counts=int((ct.cells_counted == 0).sum()),
        cells_per_ml=dict(min=float(ct.cells_per_ml.min()), median=float(ct.cells_per_ml.median()), max=float(ct.cells_per_ml.max())),
        independent_glaciers=1, ancestry=["evaluation target of earlier held-out and H1 tests"],
        leakage_notes="reused as an evaluation set before; not untouched",
        endpoints_available=["algal abundance (component truth for the algal term only)"]))

    datasets.append(dict(
        id="S6_2017_ARF_and_published_retrievals", title="ARF (HCRF/albedo) and published BioSNICAR retrievals for S6 2017",
        source="biosnicar-py ARF_master.csv and Spectra_Metadata.csv (MIT)", version="revision fe74eeef",
        license="MIT", files=[file_entry(R("data/empirical/biosnicar_field_ARF.csv"), "table"),
                              file_entry(R("data/empirical/cook2020_field_metadata.csv"), "table")],
        site="S6", period="2017", spectral_quantity="ARF = HCRF/albedo (derived); retrievals are MODEL OUTPUT (algae_cells_inv_model_S2 etc.)",
        geometry_illumination="n/a", units="dimensionless; cells per mL",
        wavelengths=wide_spectra(arf, "wavelength_nm", r_ids, window=(400.0, 1300.0)), documented_uncertainty="none",
        verified_here=dict(structure=True, uncertainty=False),
        measured_components="none (derived and model-output quantities)", colocation="29 of 51 ARF columns are count plots",
        temporal_match="same campaign", n_samples_total=len(r_ids), n_published_retrievals=int(len(meta)),
        independent_glaciers=1, ancestry=["k prior for every physics retrieval in this repository"],
        leakage_notes="published retrievals are a REFERENCE, never a truth target; the fitting procedure that produced them was not verified here",
        endpoints_available=[]))

    # ---------------------------------------------------------------- S Greenland 2021
    s4 = pd.read_csv(R("data/empirical/chevrollier2023_tableS4.csv"))
    sg = pd.read_csv(R("data/empirical/chevrollier2023_hcrf.csv"))
    g_ids = [c for c in sg.columns if c != "wavelength_nm"]
    s4_2021 = s4[s4.sample_id.map(sgreenland_day).notna()]
    s4_s6 = s4[s4.sample_id.map(sgreenland_day).isna()]
    datasets.append(dict(
        id="SGreenland_2021_hcrf_counts", title="Southern Greenland 2021 HCRF spectra and counts (Chevrollier et al. 2023 / 2026 data)",
        source="Zenodo 10.5281/zenodo.18826013 (HCRF netCDF, as recorded in data/empirical/SOURCES.md) and Chevrollier et al. 2023 Table S4 (parsed from the supplement text)",
        version="repo copy", license="CC BY 4.0",
        files=[file_entry(R("data/empirical/chevrollier2023_hcrf.csv"), "spectra"), file_entry(R("data/empirical/chevrollier2023_tableS4.csv"), "table")],
        site="southern Greenland (site name not recorded in the repository files; one site)", period="2021-08-05 to 2021-08-06",
        spectral_quantity="HCRF, nadir (ASD FieldSpec 4; directional, NOT albedo)",
        geometry_illumination="nadir; solar zenith per spectrum in Table S4 (sza_deg)", units="dimensionless reflectance factor",
        wavelengths=wide_spectra(sg, "wavelength_nm", g_ids),
        documented_uncertainty="counts: number counted not published (log-normal SD 0.10 dex assumed in earlier work); scraped depth 1-6 cm not recorded per sample",
        verified_here=dict(structure=True, uncertainty=False),
        measured_components="ice-algal and snow-algal cells (Table S4). No dust, BC, cryoconite or pigment per sample",
        colocation="20 spectra with matching Table S4 rows", temporal_match="same day",
        n_samples_total=len(g_ids), table_s4_rows=int(len(s4)), table_s4_rows_2021=int(len(s4_2021)),
        table_s4_rows_that_are_S6_2017_plots=sorted(s4_s6.sample_id.tolist()),
        measured_flags=s4.measured.value_counts().to_dict() if "measured" in s4 else None,
        sampling_groups=dict(level="sampling day", days=sorted({sgreenland_day(i) for i in g_ids})),
        independent_glaciers=1,
        ancestry=["BioSNICAR's default ice-algal optics (Chevrollier 2023) were measured on algae from this research programme",
                  "used as the PRIMARY target of held-out_v2 and for sigma/tau of the reverse fold (docs/data_use_register.md)",
                  "Table S4 also carries the authors' own fitted model parameters (model_density, model_reff_um): fits to these spectra"],
        leakage_notes="not untouched; 6 Table S4 rows are S6-2017 plots without spectra in this file",
        endpoints_available=["directional reflectance", "algal abundance"]))
    for i in g_ids:
        rows.append(dict(dataset="SGreenland_2021_hcrf_counts", sample_id=i, site="SGreenland2021", campaign="Chevrollier2023",
                         sampling_group=f"SGreenland2021:{sgreenland_day(i)}", has_albedo=False, has_hcrf=True,
                         has_counts=bool(i in set(s4[s4.measured == True].sample_id)),                     # noqa: E712
                         role="transfer_check_reused", note="HCRF only; reused in held-out_v2"))

    # ---------------------------------------------------------------- Norwegian snow (snowlaps dataset)
    nm = pd.read_csv(Nw("spectra/Chevrollier_et_al_2024_TC_metadata.csv"))
    ns = pd.read_csv(Nw("spectra/Chevrollier_et_al_2024_TC_spectra.csv"), index_col=0)
    res = pd.read_csv(Nw("optimization_results/results_inversion_no_anisotropy_correction_Chevrollier2024_TC.csv"), index_col=0)
    datasets.append(dict(
        id="Norway_snow_2023_hcrf", title="Finse (southern Norway) snowfield HCRF, 180 spectra (snowlaps dataset)",
        source="github.com/openosmia/snowlaps data/spectra; Zenodo 10.5281/zenodo.14639602 (per the paper)",
        version="snowlaps commit " + str(git_head(SNOWLAPS)), license="GPL-3.0 (repository)",
        files=[file_entry(Nw("spectra/Chevrollier_et_al_2024_TC_spectra.csv"), "spectra"),
               file_entry(Nw("spectra/Chevrollier_et_al_2024_TC_metadata.csv"), "table"),
               file_entry(Nw("optimization_results/results_inversion_no_anisotropy_correction_Chevrollier2024_TC.csv"), "published_result"),
               file_entry(Nw("emulator/snowlaps_emulator_v1.h5"), "emulator_weights"),
               file_entry(Nw("scaler/minmax_scaler.save"), "scaler")],
        site="Finse, southern Norway (60.57 N, 7.46 E); snow, not glacier ice", period="2023-07-08 to 2023-08-06",
        spectral_quantity="HCRF (ASD FieldSpec 4, 25 degree foreoptic, relative to Spectralon), treated as albedo by the authors",
        geometry_illumination="nadir; SZA 38-58 degrees (metadata); 88 of 180 spectra flagged diffuse (SZA set to 50 by snowlaps)",
        units="dimensionless; RAW values include non-physical numbers in atmospheric absorption bands",
        wavelengths=dict(wide_spectra(ns, None, list(ns.columns), window=(400.0, 1300.0)), n_diffuse=int((nm.diffuse == 1).sum())),
        documented_uncertainty="not stated per spectrum", verified_here=dict(structure=True, uncertainty=False),
        measured_components="NONE (no LAP concentration measured; black carbon presence not measured)",
        colocation="n/a", temporal_match="n/a", n_samples_total=int(ns.shape[1]), n_dates=int(nm.date.nunique()),
        published_result_rows=int(len(res)),
        sampling_groups=dict(level="sampling day (12 dates), one site", days=sorted(nm.date.unique().tolist())),
        independent_glaciers=0,
        ancestry=["reference data for reproducing the snowlaps retrieval only"],
        leakage_notes="snow substrate, red snow algae: NOT a glacier-ice evaluation set",
        endpoints_available=["spectral reflectance only"]))
    tag_date = dict(zip(nm["tag"], nm["date"]))        # tags are month-first (070823 = 2023-07-08); use the metadata date
    for tag in ns.columns:
        rows.append(dict(dataset="Norway_snow_2023_hcrf", sample_id=tag, site="Finse", campaign="Chevrollier2025_TC",
                         sampling_group=f"Finse:{tag_date[tag]}", has_albedo=False, has_hcrf=True, has_counts=False,
                         role="external_method_reproduction", note="snow; no LAP truth"))

    # ---------------------------------------------------------------- laboratory and auxiliary inputs
    aux = [("Williamson2020_lab", "Williamson et al. 2020 S6 2016 pigments, MACs, HPLC spectra (calibration data)", "OGL v3",
            ["williamson2020_pigments_per_cell.csv", "williamson2020_mass_absorption_coefficients.csv",
             "williamson2020_phenolics_mac_with_error.csv", "williamson2020_phenolics_hplc_abs_spec.csv",
             "williamson2020_phenolics_per_sample.csv", "williamson2020_count_biomass_S6_2016.csv"],
            "in vivo absorption of algal pigments (per mass, phenol-equivalent) and counts, S6 2016",
            "laboratory; the only MEASURED variability of algal pigment optics in the repository (53 extracts, one site)"),
           ("Ice_SSA", "Measured specific surface area of bubbly ice", "values from publications", ["ice_ssa_measurements.csv"],
            "19 SSA values (Cooper 2021; Dadic 2013)", "ice-structure prior"),
           ("S2_SRF", "Sentinel-2 spectral response functions", "Copernicus open licence", ["esa_s2_srf_TN-15-0007_v4.0.csv"],
            "ESA SRF v4.0, 1 nm", "band-integration operator")]
    for did, title, lic, fnames, content, role_note in aux:
        datasets.append(dict(id=did, title=title, source="data/empirical/SOURCES.md", version="repo copy", license=lic,
                             files=[file_entry(R("data/empirical/" + f), "table") for f in fnames], spectral_quantity=content,
                             notes=role_note, verified_here=dict(structure=True, uncertainty=False), independent_glaciers=0))

    # ---------------------------------------------------------------- station and biology endpoints
    pr = os.path.join(MAIN, "data", "promice_raw")
    datasets.append(dict(
        id="PROMICE_KAN_hourly", title="PROMICE KAN_L and KAN_M hourly (L3), 2016-2023", source="GEUS Dataverse doi:10.22008/FK2/IW73UU (live L3 product, downloaded 2026-10-09)",
        version="dataset version 40; sha256 recorded in data/promice_raw/SHA256SUMS", license="CC BY 4.0 (GEUS)",
        files=[file_entry(os.path.join(pr, "KAN_L_hour.csv"), "table"), file_entry(os.path.join(pr, "KAN_M_hour.csv"), "table")],
        site="KAN_L (about 650 m) and KAN_M (about 1270 m), SW Greenland", period="2016-2023 (JJA used)",
        spectral_quantity="broadband shortwave up/down radiation and albedo (radiometer); NO spectral information",
        geometry_illumination="hourly means; tilt-corrected only where dsr_cor exists (KAN_M none in 2019/2021)",
        verified_here=dict(structure=False, uncertainty=False),
        notes="untracked read-only data in the main checkout; hashes taken from this run",
        measured_components="none", independent_glaciers=2,
        endpoints_available=["station broadband albedo", "absorbed shortwave (with SW down)", "ablation (disagreeing sensors in 2019)"],
        ancestry=["forcing and validation endpoint of H2/H4 (frozen protocol); KAN_M also used for clear-sky transmissivity and chi"]))
    pm = pd.read_csv(R("records/prombio_matching/manifest.csv"))
    datasets.append(dict(
        id="PROMBIO", title="PROMBIO ice-algae counts under PROMICE radiometers (GEUS)", source="doi:10.22008/FK2/IJA2Y4 (2021-2023) and doi:10.22008/FK2/D0JC0P (2024)",
        version="see records/prombio_matching/manifest.csv", license="per GEUS Dataverse",
        files=[file_entry(R("records/prombio_matching/manifest.csv"), "table"), file_entry(R("records/prombio_matching/station_days.csv"), "table")],
        site="PROMICE stations (KAN_L, QAS_L/M/U, TAS_A, ...)", period="2021-2024", spectral_quantity="n/a (counts); station albedo is the radiative endpoint",
        measured_components="ice-algal cells per mL; strata incl. 'under radiometer', 'dark ice' (targeted, not an area mean)",
        colocation="sample under the radiometer; coordinates are STATION positions", temporal_match="station-day",
        manifest_rows=int(len(pm)), manifest_columns=list(pm.columns)[:14],
        roles_frozen="docs/new_data_inventory.md: 2021 and 2023 development; 2024 final test, withheld (5 station-days; 2 usable after exclusions in H5)",
        independent_glaciers="station-days across several stations; usable test n is very small",
        verified_here=dict(structure=True, uncertainty=False),
        endpoints_available=["station broadband albedo (daily)"], ancestry=["2024 file header and some rows were seen during an import audit (docs/new_data_inventory.md)"]))

    # ---------------------------------------------------------------- not accessed
    for did, title, why in [
        ("DarkZone_2026_data", "Chevrollier et al. 2026 Science Advances data and inversion code", "paper full text and any archive were not reachable; SNICAR-fx repository holds the forward model only"),
        ("PRISMA_Bohn2022_scene", "PRISMA scene analysed by Bohn et al. 2022 (SW Greenland, 2020-08)", "not downloaded (large; no co-located field measurements per the paper)"),
        ("Himalaya_observations", "Himalayan spectra with algae and mineral measurements", "none found in this search and none in the repository"),
        ("Alps_algae_spectra", "Di Mauro et al. 2020 Alpine glacier algae spectra", "not retrieved"),
    ]:
        datasets.append(dict(id=did, title=title, status="NOT ACCESSED", reason=why, files=[], independent_glaciers=None))

    # ---------------------------------------------------------------- particle optics (not observations)
    lap = os.path.join(SCR, "ext_snicarfx", "src", "data", "light_absorbing_particles")
    datasets.append(dict(
        id="Particle_optics_libraries", title="Particle optical-property libraries (model inputs, NOT observations)",
        source="BioSNICAR (impurity npz files) and SNICAR-fx data/light_absorbing_particles", status="available as software data",
        biosnicar_impurities=["black_carbon bc_ChCB_rn40_dns1270", "snow_algae empirical Chevrollier2023", "glacier_algae ice_algae empirical Chevrollier2023",
                              "dust dust_greenland_Cook_CENTRAL_20190911"],
        snicarfx_files=sorted(os.listdir(lap)) if os.path.isdir(lap) else [],
        licence_note="BioSNICAR MIT; SNICAR-fx repository ACSL v1.4 (restricts commercial use) - reuse needs a licence decision",
        files=[], independent_glaciers=None))
    return dict(datasets=datasets), rows


def main():
    os.makedirs(OUT, exist_ok=True)
    body, rows = build()
    manifest = dict(
        generated_utc=dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        repo_commit=git_head(REPO), biosnicar_commit=git_head(BIOSNICAR), snowlaps_commit=git_head(SNOWLAPS),
        python=sys.version.split()[0], numpy=np.__version__, pandas=pd.__version__, **body)
    with open(os.path.join(OUT, "data_manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=1, default=str)
    cols = ["dataset", "sample_id", "site", "campaign", "sampling_group", "has_albedo", "has_hcrf", "has_counts", "role", "note"]
    with open(os.path.join(OUT, "sample_groups.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    missing = [(d["id"], f["path"]) for d in manifest["datasets"] for f in d.get("files", []) if not f["exists"]]
    print(f"{len(manifest['datasets'])} datasets, {len(rows)} sample rows; missing files: {missing or 'none'}")
    for d in manifest["datasets"]:
        if d.get("n_samples_total") is not None:
            print(f"  {d['id']:34s} n={d['n_samples_total']:4d}  glaciers={d.get('independent_glaciers')}")


if __name__ == "__main__":
    main()
