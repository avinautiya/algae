"""
Manual ingestion of Stibal et al. (2017) field data (S6, 17 June - 11 Aug 2014) for the field validation.

Stibal M. et al. (2017) Algae drive enhanced darkening of bare ice on the Greenland ice sheet.
Geophys. Res. Lett. 44:11463-11471, doi:10.1002/2017GL075958.

The supporting-information dataset (grl56634-sup-0002-2017gl075958_data_si.xlsx) is served only by
Wiley/AGU, which blocks automated downloads (HTTP 403). To use it:

  1. Download it in a browser from the article's "Supporting Information" section (or ask the authors),
     plus any spectra file, and put the files in data/empirical/stibal_2017/.
  2. Copy data/empirical/stibal_2017/mapping.example.json to mapping.json and edit the file, sheet and
     column names to match what you downloaded (open the xlsx and look - nothing here guesses them).
  3. Check:   python phase2/stibal2017.py          (prints what was read, or exactly what is missing)
  4. Run Phase 4 as usual: field_samples() picks the dataset up as 's6_2014' automatically.

A sample enters the validation only if it has BOTH a cell count and a spectrum (350-2500 nm, or at
least 400-900 nm so the four Sentinel-2 bands are covered). Counts must be cells per mL of meltwater;
set "cells_scale" to convert other units. Spectra may be albedo (cosine receptor) or HCRF; for albedo
the reflectance factor k (HCRF/albedo) is not applicable, and field_validation uses k ~ N(1, sd_k)
with sd_k from mapping.json ("albedo_k_sd"; a stated instrument-calibration assumption).
"""

from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DIR = os.path.abspath(os.path.join(HERE, "..", "data", "empirical", "stibal_2017"))
REQUIRED = {"counts": ("file", "sample_col", "cells_col"), "spectra": ("file", "orientation", "quantity")}


class StibalDataError(RuntimeError):
    pass


def available(directory: str = DIR) -> bool:
    return os.path.isfile(os.path.join(directory, "mapping.json"))


def _read_table(directory, spec, what):
    path = os.path.join(directory, spec["file"])
    if not os.path.isfile(path):
        raise StibalDataError(f"{what}: file not found: {path}")
    if path.lower().endswith((".xlsx", ".xls")):
        try:
            return pd.read_excel(path, sheet_name=spec.get("sheet", 0), header=spec.get("header_row", 0))
        except ValueError as e:                        # wrong sheet name: list the real ones
            sheets = pd.ExcelFile(path).sheet_names
            raise StibalDataError(f"{what}: {e}. Sheets in {spec['file']}: {sheets}") from None
    return pd.read_csv(path, header=spec.get("header_row", 0))


def _need(df, col, what):
    if col not in df.columns:
        raise StibalDataError(f"{what}: column {col!r} not found. Columns are: {list(df.columns)[:40]}")


def load(directory: str = DIR):
    """Returns (table, spectra) in the schema of empirical_data.field_samples():
    table columns dataset, sample, cells, cells_counted, sza, quantity; spectra = DataFrame with
    wavelength_nm and one column per sample."""
    from empirical_data import _solar_mu
    with open(os.path.join(directory, "mapping.json")) as fh:
        mp = json.load(fh)
    for block, keys in REQUIRED.items():
        for k in keys:
            if k not in mp.get(block, {}):
                raise StibalDataError(f"mapping.json: '{block}.{k}' is required")
    c = mp["counts"]
    cnt = _read_table(directory, c, "counts")
    for col in (c["sample_col"], c["cells_col"]):
        _need(cnt, col, "counts")
    cnt = cnt[[x for x in (c["sample_col"], c["cells_col"], c.get("datetime_col")) if x]].dropna(
        subset=[c["sample_col"], c["cells_col"]])
    cnt["sample"] = cnt[c["sample_col"]].astype(str).str.strip()
    cnt["cells"] = pd.to_numeric(cnt[c["cells_col"]], errors="coerce") * float(c.get("cells_scale", 1.0))
    cnt = cnt[np.isfinite(cnt.cells)]

    s = mp["spectra"]
    sp = _read_table(directory, s, "spectra")
    if s["orientation"] == "wavelength_rows":            # one row per wavelength, one column per sample
        wcol = s.get("wavelength_col", sp.columns[0])
        _need(sp, wcol, "spectra")
        spectra = sp.rename(columns={wcol: "wavelength_nm"})
        spectra.columns = ["wavelength_nm"] + [str(x).strip() for x in spectra.columns[1:]]
    elif s["orientation"] == "wavelength_columns":       # one row per sample, one column per wavelength
        scol = s.get("sample_col")
        _need(sp, scol, "spectra")
        wl_cols = [x for x in sp.columns if x != scol and _is_number(x)]
        spectra = sp.set_index(sp[scol].astype(str).str.strip())[wl_cols].T
        spectra.insert(0, "wavelength_nm", [float(x) for x in wl_cols])
        spectra = spectra.reset_index(drop=True)
    else:
        raise StibalDataError("spectra.orientation must be 'wavelength_rows' or 'wavelength_columns'")
    spectra["wavelength_nm"] = pd.to_numeric(spectra.wavelength_nm, errors="coerce")
    spectra = spectra.dropna(subset=["wavelength_nm"])
    for col in spectra.columns[1:]:
        spectra[col] = pd.to_numeric(spectra[col], errors="coerce") * float(s.get("scale", 1.0))
    if spectra.wavelength_nm.min() > 450 or spectra.wavelength_nm.max() < 880:
        raise StibalDataError("spectra must cover at least 450-880 nm (Sentinel-2 B2-B8)")

    site = mp.get("site", {"lat": 67.08, "lon": -49.4})
    matched = cnt[cnt["sample"].isin(spectra.columns)]
    if matched.empty:
        raise StibalDataError(f"no sample IDs match between counts ({cnt['sample'].head(5).tolist()} ...) and "
                              f"spectra ({list(spectra.columns[1:6])} ...); set the ID columns or rename")
    dtc = c.get("datetime_col")
    if dtc:
        t = pd.to_datetime(matched[dtc], format=c.get("datetime_format"), errors="coerce")
        if not c.get("datetime_has_time", False):        # date only: solar noon
            t = t.dt.normalize() + pd.to_timedelta(12 - site["lon"] / 15.0, unit="h")
        mu, _ = _solar_mu(t, site["lat"], site["lon"])
        sza = np.degrees(np.arccos(np.clip(mu, -1, 1)))
    else:
        raise StibalDataError("counts.datetime_col is required to compute the solar zenith angle")
    table = pd.DataFrame(dict(dataset="s6_2014", sample=matched["sample"].to_numpy(), cells=matched.cells.to_numpy(),
                              cells_counted=np.nan, sza=sza, quantity=s["quantity"]))
    keep = ["wavelength_nm"] + table["sample"].tolist()
    return table, spectra[keep], dict(albedo_k_sd=float(mp.get("albedo_k_sd", np.nan)))


def _is_number(x):
    try:
        float(x)
        return True
    except (TypeError, ValueError):
        return False


if __name__ == "__main__":
    d = sys.argv[1] if len(sys.argv) > 1 else DIR
    if not available(d):
        print(f"No {os.path.join(d, 'mapping.json')} - see this module's docstring and {d}/README.md")
        sys.exit(1)
    try:
        tab, spec, meta = load(d)
    except StibalDataError as e:
        print(f"Stibal 2017 ingestion failed: {e}")
        sys.exit(2)
    print(f"{len(tab)} samples with count + spectrum; quantity {tab.quantity.iloc[0]}; "
          f"SZA {tab.sza.min():.1f}-{tab.sza.max():.1f} deg; cells {tab.cells.min():.0f}-{tab.cells.max():.0f} mL^-1")
