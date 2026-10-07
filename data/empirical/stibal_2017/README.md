# Stibal et al. (2017) — manual ingestion

Stibal M. et al. (2017), *Algae drive enhanced darkening of bare ice on the Greenland ice sheet*, Geophys. Res. Lett. 44:11463–11471, doi:10.1002/2017GL075958. Site S6, 17 June – 11 August 2014.

The data are only distributed as Wiley/AGU supporting information (`grl56634-sup-0002-2017gl075958_data_si.xlsx`). Every automated route returned HTTP 403, so the files have to be added by hand.

1. Download the supporting-information file(s) in a browser, or request them from the authors. Put them in this folder.
2. Copy `mapping.example.json` to `mapping.json` and edit it:
   - the file, sheet and column names;
   - the count unit: cells per mL of meltwater, or set `cells_scale`;
   - the spectrum layout and whether the spectra are `albedo` or `hcrf`.
3. Check the ingestion: `python phase2/stibal2017.py`. It prints the matched samples, or exactly which file, sheet or column is missing.
4. Rerun Phase 4. `empirical_data.field_samples()` adds the matched samples as dataset `s6_2014`, and the field validation reports them as a third dataset.

**Assumption to state if used.** For albedo spectra (cosine receptor), the HCRF/albedo factor k does not apply. k is then given a N(1, `albedo_k_sd`) prior, where `albedo_k_sd` represents instrument calibration. The default 0.02 is an assumption, not a measurement.

Only samples with **both** a count and a spectrum covering at least 450–880 nm are used. Nothing is guessed: if a column name is wrong, ingestion stops with the list of available columns.
