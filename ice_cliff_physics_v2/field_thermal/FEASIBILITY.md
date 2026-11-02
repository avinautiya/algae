# Public debris-temperature component validation: feasibility audit

No predictive model was fitted or scored in this audit. All measurements,
coordinates and support statements come from the verified source files or the
explicitly cited author source. This does not validate the new wet-cliff
mechanism or whole-glacier melt prediction.

## Source verification

The byte-pinned Zenodo14625691 ZIP is9,948,630bytes, MD5
f740fc1d0e6ff19c1b8049b427db0c3a. Selected extraction preserves34 CSV/text members
with original ZIP CRC32, uncompressed size and SHA256. Only thermal profiles,
depths/metadata and mass-balance/debris-thickness tables were extracted. Safe
reproduction is in extract_source.py.

Existing bounded-range metadata for the later16389806 archive were checked:
all34 selected members have the same CRC32 and uncompressed size. This is subset
identity evidence, not a full checksum verification of the190MB latest archive.
It does not supply earlier2023 records that were missing in the old archive.

## Temperature support

There are12 profile files across2020–2023. All records have an explicit UTC
Datetime and an ADT companion exactly8hours earlier. Each file is a strictly
increasing uninterrupted5-minute series with no missing/nonfinite numeric temperature
columns. Header depths match separate depth CSVs exactly. Source instrument
accuracy is±0.2°C. Temperatures and geometry are observed; the depth-relative
measurements are not direct surface heat flux or water-contact measurements.

| Profile | Actual start/end UTC | Span days | Columns at depths cm |
|---|---|---:|---|
|2020a|27Jul08:00–6Sep20:25|41.517|4,11,18,25|
|2020b|27Jul08:00–6Sep20:25|41.517|5,15,22*|
|2021a|6Jun20:00–20Jul20:00|44|4,9,14,19|
|2021b|18Jul20:00–11Sep20:00|55|3,15,27,39|
|2021c|6Jun20:00–12Jul20:00|36|5,10,15,20|
|2021d|6Jun20:00–4Aug20:00|59|7,12,17,22|
|2022a|31May08:00–12Sep08:00|104|3,7,11,15|
|2022b|31May08:00–24Aug08:00|85|6,9,12,15|
|2022c|31May08:00–12Sep08:00|104|4,8,12,16,20,24,28,32|
|2023a|10Jun20:00–12Sep20:00|94|2,7,12,17,22|
|2023b|10Jun20:00–12Sep20:00|94|9,16.5,24,31.5|
|2023c|10Jun20:00–12Sep20:00|94|2,9.5,17,24.5|

The2021d source ends4August, although the summary describes10September. Raw time
support must override that summary date. Some profiles include occasional
negative temperatures; these are retained. A linear thermal-transfer comparison
can report them with a phase-physics caveat, but the current temperate wet-debris
kernel must not be represented as a freeze/thaw solver.

### Critical synthetic-boundary correction

*2020b has only TWO actual thermistors. Its22cm column is exactly0 throughout and
is an author-inserted basal boundary, not a third sensor. The author repository's
profile_info.csv lists two instrument IDs at5/15cm and Base at22cm. Its processing
class explicitly constructs Base as a zero series and describes it as an
unsensed interface assumption. The primary2022 article independently reportsN2.

Therefore the primary measured-endpoint, withheld-interior comparison can use
11 profiles and27 strictly interior probe series. Exclude2020b from that primary
benchmark: two real endpoints leave no measured interior. It can appear only in
an explicitly assumed-basal-boundary secondary experiment. Do not count its zero
column as field truth or as a successful prediction.

## Primary temperature comparison that is now feasible

Before fitting, fix a measured-subcolumn experiment: shallowest and deepest REAL
probes provide time-varying boundary data; strict interior probes are targets.
First14days may calibrate an effective diffusivity and the next14days test it.
The final frozen protocol specifies initial state from the FIRST training record
and continuous integration into test time without an interior-state reset.
All11 eligible files have sufficient contiguous support. Keep all qualified
profiles and strong identical-input baselines. Do not choose periods or exclude
profiles using predictive errors. The independent protocol agent owns the final
freeze and scoring specification; this paragraph is a feasibility assessment.

This validates a forced thermal-transfer component between measured depths,
not atmospheric forcing, the unsensed topmost debris, the basal ice boundary or
melt flux. Observed future endpoint temperatures make this a retrospective
conditional hindcast. Probes and successive timestamps are correlated; thousands
of5-minute predictions are not thousands of independent field sites.

## Real ablation links and limits

The summary explicitly associates2020b with P_AWS_1 and2023a/b/c with
P11/P12/P13. These are author-declared links, not a nearest-neighbor search.
The raw point_mass_balance.csv provides:

- P_AWS_1:25July–6September2020,−360kg/m²,21cm stake debris.
- P11:8June–2September2023,−1035kg/m²,24cm stake debris.
- P12:same2023 interval,−873kg/m²,29cm stake debris.
- P13:same2023 interval,−688kg/m²,25.5cm stake debris.

Profile thicknesses are respectively22/24/32/26cm. In particular P12's3cm
thickness difference is material; do not choose the more favorable thickness
using held-out melt. Initial coordinate distances are about2.83/0.33/0.11/0.56m
respectively, after normalizing source stake longitude from positive degreesWEST
to signed decimal longitude. Source stake dates have no established clock/time of
day here. Raw point values and rounded summer/annual summaries are not independent
observations and must not be double counted.

### Original2020 raw data partially resolve a forcing gap

The primary2022 paper links a public author repository. Five small files were
retrieved at pinned commit4863e0bec8e496452421f97d616b03925451f215 and retained in
original_2022/. Two real2020b sensors cover18June–9SeptemberUTC continuously at
5minutes. Original headers explicitly stateGMT−08:00; adding8hours exactly matches
ALL11,958 overlapping processed values at both5cm and15cm.

Thus2020 ablation's25July–6September temperature interval is available from the
original records. Exact probe-depth transition is not established: source June
configuration says0.5/18cm, while July/August configuration and the processed
27July onward series say5/15cm. For the leading~2days, a depth-separation range
0.10–0.175m and visit-time scenarios could be retained explicitly. The rest of
the processed window has confirmed5/15cm source geometry. Do not silently infer
the transition date from a favorable melt fit.

For2023 the measured profile still starts10June, later than the8June mass-balance
start date. That leading interval remains genuinely unobserved. Do not compare a
partial temperature integral against the full-season observed melt as if support
were complete. A bounded leading-input scenario could be secondary evidence;
without a defensible bound, the supported melt cannot provide a full-period
point-accuracy claim. The fixed-time temperature component is unaffected.

## Identifiability and calibration

Temperature transfer under pure diffusion constrains effectivealpha=k/C. Scaling
k and volumetric heat capacityC together preserves temperature while changing
conductive heat flux. Hence temperature skill alone does not establish melt-flux
accuracy. The source paper uses ablation to constrain this scale; its published
full-periodk values must not be reused as independent parameters on the same
observations. The newer summary's other year-specifick values likewise have no
independent calibration provenance established by this audit.

A possible secondary study would calibrate scale only on2020 ablation plus
source-supported temperature/depth scenarios, then evaluate2023 melt without
refitting to those three melt values. It would still need the2023 leading-time,
thickness, nonconductive-flux and source-target matching sensitivities. No such
calibration or scoring was performed here. Pure conduction, debris convection
and thermal memory are established literature topics; validation of this
component is not evidence of new wet-cliff physics or historical originality.

Primary sources: dataset metadata and raw tables; Petersen etal.2022,
https://doi.org/10.1029/2021JF006520 ; author data/config,
https://github.com/ericivanpetersen/supraglacial_debris_heat_flux_jgr/tree/4863e0bec8e496452421f97d616b03925451f215 .
