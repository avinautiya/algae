# Independent thermal-state information and Himalayan test feasibility

No Tower2 camera values were opened. No new melt prediction was scored. Source process features and small public Himalayan temperature files were audited.

## Pirámide cannot supply the proposed cold reference

`PIRAMIDE_COLD_SIGNAL_AUDIT.json` is reproduced from source-bound Tower1 thermal CSV only. In the672training samples, DT5 ranges−0.006 to+0.907°C; there are **zero** observations below−0.4°C. Across the entire4760sample process record its minimum is−0.353°C, still inside source±0.4°C accuracy. The record therefore supplies no independently resolved below-freezing ice-reference event. In training,320readings exceed+0.4°C. A single constant offset cannot make the whole trainingDT5 record equal0 within±0.4°C; the corresponding allowable offset interval is empty.

This does not prove warm glacier ice. Last-probe position/contact, sensor bias and local thermal gradients remain alternatives. The records are debris-probe temperatures; none measures an ice profile below the interface. With separately bounded±0.4°C reading errors, a two-probe temperature difference has±0.8°C uncertainty. Across the0.125m basal probe spacing that is±6.4K/m before depth uncertainty. A30minute temporal difference has worst-case±1.6°C/hour reading-error contribution. These are bounds, not standard deviations; stable offsets would cancel in temporal differences, but stability was not independently calibrated.

The source [Pirámidepaper](https://doi.org/10.1029/2025JF008360) provides nominal geometry/accuracy and lithology estimates, not a separately calibrated ice profile or absolute debris calorimetry. The existing k/C/R scaling ambiguity remains for a melting-point Robin boundary. Source [Rowan et al.2021](https://doi.org/10.1017/jog.2020.100) also documents an initially interfacial Khumbu thermistor migrating away from ice and warming above0°C. Positive nominal basal temperatures alone therefore cannot establish new contact physics.

## Actually downloaded and verified: NgozumpaNG2

[Zenodo3935685](https://doi.org/10.5281/zenodo.3935685) supplies a150,487byte native text record. Its published MD5 matches. `NG2_NATIVE_THERMAL.csv`, `NG2_THERMAL_ELIGIBILITY.json` and `audit_ng2_native.py` provide a reproducible binding:

-1940profiles,6December2014 at01:05 through3April2016 at19:05.
-Eleven probes at0.01/0.2/0.4/0.6/0.8/1.0/1.2/1.4/1.6/1.8/2.0m.
-Source states2.0m debris and installation of the deepest probe at the debris–ice interface; accuracy±0.25°C.
-Native six-hour steps contain±1second timestamp jitter; exact labels are preserved, with no invented UTC.
-The interface reaches−6.18°C. All measured temperatures in the lower1–2m remain below−1°C for approximately107days from10December2014 through27March2015.

This is a credible **cold thermal-response** record with signal well beyond stated sensor accuracy. Its lower1–2m subdomain avoids relying on the source paper's problematic near-surface probe. It does not contain a probe below the ice interface, independently measured absolute bulk heat capacity, basal heat flux or an ice-ablation target. No melt-accuracy claim follows from these temperatures alone.

The long cold interval was identified from inspected temperature features; it is development evidence, not a blind test chosen independently of observations. A chronological split or a different source must be frozen before any prediction comparison. Six-hour sampling supports seasonal/cold-wave behavior, not claims about resolved hourly melt pulses. Subgrid phase state, initial ice cold content and effective geometry remain assumptions.

## Imja: readable temperature data, unsuitable cold-reference shortcut

[Zenodo3906947](https://doi.org/10.5281/zenodo.3906947) provides a1.48MB temperature CSV, verified by MD5 and decoded to `IMJA_NATIVE_THERMAL.csv`. Kelvin units are explicitly converted to Celsius. Two undated trailing rows with temperature values are excluded without inventing dates; the exclusion is recorded.

The file contains8429dated half-hour rows through22November2014, beyond the publication's9November field end. AtIM14 the lowest source probe is0.24m while the paper lists0.26m debris. Through8November this probe never goes below−0.3°C (minimum+0.1°C); its large cold excursion occurs on9November night and later. Without instrument-retrieval timing, that cannot be used as a demonstrated cold-ice response. The separate106byte melt file remains unfetched/unopened. Its metadata describes a long18May–9November integration, not multiple independent chronological calibration/test intervals.

## Other primary sources

| Source | What is public | What prevents an immediate absolute-C/melt test |
|---|---|---|
| [NgozumpaNG1](https://doi.org/10.5281/zenodo.3935687) | Yearlong upper-debris temperatures | Deepest probe0.77m versus estimated2.2m cover; no ice-interface measurement. |
| [LirungChand](https://doi.org/10.5281/zenodo.4446387) | Verified50.3KB thermalCSV, three short seasonal campaigns; melt file remains unopened | Site/depth changes, roughly13day campaigns, no separate below-ice temperature profile or measured absoluteC identified. |
| [KhumbuPANGAEA883071](https://doi.org/10.1594/PANGAEA.883071) | Seven thermal site-year files and calibration/depth metadata | Probe migration is documented; period-specific interfaces must be audited. Metadata's claimed NepalGMT+4.5 conflicts with standard Nepal clock: preserve source labels until resolved, do not silently correct. |
| [ChangriNup3048780](https://doi.org/10.5281/zenodo.3048780) | Temperature profiles plus independent HFP01 interface heat flux,±3% stated | Already inspected elsewhere in this project; plate heat is not automatically latent melt. Sensors may migrate. Useful absolute heat-scale process calibration with honest prior-exposure labeling. |
| [Kennicott2025](https://doi.org/10.1002/gdj3.70032) | Thermal profiles and coincident stakes | Published bulkC is inferred from some melt observations. It cannot serve as independentC for those same observations. |

The Lirung2018 source is [Chand and Kayastha](https://doi.org/10.3724/SP.J.1226.2018.00357). The separate larger2021 Lirung study explicitly makes on-glacier AWS/debris records available through authors rather than an unrestricted archive; its public UAV/meteorology products do not automatically supply the missing paired thermal/ablation observations. No authors were contacted.

## Actual below-ice records: EverDrill

The [BAS2021V2record](https://doi.org/10.5285/32ecd5f4-1f00-4eeb-bd1e-5eaffb60f556) contains six Khumbu borehole strings,2017–2019, distributed as23.1MB ASCII files. Its metadata documents ice-bath calibration with expected uncertainty near0°C around±0.05°C, explicitly not a hard maximum. It flags cable failures and invalid Site4sensor subsets. No full collection was downloaded here.

The associated [Miles et al.2018 primary paper](https://doi.org/10.1038/s41598-018-34327-5) provides measured subsurface ice profiles and drilling-equilibration treatment. These are independent evidence for Himalayan cold ice but are not co-located/contemporaneous with NG2or the2014–16Khumbu debris records. Do not paste one glacier site's profile into another as a measured initial state. Seasonal cold-wave propagation has already been studied using these records, so that generic topic is prior art.

## Executable next physical test and exact claim ceiling

Use NG2's observed1.0m temperature as a boundary for a1.0m-thick lower-debris domain, with its observed initial1–2m profile. Compare ordinary diffusion with a cold-ice-coupled boundary and explicit uncertainty in the hidden ice initial/deep state. Fit only on a frozen training section; evaluate future interior and interface temperatures. Keep zero-temperature/fixed-boundary and strong statistical thermal baselines. Publish all predictions before evaluation errors. This supplies a direct thermal-response test with Himalayan relevance, not ice-melt validation.

Known ice effusivity can theoretically provide an absolute scale only after its initial/deep boundary problem is constrained. In a homogeneous idealization, e_d=C_d*sqrt(alpha_d); if thermal observations identify eta=e_i/e_d andalpha_d whilee_i is known, then C_d=e_i/[eta*sqrt(alpha_d)] andk_d=e_i*sqrt(alpha_d)/eta. These algebraic relations do not prove thateta is identifiable from finite noisy data. An arbitrary hidden ice initial-temperature profile may mimic the boundary response, and phase/stratification add nuisance parameters. Quantify that ambiguity rather than fixing an unmeasured ice profile to manufacture a result.

No inspected source currently supports a fresh, co-located complete package of calibrated ice profile, independently measured debrisC, native weather and multiple untouched ice-ablation intervals. The strongest immediate advance is the executable NG2thermal test plus a transparent identifiability study. A superior environmental melt predictor and historical mathematical originality remain unestablished.
