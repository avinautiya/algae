# Transfer test: sites, measurements, data history (frozen before any extension result exists)

Nothing in this file reports a result of the chemistry-variable extension (B), which has not been built.

## What each candidate dataset contains and what has already touched it

| Dataset | Independent glaciers | Component measurements | Spectral quantity | Prior influence on this repo | Role now |
|---|---|---|---|---|---|
| S6 2017 (Cook 2020 archive) | 1 (S6) | algal cells, 47 plots, colocated | hemispherical albedo (ALDN/ALUP, cosine receptor) + nadir HCRF | k prior, sigma and radius prior, H1 target, earlier held-out tests | development (leave-one-day-out) |
| S Greenland 2021 (Chevrollier 2023 Table S4 + ASD) | 1 | algal cells, 18 plots, 2 days | HCRF | `heldout_v2` reverse fold; BioSNICAR default ice-algal optics were measured on algae from this campaign | **reused**, transfer check only |
| Dark Zone archive (Zenodo 18826013) | 2 sites (S21 2021, S22A 2022) | **none per spectrum** (paper: dust n = 10, mineralogy n = 17, cryoconite optics from 2 spectra at S21) | HCRF | 20 of the 105 spectra are identical to the S Greenland 2021 set (checked, max difference 5e-6); S21 supplies the algae and cryoconite optics used by the competitor | reconstruction/consistency only; no concentration validation |
| PROMBIO 2024 station-days | several stations | ice-algal cells (few usable station-days) | none at spectral level in hand | header and some rows were seen in an earlier import audit | **withheld final test**; not touched here |
| PROMICE KAN_L / KAN_M | 2 | none | radiometer broadband albedo, energy balance | forcing and endpoint of frozen H2/H4 protocol | only dataset with absorbed-energy/melt-relevant observations; no colocated spectra or components |
| Norway snow 2023 | 0 (snow) | none | HCRF | reproduction of snowlaps only | not a glacier-ice dataset |
| Himalaya / Alps spectra | none in hand | none | none | none | **no data**; no Himalayan claim possible |

## What an extension transfer test would need, and what exists
1. **A site not used for optics, priors or earlier tuning**, with spectra *and* component measurements. In hand: none. S Greenland 2021 and S21 are tied to the optics origin; S6 is the development site; PROMBIO has no spectra here.
2. **Component truth for more than the algal term.** In hand: algal counts only (S6, S Greenland 2021). No per-plot mineral, cryoconite or pigment measurement exists, so a test of mineral or cryoconite variability against truth is not possible.
3. **A quantity matching the endpoint**: albedo for albedo (S6 only), HCRF for HCRF (with k free), absorbed energy needs radiometer albedo and incoming flux (PROMICE, no spectra).

## Consequence
* The only honest transfer-style comparison available is **A versus B on S Greenland 2021 / S21 spectra, reconstruction of held-out spectral regions and consistency of retrieved algal amount with 18 counts**, labelled retrospective and reused. It cannot show that B generalises to new glaciers.
* Untouched observations reserved: PROMBIO 2024, `records/heldout_v2` (immutable). They stay reserved until B is specified, frozen and shown on development data to be worth testing.
