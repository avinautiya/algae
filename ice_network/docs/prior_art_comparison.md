# Prior art: the equations SaDS and SHED actually use, and the capability this project would add

Status of this document: a comparison of **published formulations**. Nothing here shows that our method is new or better. Items marked [not read] were not accessed.

## Sources read and what was and was not obtained

| Source | What was read | Access / pin | Not obtained |
|---|---|---|---|
| SaDS: Hill and Dow, JGR Earth Surf. (2021), doi 10.1029/2021JF006309 | **Chapter 4 of T. Hill's University of Waterloo MSc thesis, "Mathematical modelling of supraglacial meltwater production and drainage"**, which states it contains a modified version of this manuscript (UWSpace item 10012/17307; 19.5 MB PDF read in full for section 4.2; Chapter 3 numerical methods not read) | open thesis; paper abstract via Crossref (CC BY-NC 4.0) | the published paper full text (publisher returns 403); the thesis text may differ from the published paper |
| SaDS archive, Zenodo 10.5281/zenodo.4923858 (CC BY 4.0), `sads_results.zip` 7,937,218,316 bytes | **only the zip central directory (17 kB, by HTTP range request) and 5 README files**. 136 entries: 38 `.mat` model-output files (~0.1-0.5 GB each), 24 MATLAB `.m` plotting/analysis scripts, shapefiles, figures. **It contains no model solver source.** The README states the archive is "analysis scripts and model outputs" and refers to a GitHub repository (not accessible here; unauthenticated clone refused) | not downloaded (7.9 GB) | solver source; whether a separate code release exists was not established |
| SaDS-Greenland, Zenodo 10.5281/zenodo.7776531 (Hill and Dow, TC 2023) | record metadata only: one 3.35 GB `SaDS-Greenland.tar.gz`; a gzip stream cannot be listed without reading it | not downloaded | contents (may or may not contain the solver) |
| SHED v1.0: Gantayat et al., GMD 16, 5803 (2023) | full text (open), sections 4.2-4.3 and code/data availability | pinned: Zenodo 10.5281/zenodo.7633220 (CC BY 4.0, Fortran-77, 33 files, ~160 kB); `channel_inc1.f` md5 4fd00df8..., `lake_chn_inc_ov1.f` md5 5ca839dd... downloaded to an isolated directory and **only the header and control flow of `lake_chn_inc_ov1.f` were read** | line-by-line implementation of the incision update |
| Schoof, Cook, Kulessa, Thompson, J. Fluid Mech. 961, A4 (2023), doi 10.1017/jfm.2023.130 | open full text, abstract and model-section text searched for heat treatment | open access (CC BY) | derivation details beyond the abstract |
| Mejia and Gulley et al., TC 19, 6445 (2025) | abstract, data availability | open | methods in detail |
| Tero et al., Theory Biosci. 127, 2008 (Physarum) | **metadata only (Crossref)**; the mathematics below is quoted from memory of the general model and is therefore **[not read]** | - | full text |

## SaDS: exact equations (thesis section 4.2)
Two coupled systems on an unstructured triangular mesh: a distributed sheet on elements and channels on **every interior mesh edge** ("we treat every interior edge in the triangular mesh equally as a possible location for a supraglacial stream and let the channel network evolve naturally").

* Hydraulic potential, sheet and channel: phi_s = rho_w g z_s + rho_w g h_s (4.1); phi_c = rho_w g z_c + rho_w g (h_c - H_c) (4.6) (z_c lip elevation, h_c water depth, H_c incision depth).
* Sheet mass: dh_s/dt + div q_s = m_s + Xi_s/(rho_w L) (4.2); q_s = -k_s h_s^{alpha_s} |grad phi_s/(rho_w g)|^{beta_s - 2} grad phi_s/(rho_w g), alpha_s = 5/4, beta_s = 3/2 (4.3); Xi_s = |q_s . grad phi_s| (4.4).
* Channel mass: d sigma/dt + d q_c/ds = w_c m_c + Xi_c/(rho_w L) (4.7); channel melt source Xi_c = q_c d phi_c/ds (4.8), "applied along the entire channel perimeter, rather than the wetted perimeter only".
* Channel flux: q_c = -k_c w_c h_c^{alpha_c} |(1/rho_w g) d phi_c/ds|^{beta_c - 2} (1/rho_w g) d phi_c/ds (4.11), alpha_c = 5/3, beta_c = 3/2, k_c = 10.
* Channel growth: dA_c/dt = Xi_c/(rho_i L) + (rho_w/rho_i) w_c (m_c - m_s) (4.9); rectangular section with fixed width-to-depth ratio, w_c = r H_c, r = 3; dH_c/dt = (1/2)(rho_w/rho_i)(m_c - m_s) + (1/(2 w_c)) Xi_c/(rho_i L) (4.13). Lip lowering dz_c/dt = -(rho_w/rho_i) m_s (4.10).
* Exchange sheet to channel: fraction f_e of the normal sheet flux, piecewise linear in channel fill (4.16).
* Prescribed channel melt m_c "set to zero in the current work": **no water-temperature variable, no advected sensible heat**. Dissipation is assumed to melt ice immediately (ice at the melting point).
* Numerics: explicit RK4, dt 30-240 s, finite volume, mass-conserving at nodes, moulins and free-flux nodes remove water.

## SHED: exact equations (GMD paper section 4.3)
* Routing between cells: travel time from Manning's equation (bare ice) or Colbeck's equation (snow/firn); cells filled to lake storage; **channels are not evolving graph edges**: a lake overflows through a **single outlet** at the lowest cell of its catchment boundary (after Kingslake et al. 2015).
* Lake depth and outlet incision (equations 6-9 of the paper; the HTML text extraction scrambled the fractions, so only their structure is quoted and **the exact prefactors are not transcribed**): dH_L/dt is a lake-geometry factor times (Q_in - beta zeta^{1.5}), with zeta = H_L - H_c the water level above the channel bed (a broad-crested-weir / Darcy-Weisbach stage-discharge law, channel width fixed at w = 5 m, roughness f_R = 0.25, bed slope S from the DEM); the channel bed height above the lake bottom falls as dH_C/dt = -alpha zeta^{1.5}, with alpha built from rho, f_R, L (latent heat), rho_i, g and S, i.e. **frictional (dissipative) melt only**. DEM cells along the channel are then lowered.
* No water temperature; no two-outlet competition; no hydraulics with backwater (stage-discharge law of a weir with a fixed outlet).

## Schoof et al. 2023 (JFM)
St Venant equations in one channel fed by a lake, coupled to channel elevation evolving by advection (ice flow), uplift and downward melt, "dissipation-driven"; after reduction a stream-power hyperbolic model with a shock; a single channel; no advected sensible heat.

## Who already provides what

| Capability | SaDS | SHED | Schoof 2023 | This project (version 1 target) |
|---|---|---|---|---|
| Channels that grow by melt and change routing | **yes**, all mesh edges candidates, no prior channel locations | partial: single outlet, bed lowering, DEM updated | single channel | yes, on a declared candidate graph (cannot create paths outside it) |
| Water-level (potential) dependent flow, ponding | **yes** (potential includes depth) | weir law for the lake outlet | **yes** | yes (quasi-steady GVF with level boundaries) |
| Mass-conserving network | **yes** | not stated as exact | n/a | yes, by construction + audit |
| Dissipation-driven incision | **yes** | **yes** | **yes** | yes, as an explicit heat source |
| **Advected sensible heat from warm water as a melt source** | no (m_c prescribed, set to 0) | no | no | **yes (the difference to test)** |
| Explicit wetted-perimeter partition of melt between bed and walls | no (perimeter-wide, fixed aspect ratio) | no | wetted-perimeter relations | yes (equal heat flux on the wetted perimeter) |
| Competing outlets with level coupling in a graph | implicitly (mesh) | no | no | yes (small graphs) |
| Surface lowering narrowing channels, ice flow | yes (m_s) / n/a | n/a | yes | **no** in version 1 |
| Hydraulic jumps, momentum | no | no | St Venant (yes) | no |
| Validation against independent discharge | synthetic SHMIP and a Greenland catchment (moulin inputs) | proglacial discharge at a gauge, lake area | none (theory) | none yet |

## The exact capability this project would add, stated as hypotheses
H-A. **Thermal forcing by advected water temperature** as a physical incision mechanism in a graph of storage nodes and open-channel edges, with an energy budget that separates transported sensible heat, latent heat and dissipative heating. SaDS set this term to a prescribed constant (zero). Other thermal-erosion models for lake drainage or englacial channels may exist and **were not searched or read here** (e.g. work on warm-water erosion of ice spillways, and englacial channel evolution models); H-A is therefore only "not in the three formulations read".
H-B. A **small, conservative graph reduction** that reproduces the routing of a resolved reference (SaDS-like mesh solution or a 2-D hydraulic solution) at far lower cost. This requires such a reference; none is available locally.
H-C. Route-switching prediction between **competing outlets**. SaDS can in principle produce such switches on a mesh; SHED and Schoof cannot.
None of these has been demonstrated. In particular, the claims "dynamic routing", "channel incision" and "network modelling" are **not new**.

## Where the biological analogy enters, and what it is not
Borrowed element **[not read; general structure of the Physarum model]**: an edge conductivity D_ij adapting to flow, dD/dt = f(|Q_ij|) - gamma D, with Kirchhoff conservation at nodes. The Roethlisberger-type incision law (growth proportional to dissipation, Xi = Q dphi/ds) is a physically derived adaptation law of the same structure, known long before the biological model; SaDS's (4.7)-(4.13) is an instance. No decay term gamma D, Murray's law or cost-minimisation is used or assumed here; channels are not pruned (they become hydraulically inactive and keep their geometry in this version). Therefore "bio-inspired" is a conceptual origin only, and the model, if built, is an established physical drainage network. The newly implemented parts would be the graph hydraulics, the thermal transport, and the audit, not the adaptation idea.

## Capability absent from each implementation, and the test that would demonstrate it
Basis: the published equations quoted above (SaDS thesis chapter 4; SHED paper section 4.3). **Neither SaDS nor SHED was run**; the SaDS solver source was not obtained (the Zenodo archive holds outputs and plotting scripts) and SHED's Fortran was only skimmed. "Absent" means "not in the formulation as published". A test is "demonstrating" only if it runs the existing method (or a faithful re-implementation of its stated equations) next to ours on the same case; none of these has been done.

| Existing method | Capability absent as published | Test that would demonstrate the difference | Status |
|---|---|---|---|
| SaDS | Advected sensible heat of warm water as a melt source (channel melt m_c prescribed, set to 0; ice at the melting point; dissipation is the only melt) | Warm-lake drainage in a flat channel where dissipation is negligible: incision rate versus inflow temperature T_in = 0, 0.1, 0.3 C. Dissipation-only (SaDS-equation, eq. 4.13) predicts rate independent of T_in; the energy balance predicts rate proportional to h_T T. Requires a SaDS run or a re-implementation of eqs 4.7-4.13 for one channel | not run; our T_in = 0 limit and the heat-driven rate are in the tests of this milestone, the SaDS comparator is not |
| SaDS | Channel melt applied on the **wetted** perimeter with bed/wall partition from the energy balance (SaDS applies dissipation "along the entire channel perimeter", fixed width-to-depth ratio r = 3) | Constant discharge through a narrow-deep and a wide-shallow section: compare predicted bed lowering versus widening against a resolved (CFD or laboratory) melting channel | not run |
| SaDS | Critical-depth (weir/brink) control and entrance loss: its flux law (4.11) is a potential-gradient power law | Reservoir draining through a short steep channel: analytic h(t) (the weir law); SaDS power law gives a different exponent | our analytic test passes; SaDS side not run |
| SaDS | A thermal/energy budget with separate sensible, latent and dissipative terms | Report of the three terms for the same case | ours only |
| SHED | **More than one outlet per lake and level-dependent splitting**: a lake overflows through a single outlet at the lowest cell of its catchment boundary | Lake with two sills 0.05 m apart; SHED selects the lower outlet only, while level-coupled hydraulics sends flow through both once the stage passes the second sill. Our branching and symmetry tests give the level-coupled answer; SHED not run | partial (ours only) |
| SHED | Advected sensible heat (incision is frictional only) and channel width evolution (fixed 5 m) | Same warm-lake test as for SaDS | not run |
| SHED | Water or energy budget diagnostics for the incision module | report | ours only |
| Both | Independent discharge validation of the channel evolution | No dataset with independent discharge and repeat channel geometry has been identified; Qiyi discharge is estimated from geometry and modelled melt, not measured | open |

The capability that would still need a decisive demonstration is the first row: heat-driven incision. If a SaDS run (or reimplementation) with m_c set to a nonzero prescribed value reproduces the same sill lowering as our advected-heat model on the warm-lake case, the difference collapses to "SaDS has a hook for a prescribed melt term" and the contribution reduces to deriving that term from a water-temperature balance.


### Correction to the SaDS dissipation-only comparison (from the budget table, `records/budget_tables.json`, case C)
SaDS assumes all dissipated energy melts ice immediately (Xi_c/(rho_i L), eq. 4.9). In this code, with inflow at 0 C, only the channel share of the dissipation is available to the ice (8.0 MJ of the 43.5 MJ total in the 6 h case; the rest is lost at the exit and heats the outlet water), and the finite heat-transfer coefficient lets the reach carry part of that heat away as a small temperature rise: **1.4 MJ of latent heat (4.2 kg of ice) versus the 24 kg a SaDS-type "all channel dissipation melts ice" assumption would give for the same 8.0 MJ**. So the T_in = 0 limit of this model does **not** reproduce the SaDS rate unless the thermal e-folding length rho_w c_w Q/(h_T P) is short compared with the reach length. This is a structural difference to be tested, not an agreement.
