# Milestone-0 completion report — hydraulic control of ice-cliff survival (foundations)

Branch `claude/ice-cliff-survival-m0`, base `b566446a749303dd55b00998cb288a69d9b9ce23` (verified as the tip of `claude/ice-network-v2` at start; no other calculation was running; no existing worktree, checkpoint or frozen result touched).

## Implemented and verified (41 tests, `python -m ice_cliff_survival verify`, < 1 s)
* Solid-volume grain accounting: concentration bases (c_mix/c_ice), porosity conversion, mantle source with delivered/undelivered split, englacial release equal to the moving-boundary relative flux (tested against independent volume bookkeeping), exclusive source accounting (double-counted mantle in S_ext removed and reported).
* Conservative toe stock with actual (clipped) export recorded; independent budget reconstruction on random forcing; constant-forcing analytic limit; zero-stock admissibility.
* Exposure geometry: wedge height/area against an independent polygon (shoelace) calculation; critical covering volume; the same solid volume covers different fractions for different repose angles (a volume threshold is geometry-dependent).
* Interface kinematics: Stefan melt speed against the exact Neumann solution; latent energy equals the independently quadratured conductive flux (1e-6); density-jump mass balance (u_w·n = (ρ_i/ρ_w − 1) m_n) and its equal-density limit.
* Scalar falsifiers (§6–§7 of MODEL_FOUNDATIONS.md): uniqueness/stability, ds*/dT and dR*/dT against finite differences on 60 random admissible families; cumulative response ≥ 0 and the exact identity D∫ΔR = Δs + ∫ΔC on 12 ordered-temperature cases; sign reversal when C_s ≥ 0 is violated; clipping refused; accumulation hitting-time upper bound on random bounded forcing.
* Clast-scale transport: reproduces Petersen et al.'s published mobile clast sizes (half-sphere 9–30 cm, cube 14–51 cm) to < 1 cm; v² scaling; buoyancy direction; clearing capacity zero below threshold.
* Observation operators and ingestion: horizontal→normal recession on a planar face, irregular interval rates with uncertainty, real backwasting loader (13 series, 304 observations, monotone D), gauge discharge in the published range.
* Deterministic restart with config/forcing/closure fingerprints (restart refused on mismatch); environmental runs refused (`UnsupportedRun`).
## Data (pinned, MD5 and size verified; raw not committed)
Both archives downloaded and audited (`data_audit/DATA_AUDIT.md`, `AUDIT.json`, `download_manifest.json`). Findings: 304 non-blank stake observations (309 rows incl. 5 blank) from ≤ 6 independent cliffs in 13 site-years, no regularly sampled series; gauge discharges reproduce 0.061–0.210 m³ s⁻¹, one listed experiment (gauge 4, 17 Jul) has no data folder; **slope_min > slope_max in every row of ice_cliffs_stats.csv** (10,306 cliffs).
## Assumed (not derived or validated)
Planar face, fixed repose angle, pore water not ice contact, clearing parameters (ζ, L_clear), temperate ice, force-balance mobilisation (no lift, partial-submergence only flagged), all closure coefficients in demos (synthetic).
## Not implemented
Coupled flow-heat/moving-ice/grain model; debris blockage of the stream; thermal buoyancy; collapse; ponds; any environmental run; any survival validation; the held-out comparison (data gap). CFD status: `INTEGRATION_CONTRACT.md` — passive-heat free-surface solver and 1-D Stefan solver exist separately, curved beds/jump/outlet unresolved, no coupled test.
## Literature depth
Petersen et al. 2024 read in full for the requested sections; Evatt and Sato partially (open HTML); Buri 2016, Douglas & Lamb 2024, Röhl 2006, Verhaegen & Huybrechts 2026 **not obtained** (paywall/403), listed as unread in `PRIOR_ART.md`.
## Next bounded campaign
`DECISIVE_COMPARISON.md` (specification only). No computation was launched: the required coupled physics does not exist.
