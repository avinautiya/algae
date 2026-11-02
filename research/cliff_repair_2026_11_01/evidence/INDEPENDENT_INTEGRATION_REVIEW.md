# Independent integration review of the v5 phase source and matched-ROI exporter

The reviewer did not edit the owner's solver, ROI code or outputs, and did not launch a flow campaign. `review_integration.py` independently inspects the three two-second physical smoke outputs and executes small synthetic invalid-bundle checks. Its JSON output records the exact source hashes reviewed. This is code/arithmetic review, not field validation or independent physical confirmation.

## Outcome

**The concrete schema/preflight defects identified in the first review have been repaired in the reviewed final source. No remaining blocker was found for releasing this bounded numerical prototype with the scope stated below. It is not ready to claim a full debris-covered cliff model, environmental accuracy or scientific novelty.**

The inspected final files are `roi/matched_roi.py`, `roi/run_v5_verification.py`, `phase/swe2d_v5.py` and its vendored `phase/swe2d_v4_pinned.py`. The three regenerated bundles under `roi/results/final_physical_smoke/` match all four current source hashes. Initial prototype bundles had an earlier record schema and stale solver hash; they cannot be relabelled as the final run. They should remain archived separately rather than overwritten.

## Bugs found and independently rechecked after repair

1. The initial loader could accept `geometric_melt_m3=NaN`: `abs(volume-NaN)>tolerance` is false. The final loader rejects nonfinite JSON constants and validates the volume's numeric domain.
2. The initial loader did not reconstruct the footprint from bounds, allowing mislabeled bounds with an otherwise correct mask. The final loader checks the bounds, grid alignment, physical area and exact reconstructed rectangular mask.
3. Initial imports could carry a forged field-evidence label or nonfinite heat/residual metadata. The final loader refuses those records and retains the conditional-numerical classification.
4. The initial comparator could accept NaN area. The final comparator applies the common strict record validator before comparing support.
5. The initial exporter wrote the NPZ before discovering a nonfinite configuration, leaving an orphan after an input refusal. The final exporter preflights configuration, metadata, arrays and record before publishing. Independent reproduction confirms that this invalid input now leaves neither output.
6. Initial run provenance omitted several constructor inputs and the inherited kernel hash. The final spec and bundle metadata include the explicit Q, water temperature and geometry overrides, physical duration, timestep/grid setting, and hashes of the v5 solver, pinned v4 dependency, exporter and runner. Omitted constructor defaults are recoverable from the authenticated pinned source; no claim of environmental truth follows.

All five forged-record/area examples in the independent check are refused. The nonfinite configuration is refused without an NPZ orphan. The reproduction is retained in `INDEPENDENT_INTEGRATION_REVIEW.json`.

## Independent arithmetic of actual outputs

The coarse, fine and fine/half-step runs each cover the same0.08m² fixed-XY ice-floor footprint and2s physical duration. Recomputed ROI geometry agrees exactly with its recorded volume. Geometry minus heat-implied ROI melt is about5.05e−19,−9.98e−19 and5.38e−19m³ respectively. Whole-interior geometric bed loss agrees with integrated phase ice loss to about1e−18m³. Generated liquid volume agrees with0.917×ice loss to about3e−23m³; latent energy agrees with917×334000×ice volume within1.5e−14J. All final bundle loads succeed.

These are consistency checks between the executed geometry and the same declared thermal source. They establish no independently measured heat-transfer closure or environmental melt accuracy. The two-second grid/timestep differences do not establish event-timing convergence, long-time splitting convergence, behavior near contact loss or a recovery boundary.

## Source and time bookkeeping

The v5 physical path advances the hydraulic/thermal kernel fordt and converts the same accumulated heat window into phase change over thatdt. Bare stationary temperate ice supplies water at the melting point with zero born horizontal momentum: depth increases byrho_i/rho_w×melt depth, whilehu,hv andhT remain unchanged during the source step. The ledger records phase mass, sensible heat, latent heat, positivity sources, boundary exports and explicit resets separately.

The covered-source refusal occurs before the physical hydraulic window, preventing the earlier half-applied covered step. Uniform redistribution is refused in physical mode. Reusing a consumed heat window is refused. The legacy quasi-steady hydraulic-settling/morphological API remains a different approximation; its pseudo-time heat is not a physical-time energy certificate. Any report should preserve that distinction.

The returned thermal/latent ledger is not a complete mechanical-energy/ALE ledger. Gravitational potential energy, moving-interface pressure work, vertical momentum, cold ice, pore storage and debris enthalpy remain outside the solved system. Positive phase mass closure must not be described as validation of those omitted equations.

## Remaining release cautions

- **Coverage:** the native physical path currently refuses phase change beneath solid/porous debris without a justified drainage/storage model. It cannot yet execute the project's main covered-ice mechanism as a physically verified fully coupled phase model. The explicit quasi-steady assumed transfer is not a substitute for that missing closure.
- **Two-file publication:** input preflight is atomic with respect to invalid metadata, and hard-link publication prevents overwriting existing targets. NPZ and JSON are still two separate filesystem objects. A process crash or I/O failure between them can leave an orphan NPZ. Treat the manifest as the commit marker; absent/invalid manifest means no complete result. Do not describe the pair as crash-atomic or automatically overwrite an orphan.
- **Observation time:** generic `observe` receives its time label from the caller. The final runner explicitly checks `sim.physical_time_s==END_S`, which makes these records consistent. Any future caller must bind capture/start/end and physical time just as explicitly; the geometry helper alone does not infer an observation window.
- **Authenticity:** source, array and configuration hashes establish declared byte integrity. A user could still fabricate an internally consistent bundle; no exporter can certify that it was run with those inputs or that the inputs represent nature.
- **Full-v4 gate:** these v5 smoke results do not supply the missing full-v4 legacy matched-ROI campaign/scorer. The evidence registry's v4 gate must remain blocked. A separately versioned v5 production result would need its own specification and evaluation.

No additional expensive campaign is justified by this review alone. The next scientific step remains the explicitly defined, prior-art-distinct recovery question with stronger initially matched baselines, after its physical model support is established.
