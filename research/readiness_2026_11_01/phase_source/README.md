# Isolated density-consistent meltwater source

This operator repairs one specific source-step omission; it is not integrated into the repository SWE solver. Do not claim it fixes the full model until a new version is coupled and verified.

For stationary temperate ice melting vertically bydz, conservation gives `delta_h = rho_i/rho_w * dz`, `delta_zb=-dz`. With new water atTm and zero horizontal ice velocity, `hu`, `hv` and `hT` remain unchanged. Consequently water temperature and horizontal speed dilute, while free-surface elevation changes by `-(1-rho_i/rho_w)*dz`. The equal-density limit keeps the surface fixed in an isolated no-drainage cell.

The heat responsible for melting must ALREADY have been removed by the thermal solve. This source must not remove latent heat a second time. Apply bed lowering exactly once, not the old lower_bed followed by a second lowering here. Supply only actual control-volume melt; do not count ghost cells. Track generated volume and its eventual boundary export in the full ledger.

The caller must check bed-update-interval convergence. Applying a full morphological interval's source as one lump may create splitting errors even when the isolated step conserves mass. Prescribed drainage may justify another operator only with an explicit independently checked flow/energy contract. A globally small phase source does not prove negligible feedback in a nearly dry branch.

Unsupported: cold-ice enthalpy, freezing, moving-ice horizontal momentum, seepage, vertical momentum/pressure work, and complete ALE or SWE coupling. No field-validation or novelty claim is made. Eight tests establish source-step identities and refusal behavior, not environmental accuracy.
