# CFD/Stefan/grain/exposure integration contract and status (inspected at base b566446; nothing newer on `claude/ice-network-v2` at inspection)

## Status of the existing pieces (file-level)
| piece | what exists | what it is not |
|---|---|---|
| `ice_network_v2/foam/interThermalFoam/interThermalFoam.C` | interFoam + one passive sensible-heat equation in both phases (`d(Cv T)/dt + div(T[Cw αφ + Ca(φ−αφ)]) − div(kα∇T) = 0`; header: "passive scalar, no buoyancy, no phase change") | no thermal buoyancy, no phase change, no ice, no grains |
| `ice_network_v2/foam/stefan1DFoam/` | 1-D ALE sharp-interface Stefan solver, density jump included, verified to first order against the exact solution (`stefan_of_results.json`); header comment still describes the second-order flux, the default is `fluxMode discrete` | 1-D only; no flow, no free surface, no coupling |
| `MILESTONE2_REPORT.md`, `VERIFICATION_RECORD_CORRECTIONS.md` | free-surface hydrostatic rest passes only on an orthogonal stepped bed (rest-state currents grow with refinement 2.0e-5→7.6e-5 m/s), curved beds fail, the conjugate-depth bore check is unresolved (depth error −33/−21/−12 %), outlet level conditions failed, frozen-flow diffusion is an isolated transport test; Stefan: discrete ledger vs physical sensible+latent residual reported separately | coupled free-surface flow + moving ice was not built |
| `ice_network_v2/predict.py` | depth-averaged reach model (1-D HLL shallow water, Dittus–Boelter, one interface power), restart-fingerprinted | not CFD; no grains |
Separately passing flow, heat and Stefan tests do **not** verify their coupling; no coupled test exists.

## Per-step exchange contract (extensive quantities, SI, same reporting surfaces for every ledger)
| interface | quantity | definition | consistency requirement |
|---|---|---|---|
| flow/heat → ice | P_w(t) [W], wetted area A_w [m²] | ∮_Γ k_w ∂_nT_w dA over water-contact patches | **same P_w removes sensible heat from water and melts ice** (`ExchangeRecord.water_energy_residual`) |
| ice → flow | recession m_n per patch, ice volume lost V_ice, meltwater volume (ρ_i/ρ_w)V_ice, v_Γ·n | m_n from Stefan balance (temperate: latent only) | `ExchangeRecord.validate` mass check; added water counted once (ALE flux **or** source, not both) |
| flow → grains | v, h, Q near the toe, per patch | time averages over the clearing time | inputs to the clast-scale force balance (`transport.py`) |
| grains → flow/ice | covered patch list/area; blockage | patches with z < z_f are removed from the water-contact set | unimplemented blockage/roughness raises `UnsupportedRun` |
| ice → grains | S_eng = c_mix A r_n | ice volume lost converts to grains via c_ice = c_mix/(1−c_mix) (`cfd_contract.volume_to_grain_inputs`) | frozen grains share u_ice; relative flux uses moving boundary |
| solver → ledgers | sensible_in/out, storage change, dissipation (once), heat_to_ice, latent = ρ_iL_fV_ice | per step | discrete residual and physical sensible+latent residual reported **separately** |
Flux export format: CSV `t0,t1,q_wall_W,wetted_area_m2,melt_volume_ice_m3,meltwater_volume_m3,sensible_in_J,sensible_out_J,storage_change_J[,dissipation_J]` parsed by `cfd_contract.parse_exchange_csv`; analytic fixture `neumann_fixture`.

## Grain treatment for clasts comparable to the water depth
Measured streams are 12–35 cm deep and carry clasts of 9–50 cm over an **ice** bed (Petersen 2024). A deep-bed transport law (Meyer-Peter–Müller, Shields bedload, Hjulström) assumes a sediment bed much deeper than the grain, flow depth ≫ grain size and bed shear on a mobile layer — none holds. Implemented instead (`transport.py`): (i) a **clast-scale force balance on smooth ice** (drag vs friction μ ≈ 0.05; reproduces Petersen's published 9–30 cm and 14–51 cm to < 1 cm in tests; buoyancy optional, lift not included, partial submergence interface only flagged); (ii) an **explicit clearing model** (mobile fraction × ζv/L_clear) whose ζ and L_clear are *assumed inputs with no default*. A resolved-clast CFD (discrete clasts in the flow) is the only way to remove assumption (ii) and is not attempted here.

## Blockers for a coupled run
1. 2-D/3-D interface-normal mesh motion with VOF free surface (custom dynamic-mesh BC, conservative remap) — not written.
2. Verified free-surface outlet/level boundary and rest-state parasitic currents — unresolved.
3. Thermal buoyancy/dissipation heating bookkeeping with a single count — not implemented.
4. Clast-flow interaction closure or resolved clasts — not implemented.
5. Multi-epoch cliff observations with co-located discharge/temperature for any survival validation — not in the pinned archives (see DECISIVE_COMPARISON.md).
