# Integration handoff for Claude

## Start with the available code, not another package search

This branch contains both missing research packets under `research/`. Read `research/START_HERE.md`. The foundation delivery commit is `6ddbb449d82a0132b89bc8d766eb37cbea06b1ea`. An earlier delivery commit, `29a72d3`, was defective and replaced; do not cherry-pick it. The corrected foundation diff contains only added research files.

This module is an independent, dependency-free contribution. It does not require your CFD installation to run, and it does not change existing solvers, chemistry, job runners or frozen campaigns. Run the commands in README.md before integration. Do not replace the ongoing V2 milestone with a new large campaign.

## Actual source audit

Audited base: `405805a5ba359be9127befcae90aba043990c38c` on `claude/ice-network-v2`; corresponding files were read through the GitHub connector and their blob hashes checked against this commit.

| Existing file | Audited functions/behavior | Integration consequence |
|---|---|---|
| `ice_network/src/ice_network/geometry/storage.py`, blob `5c6e56448d07ab69c81e48044778000d1f948c63` | `volume_from_depth`, `depth_from_volume`, `area_from_depth` implement analytic power-law hypsometry. | Lake storage already exists. Add table geometry as an optional alternative, keeping the existing analytic path and its tests. |
| `ice_network/src/ice_network/evolution/network.py`, blob `b94865a27d5b481c0bc1bab47ab3f52c254bd055` | `reset_state` stores `V` and `S=VT`; `lake_stage` uses the analytic storage inverse; `rhs` transports volume and heat between nodes. | Convert `E = rho_w cp S` explicitly at an adapter boundary. Do not treat `S` as joules or introduce a second lake reservoir. |
| Same `network.py` | `edge_state` returns source discharge, downstream temperature, channel heat and added meltwater; `rhs` deposits meltwater downstream. | The lake port carries source discharge and source enthalpy. Channel melting is handled after that port, exactly once. |
| Same `network.py` | `budgets` audits extensive mass and energy; `run` uses LSODA inside a geometry split and exposes `atol_energy`. | Preserve its budget conventions; do not silently replace its integrator with this demonstration integrator. Compare adapters in controlled cases first. |
| `ice_network/src/ice_network/constants.py`, blob `70946ecc84e4b542364359ad908d9acf76a686ad` | Water density 1000 kg/m³, ice density 917 kg/m³, cp 4217 J/kg/K, latent heat 334000 J/kg. | This module matches those values. Pass constants explicitly if a CFD setup uses another property set. |

This is an interface-level source audit, not certification that the existing network is accurate in all hydraulic regimes. The new CFD solver was not available for inspection or adapter testing when this package was written.

## How to use a real CFD boundary

1. The lake supplies stage `H(V)` and an appropriate donor temperature. Whether stage is the admissible hydraulic boundary depends on characteristic directions and the actual inlet formulation.
2. Integrate the **actual numerical boundary flux** during the channel step: liquid volume and sensible enthalpy, using flow-weighted temperature. In an air–water solver, extract the liquid flux consistently with its phase-fraction/mass flux; do not integrate air volume as lake water.
3. Construct one `PortTransfer(volume_m3, sensible_energy_J)` for that exchange. Use opposite signs in the two domains, not independently interpolated duplicates. For backflow, use the channel donor temperature; the lake temperature must not replace it.
4. Apply external inflow and power to lake storage consistently, then apply the transfer using `exchange_port`. Reject/shorten the whole coupled step if mass or heat is unavailable. Do not clip the lake while allowing the channel to keep the rejected flux.
5. For an implicit partitioned solve, iterate stage/flux to a declared tolerance using restorable trial states. Roll back both domains between trial iterations; do not repeatedly accumulate the same transfer. Refine coupling time and iteration tolerance separately. A same-flux ledger alone does not prove partitioned-coupling stability.
6. Account for actual channel storage, channel meltwater and downstream export in the combined ledger. The demo's zero channel storage is an explicit approximation and must not leak into an unsteady CFD adapter.

`advance_lake` is already a complete lake update for its own monotone rating. **Do not call `exchange_port` on its returned lake state with its returned port again**: the exchange has already been applied. Use `exchange_port` for external-solver fluxes instead.

## Required integration checks before calling the coupled CFD model finished

* Closed coupled domain with no sources: mass and thermal energy remain constant within measured numerical error.
* Known-volume warm pulse across the boundary: the same volume/enthalpy leaves one domain and enters the other.
* Backflow: donor temperature switches with flow direction, without unexplained energy creation.
* Transient channel filling: lake loss equals channel storage gain plus export, with melting accounted for separately.
* Geometry/remeshing event: the combined budget includes conservative remap and phase-change mass exactly once.
* Lake-to-channel coupling time refinement: stage, heat export and hydrograph converge independently of spatial refinement.
* Appropriate hydraulic influence: a downstream change cannot be assumed to affect the lake through a supercritical reach; show the actual control/connection.

No integration into `Network` or the new CFD solver has been performed by this package. The tests verify this module and its synthetic driver only. Preserve this limitation until a real adapter passes these checks.

## Immediate continuation

1. Import the clean foundation commit and module commit without replacing existing files.
2. Run both foundation scripts and the module's unittest command.
3. Run the synthetic lake CLI once in a new output directory.
4. Continue your existing V2 baseline CLI and pinned CFD installation attempt. If installation fails, record the actual error and keep independent tasks moving.
5. Integrate table geometry and the flux interface only when the corresponding solver boundary is ready. Leave hydrofracture, debris and full ice mechanics out of this milestone.
