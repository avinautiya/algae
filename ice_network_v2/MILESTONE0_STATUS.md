# V2 milestone — status after first pass (blocker-first)

## Protection
* Isolated worktree `.worktrees/ice-network-v2`, branch `claude/ice-network-v2` created from `claude/ice-network-m1` (b891d84). Other worktrees untouched: `ice-network` (claude/ice-network-m1), `inversion-m1`, main checkout. Frozen tags/results (`ice-eval-v1-frozen`, `ice-capture-v1`, sequence/event records) unmodified. No background jobs were running at start (process scan clean). Resources: 4 cores, 15 GB RAM, ~14 GB free disk.

## BLOCKER 1 — the V2 packet is not available to this session
Searched the whole filesystem and every remote branch for `PROJECT_V2.md`, `MATHEMATICAL_UPGRADES.md`, `CFD_AND_VALIDATION.md`, `REFERENCES_AND_DATA.md`, `prediction_basis.py`, `verify_foundations.py`: none exist. The only uploads in the session are the four PDFs from the previous milestone. Consequently item 2 (run `prediction_basis.py` and `verify_foundations.py` unchanged, read the four documents, source/novelty comparison table) and the Ogier data-archive links of item 8 (listed in REFERENCES_AND_DATA.md) **cannot be done**. I have not reconstructed or re-written the packet from memory. Action needed: attach the packet (or commit it to a branch).

## BLOCKER 2 — no CFD software in the environment
`openfoam`, FEniCS, Firedrake, Basilisk, Nek, Dedalus, FiPy: none installed (pip/apt scan; `apt-cache policy openfoam` empty). Outbound HTTPS exists through the proxy (CONNECT to dl.openfoam.com and zenodo.org answered), but installation, version pinning and a module-compatibility demonstration (VOF interface capturing + conjugate/enthalpy energy equation + phase-change mass flux + moving boundary) have **not** been attempted yet, so no CFD claim exists. Not substituted by the depth-averaged solver.

## Item 1 — audit of what exists (file/function references from inspection; to be completed)
| Piece | Network model (`ice_network/src/ice_network`) | Resolved reference (`ice_network_reference/solver.py`) |
|---|---|---|
| Hydraulics | quasi-steady GVF with critical control: `hydraulics/channel.py` `gvf_upstream` (77), `edge_flow_forward` (122), `edge_flow` (225) | 1D shallow-water HLL: `_hll` (26), `channel_step` (59) |
| Hydraulic jump | not resolved (steady GVF energy-based, `_req` 115) | captured by HLL shock |
| Heat transport | `thermal/channel_heat.py` `reach_heat` (66), `segment_heat` (46), `heat_transfer_coeff` (21, Dittus–Boelter) | same coefficient `htc` (52) inside `channel_step` |
| Geometry evolution | `evolution/network.py` (299 lines; `Network.run`, geometry update from ice heat flux), lake storage `geometry/storage.py` | `channel_step` with `melting` flag; `RefNetwork.run` (198) |
| Mass addition (melt) | meltwater not added to flow in RTM; check in network | check `channel_step` (to be read) |
| Conservation | per-run water/mass/energy residuals in campaign records | same |
| Graetz laminar benchmark | — | `graetz_benchmark.py` (`graetz` 18, `reduced` 52) |
Both are **depth-averaged, 1D (along-reach)** models with a bulk-temperature closure. Neither resolves vertical structure, free-surface/air–water interface, thermal boundary layers or phase-change normal velocity, so neither is a vertically resolved CFD model and neither can validate the closure it embeds. They are the baseline for item 3 and cannot serve as the CFD substitute. Items 1's remaining entries (mass-addition, friction and forcing details of `channel_step`, restart support) still need a line-level read.

## What was not done
Items 2, 4–10 are not started except as above; nothing in the milestone is marked complete. No parameter was changed.

## Proposed order once the packet arrives (independent of it: A, B)
A. Prediction CLI around the existing reach model (item 3): `python -m ice_network_v2.predict case.json` with forcing-vs-prediction labelling, no future-geometry inputs, fixed width ⇒ no width output.
B. Toolchain audit: try a pinned OpenFOAM build (container/apt) and demonstrate a fixed-geometry laminar Poiseuille/Graetz check before any melting.
C. After the packet: verification scripts unchanged, comparison table, data audit.
