# Versioned integration repairs — candidate branch

This release implements the concrete defects found by the independent audit. It does not reverse the failed field-superiority verdict or certify novelty, journal suitability or competition outcomes. Claude's active branch, frozen results and jobs are unchanged.

## Implemented and checked

* **Phase-source solver:** `phase/swe2d_v5.py` uses the exact pinned v4 hydraulic kernel and a new physical-time split step. The same heat and duration drive the heat sink, bed loss and density-consistent meltwater source. Momentum and sensible heat of newly born stationary-ice water are zero. Explicit resets, ghost exclusion, repeat-consumption guards and separate hydraulic/morphological time ledgers prevent false conservation claims. Nineteen tests pass.
* **Metadata gate:** `contract/` replaces the truthiness check with typed fields, explicit provenance and documented numerical requirements. It has a backwards-key-compatible adapter and explicit no-coercion migration. Twenty-two tests pass. Complete metadata never certifies environmental accuracy or independence.
* **Evidence repair:** corrected proposed release documents and a versioned claim registry distinguish nominal/marginal bias, shared histories, contact/strength/duration endpoints, current/legacy solver evidence and numerical/statistical uncertainty. Twelve regression tests pass. The full-v4 matched-melt gate still correctly refuses missing evidence.
* **Fixed-footprint observation:** `roi/` exports float64 beds and a fixed rectangular mask with exact bounds, source/configuration hashes and byte-verified records. It rejects nonfinite values, changing supports, forged heat/geometry records, incompatible metadata and partial/legacy bundles. Thirteen tests pass.

The combined verifier reports **66 passing tests**. Three short transient physical-time runs also execute the actual hydraulic kernel at two meshes and two step sizes. Their fixed-region geometric melt and heat-implied melt agree within about1e-18m³; combined ice–water mass residual is below3e-14kg. Those tiny two-second numerical checks are not cliff-event convergence, environmental accuracy or a production campaign. See `evidence/INDEPENDENT_INTEGRATION_REVIEW.md` for a separate audit of the final arrays, hashes and exploit refusals.

## Changes installed on the candidate branch

`INSTALL_FILEMAP.json` lists the actual candidate-branch changes. New phase solver files are added alongside, not in place of, `ice_cliff_survival/reactivation/swe2d.py`. The metadata adapter/core/template and corrected evidence documents replace their release counterparts only on this new branch. Original sources are preserved in this packet and Git history. The existing frozen simulations and their failed accuracy verdict are not regenerated or relabelled.

The new phase solver v5 is unrelated to the earlier Swiss point-model V5 name. Existing `cases.py` and demonstration commands retain their historical solver. Opt-in code uses `from ice_cliff_survival.reactivation.swe2d_v5 import Sim` and `sim.step_physical(dt)`.

## Execute

From this packet:

```sh
python3 verify_repairs.py
python3 roi/run_v5_verification.py --out /path/to/NEW-output-directory
python3 contract/validity_contract.py contract/METADATA_TEMPLATE_V2.json
```

The blank template intentionally exits2. Metadata completeness is not environmental approval. The observation exporter writes the manifest last as a commit marker; a filesystem crash may leave an orphan data file, which is refused rather than reused. It is byte integrity, not external authentication of a scientific run. NumPy and Numba are required for the solver tests; contract/evidence checks are standard-library based.

## Physical validity boundary

Native physical-time mode currently supports **bare, hydraulically connected, temperate ice**, with local melt. It refuses covered/porous phase sources and the artificial uniform-incision counterfactual. Where melting occurs below an impermeable deposit, water cannot simply be placed above it without an independently justified drainage or storage model. Quasi-steady legacy-compatible calls record that transfer as an assumption and do not receive a physical-time certificate. The new solver does not resolve cold ice, fracture, collapse, pore pressure, a vertical cliff face, full mechanical/ALE energy or natural boulder transport.

## What is still not solved

The historical v4 matched-footprint campaign is not present as authenticated raw output. The new exporter does not invent or retroactively certify it. Field comparisons still fail the declared superiority rule. The physical-time solver has not been tested for21-day event prediction or full cliff melt, and its computational cost is not the cost of the accelerated quasi-steady runner. No claim of faster environmental prediction follows.

The research candidate remains a precise finite-head, heat-limited recovery boundary under changing runoff. A stronger evolving-topography comparator, matched initial hydraulics and independently frozen endpoints remain necessary before a novelty claim. Software repairs remove avoidable errors; they cannot guarantee a favorable discovery or an award.
