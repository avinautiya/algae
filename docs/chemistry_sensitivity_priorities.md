# Chemistry sensitivities prioritised by downstream consequence (2026-10-09)

## Principle

The calibration (`docs/quantity_unit_audit.md`, rows 4–6) fits a band shift dE, a width w, a magnitude f and an Fe fraction φ. The glacier model also discards the computed UV structure below 350 nm (row 7).

A chemistry sensitivity therefore matters for glacier radiation only if it changes something the calibration cannot absorb:

- (i) the **relative** positions and intensities of bands within 350–750 nm;
- (ii) whether the pigment form itself absorbs in the **visible**, where the solar weight and the measured extract–chromophore difference live.

Today the visible absorption is supplied by the empirical Fe increment φ·D.

## Table

| Sensitivity | Job | Status | Absorbed by calibration? | Expected downstream consequence | Priority |
|---|---|---|---|---|---|
| Functional (B3LYP vs CAM-B3LYP) | `L2_CAM_TDA15` (running), legacy CAM TDA15 (complete) | raw: CAM blue-shifted about 0.5 eV, no visible absorption (`records/molecular_contribution/raw_vs_measured.csv`) | mostly (dE) | measured directly in Experiment A with the legacy CAM result. The rerun duplicates the legacy settings (validated checkpointing); it is not killed (about 2.5 h invested). | **measured now** |
| TDA vs full TD (B3LYP) | `L2_B3LYP_TDA15` (queued, interrupted) | — | yes (shift/intensity ratio small for these π→π* states) | low | low |
| Root count 15 → 25 | `L2_B3LYP_TDA25` | — | n/a: extra roots lie < 300 nm, and the glacier input is held below 350 nm | negligible for radiation | low |
| Geometry convergence (provisional B3LYP geometry) | `L2_OPT` → `L2_B3LYP_TDA15_RELAXED` | optimisation not converged (`grad_rms` met only) | shift: yes; relative band changes: partly | low for radiation; needed for any claim about **raw** molecular accuracy (Q1) | medium (mechanistic) |
| Protonation: carboxylate anion (COOH pKa about 3–4; cytoplasmic/vacuolar pH about 5–7, so plausibly deprotonated) | `L2_COO_OPT_TDA15` | queued | partly | moderate: can change relative band intensities and red-shift S1 | **high** |
| Phenolate (OH deprotonation; pKa likely > 7) | not defined | — | no, if it creates visible absorption | potentially high, but chemical plausibility at in-vivo pH is unestablished. Screening only after the carboxylate result. | screening candidate |
| Fe(III)–purpurogallin catecholate complex | `FE_CAT` | queued (interrupted earlier by OOM) | **no**: a computed complex could replace the empirical Fe increment D | **high**: the only route by which chemistry could predict the visible absorber. Requires state-character (LMCT/LLCT/d–d) analysis before labelling. The measured PG–Fe spectrum (Procházková 2025; band ~593 nm) is an independent Q1 target for it, because the computed complex is not fitted to it. | **highest** |
| Fe tropolonate, xTB geometry | `FE_TROP_XTB` | queued | no | screening only (xTB geometry) | medium |
| Conformers / tautomers | not defined | — | partly | unknown; no population data at in-vivo conditions. Scenarios only, never an invented Boltzmann mixture. | after the Fe and carboxylate results |
| Solvent (PCM water vs explicit / cell interior) | not defined | — | shift: yes | low after calibration | low |

## Decision

1. Reorder the runner queue so that **FE_CAT** and **L2_COO_OPT_TDA15** start before the TDA/root/geometry checks. Job definitions are unchanged; only the order changes.
2. TDA15, TDA25 and the geometry relaxation remain queued for mechanistic completeness. They are not on the critical path of the radiation question.
3. **Screening vs production.** FE_CAT and FE_TROP_XTB are exploratory (idealised complexes, one geometry each). Their spectra are compared with the measured PG–Fe increment as a **shape/position** test only. The fitted φ is never interpreted as a measured complexed fraction.
