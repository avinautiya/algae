# Prior art and overlap: threshold-aware averaging for outlet activation / drying

## 1. What was actually read

| Source | Read | Not read / limitation |
|---|---|---|
| Kingslake, Ng & Sole (2015), Ann. Glaciol. (attached PDF) | all 15 pages | – |
| Jarosch & Gudmundsson (2012), TC 6, 493 (attached PDF) | all 11 pages | – |
| Ogier et al. (2021), supraglacial drainage (attached PDF) | pp. 1–15 incl. Appendix A | pp. 16–18 (reference list / end matter) not read |
| Schoof et al. (2023), J. Fluid Mech. (attached PDF) | pp. 1–10 (introduction, model, conjecture sections as extracted) | pp. 11–59 and supplementary §5.2 NOT read; statements about it below rest on pp. 1–10 only |
| Threshold-aware averaging / hybrid systems (Filippov; Bogolyubov slow–fast averaging; stochastic-threshold fluvial incision: Tucker & Bras 2000, Lague 2005, DiBiase & Whipple 2011; reservoir theory, Moran/Prabhu) | abstract / textbook-level knowledge only; full texts not obtained (web fetch to publishers was blocked) | any claim about their detail is from memory and flagged "not verified" |

## 2. Overlap table

| Topic | Kingslake et al. 2015 | Jarosch & Gudmundsson 2012 | Ogier et al. 2021 | Schoof et al. 2023 (pp. 1–10) | This work |
|---|---|---|---|---|---|
| Storage | single lake/reservoir, stage ODE with inflow from the channel balance (§2–3) | lake volume–stage law with outlet tunnel (Eq. 5 region) | supraglacial lake/stream phases; storage in lake behind spillway (App. A) | reservoir with steady water supply | lake of area A, stage ODE `A dH/dt = I - Q1 - Q2` |
| Inflow | constant Q_in in the stability analysis (§3); variable inflow only qualitatively | seasonal prescribed forcing; flux set to zero in half-year (Eq. 10) | diurnal / seasonal melt, phases of spillway on/off | steady | arbitrary hourly series (spells + diurnal) |
| Outlets | one channel; multi-channel competition mentioned qualitatively (§6) | one tunnel | spillway channel + dry channel (phase description) | single outlet (instability conjecture) | two outlets with different sills, wet/dry state per outlet |
| Thermal treatment | no sensible heat from the lake | prescribed `dθ/dl` (Eq. 5) | Nusselt-number fits (Nu correlations) for the heat flux | dissipation-only incision | sensible heat of a mixed lake (Dittus–Boelter), lake temperature ODE, sill lowering |
| Outlet threshold | stage vs. overflow level in a linear stability sense | none | spillway activation as a phase transition (descriptive) | none | `I* = Q1(s2)`, wet iff stage > s2; crossing analysis |
| Analytical result | Λ stability parameter for oscillation | analytical period of lake drainage cycles | scaling of phase durations | conjecture for single-outlet instability | comparison-theorem bounds, CUSUM screen, exact crossing time by quadrature |
| Averaging of forcing | not addressed | not addressed (forcing prescribed) | not addressed | not addressed | **the question studied here** |

Exact equation numbers beyond Eq. 5 and Eq. 10 (Jarosch & Gudmundsson) are not listed because they were not carried through verbatim from the readings; section-level references are given instead.

## 3. Closest methodological relatives (not verified at full-text level)

* Hybrid/Filippov systems and Bogolyubov averaging: the correct averaged system for fast forcing is obtained from the occupancy measure of the forcing, not its mean, when the right-hand side is nonlinear or switched. This is the general statement that the mean is insufficient.
* Comparison theorems for scalar ODEs: if forcing is always above/below a steady threshold flow the trajectory stays on the corresponding side; standard.
* Reservoir (Moran/Prabhu) theory: the stage is a reflected process, so wet/dry occupancy is a first-passage/level-crossing problem; CUSUM (Page 1954) bounds are the standard tool for reflected random walks.
* Stochastic-threshold incision in fluvial geomorphology: erosion rate depends on the discharge distribution above a threshold, not on mean discharge. This is the closest conceptual analogue in the Earth sciences.

## 4. What this work could add (claim, to be tested)

None of the four attached papers treats averaged vs time-resolved forcing for outlet activation. The tools above (comparison theorem, level-crossing, CUSUM, occupancy-measure averaging) are standard. The realistic contribution is therefore an **application**: (i) the specific threshold `I* = Q1(s2)` and relaxation time `tau_c = A/g1` for a lake with two sill-separated rated outlets; (ii) a screen (R) for when the mean misplaces the crossing; (iii) a cheap water-storage correction; (iv) documentation of where it fails because heat and sill lowering move the threshold. Whether this is *useful* is decided by the frozen tests; novelty beyond application is not claimed unless the derivation contains something not covered by the items in §3.
