# A separate scientific question for the verified ice-channel model

## Candidate contribution, not a claimed result

Can short-lived thermal–hydraulic amplification concentrate incision in an ice-step configuration that is asymptotically stable, and does accounting for that response improve finite-horizon incision predictions?

The useful output would be a tested forecast of **where, how much and over what time incision concentrates** after a physically admissible perturbation. It would explain a limitation of an equilibrium/local heat-transfer representation if such a limitation is actually present. A graph of abstract singular values alone is insufficient.

This study can be prepared independently of the active solver milestone. No physical production runs are launched by this package. The executable response tool currently has manufactured validation examples only.

## Priority and differentiation

Finite-time/nonmodal analysis is established in fluid mechanics and sedimentary morphodynamics. In particular, [Caruso et al. 2016](https://doi.org/10.1103/PhysRevE.93.053110) analyze transient growth and unsteady river bedforms. We cannot claim to invent the mathematical method or its morphodynamic use.

[Yokokawa et al. 2016](https://eprints.lib.hokudai.ac.jp/repo/huscap/all/64599/Yokokawa_et_al-2016-Journal_of_Geophysical_Research__Earth_Surface.pdf) study cyclic ice steps and initial instability. The previous packet also identifies more recent developed-step work and turbulent ice morphodynamics. The specific opening to investigate is the role of a *transported thermal state plus moving hydraulic controls* in finite-time incision of free-surface ice steps, beyond those closest models. Further full-text forward-citation checking is required before any priority claim.

If previous work already establishes this exact mechanism, the project must identify a narrower unresolved prediction or drop the novelty claim. A search with no exact phrase matches does not establish absence.

## Mathematical formulation

Let X collect a physically defined state: water depth and velocity, thermal degrees of freedom, bed geometry, and jump position if a sharp jump is represented. A finite-volume free-surface solver may instead resolve the transition; its sensitivity then needs grid/time convergence.

Linearize about a steady or steadily translating base state, with the appropriate frame and boundaries:

```
d(delta X)/dt = A delta X
delta X(t) = exp(A t) B a
y(t) = C delta X(t)
G(t) = max_{||a||_2=1} ||C exp(A t) B a||_2^2
```

Here B defines admissible *initial-state* perturbations and C extracts the intended observable, such as quadrature-weighted bed lowering in a fixed physical region. This formula is not a general time-dependent forcing gain: forcing requires integrating the propagator against the forcing history. The current code evaluates the autonomous initial-condition problem only. A changing base flow requires a time-dependent propagator, not exp(A t) with a guessed constant A.

Scale states using declared physical quantities before forming A, B and C: e.g. depth/H, temperature/Theta and bed displacement/L_b. Include cell quadrature weights in the corresponding coordinates or maps. Otherwise arbitrary units or grid refinement can manufacture apparent amplification. Do not label this generic quadratic norm a thermodynamic energy unless it has been derived as such.

For full-state Euclidean coordinates, G(t)>1 can occur even with decaying eigenmodes. For different input/output spaces, G is a transfer gain; comparing it with 1 requires the documented physical normalization. The program does not classify eigenvalue stability. Neutral translation and conservation modes must be treated explicitly, rather than mistaking an infinitesimal pattern shift for growth of incision.

The model remains subject to the water–ice power exchange and mass accounting already derived in the V2 packet. A heat-transfer mechanism must pay for melting with supplied enthalpy, retained dissipation, atmospheric power, storage or reduced export. More concentrated incision need not mean more total melt.

## Build the operator only after the active verification gates pass

1. Select one independently verified base state in an explicitly declared regime. Begin with temperate clear-water ice; cold-ice conduction is a later nested extension.
2. Derive a tangent operator or a reduced projection from the actual coupled equations. Do not substitute a hand-chosen 2x2 matrix and call it glacier physics.
3. Verify derivatives using central perturbations over several amplitudes and independent time evolution. Maintain consistent inlet/outlet conditions, mass and thermal inventories. Report derivative sensitivity near wet/dry changes or moving jumps.
4. If using a sharp moving jump, include jump displacement and linearized jump conditions. Smearing a displaced discontinuity on a fixed grid can create norm-dependent apparent growth.
5. If projection reduces the system to <=32 variables for the supplied tool, validate the reduced response against the full model's *observable*, using held-out perturbation shapes and times. Agreement only for training trajectories does not validate the projection.
6. Define B from permissible uncertainty or controlled disturbances. A positive water-temperature perturbation cannot be introduced without accounting for its heat. For a matched-heat comparison, impose the matching constraint on B or the paired initial states.
7. Define C before scoring: e.g. spatial bed-lowering error in metres under fixed quadrature, total melt, or control lowering. Do not select C afterward to maximize an attractive effect.

## Three competing explanations and decisive controls

**Equilibrium geometry effect.** A step redistributes shear and wall heat transfer, but instantaneous hydraulic and bulk-thermal descriptors suffice. Compare against a fairly calibrated local closure on the same geometry and forcing.

**Transported thermal-state effect.** Upstream thermal development and recirculation carry information omitted by a local equilibrium closure. Test physically generated matched local states, then compare wall flux and future incision. The prior packet's constructed temperature profiles only motivate this test; they do not prove environmental persistence.

**Coupled morphological amplification.** Ice lowering modifies the flow and thermal field, which alters subsequent lowering. Compare the fully coupled tangent response with a one-way calculation: retain the effect of heat on bed lowering, but suppress the bed-to-flow/thermal feedback blocks. Do not freeze the output bed entirely, which would make the comparison tautological.

Use equal initial water/heat budgets where the question requires them. Repeat at several perturbation amplitudes; check whether the linear prediction improves as amplitude decreases and whether the physically relevant response persists at finite amplitude.

## Forecast endpoints and falsification

Before a later campaign, fix geometry/forcing holdouts and score:

* Bed-lowering profile at a stated horizon, in metres and quadrature-weighted error.
* Location and magnitude of the incision maximum, with spatial resolution uncertainty.
* Integrated melt and downstream heat export, to distinguish redistribution from additional melting.
* Control lowering and discharge only where the evolving control is causally connected to a lake/network boundary.

Compare equilibrium closure, thermal-state model, one-way morphology and fully coupled model with matched fitting information. Selected cases need independent discretization or materially different numerical settings; a model compared only with its own linearization is insufficient.

Reject the proposed mechanism as a useful predictive addition if its response is below numerical/closure uncertainty, vanishes in the relevant regime, does not survive operator/norm refinement, or fails to improve the predeclared held-out observable. No universal percentage improvement is asserted in advance.

## Environmental interpretation

Temporary growth over a water-transit time may have negligible melt impact if the morphology time is much longer. Estimate and compare water transit, thermal relaxation, jump fluctuation, incision, forcing and ice-deformation timescales using physical properties. Retain physical latent heat. A manufactured matrix's time coordinate is not a glacier timescale.

Use dimensionless overlap and selected turbulent checks before transferring laminar behavior to glacier reaches. Report uncertainty in roughness and thermal closure. Finite-time bed response has potential environmental relevance only when it materially affects location/rate of incision or a demonstrated drainage control. It does not establish a regional melt or flood-risk forecast by itself.

## Deliverable once the mechanism is established

A reduced, conservative predictor supplied with a transient-response regime diagnostic, tested on independent channel geometries and forcing histories. The exact mathematical method need not be unprecedented for a physical finding to be important, but the mechanism and measurable predictive benefit must be demonstrated. Current code prepares that analysis; it has not delivered that glacier finding yet.
