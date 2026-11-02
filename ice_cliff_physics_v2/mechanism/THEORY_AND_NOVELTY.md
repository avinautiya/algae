# Wet-debris thermal storage: an adversarial mechanism audit

## Decision

A conservative wet-debris boundary condition is a necessary environmental extension, not inherently a major new environmental discovery. Thermal inertia, debris moisture and coupled heat/water transport already have substantial precedents. The narrow candidate question is whether **a transient thermal buffer changes an irreversible hydraulic-contact/recovery outcome after thermal transients have relaxed**, through evolving ice geometry. A result confined to a warmer initial debris reservoir or a temporary stored-energy difference is not that mechanism.

Existing project packets already implement geometric recovery memory, water/energy exclusion bounds and a quasi-steady thermal-clock falsifier. Do not present those as newly discovered by this upgrade. The older clearance–reconnection question remains more defensible than a broad first thermal-memory claim.

## Passive network and exact identity

Consider one temperate-ice contact with bare conductance Gb and water/debris conductance Gd, water heat-capacity rate F=rho_w cp Q, debris heat capacity C, ice/debris conductance K, prescribed inlet superheat Tin and dry-surface atmospheric power qa. Every conductance has units W/K; C has units J/K. Geometry and coefficients are fixed for this derivation.

Set G=Gb+Gd, alpha=Gd/G, phi=(1-exp(-G/F))/(G/F). The mean stream temperature over the common contact is

    Tbar = phi Tin + (1-phi) alpha Td.

Thus qd=a Tin-b Td and qb=Gb phi Tin+d Td, where

    a=Gd phi,
    b=Gd[1-(1-phi)alpha],
    d=Gb(1-phi)alpha,
    b-d=a >= 0.

The debris equation is

    C dTd/dt = qa + a Tin - (K+b) Td.

The instantaneous-equilibrium temperature is Tss=(qa+a Tin)/(K+b), and total ice power is qb+K Td. Subtracting the full power from the otherwise identical zero-storage equilibrium model and integrating gives, for arbitrary time-varying Tin and qa,

    Qice,eq - Qice,dynamic = Gamma C [Td(end)-Td(start)],
    Gamma=(K+d)/(K+b),   0 <= Gamma <= 1.

This identity requires constant network coefficients. It follows directly from a linear energy balance and is not a historical-originality claim. It is useful as an independent code falsifier. Any enduring cumulative order effect in a fixed network after equal terminal storage violates the assumptions or the implementation. Changing ice geometry/contact breaks the fixed-coefficient reduction and is the mechanism that must be demonstrated separately.

The maximum fixed-network phase-energy discrepancy associated with a supplied terminal temperature excursion DeltaT is Gamma C DeltaT. Its equivalent lowering over a declared ice footprint A is

    delta_m <= Gamma C DeltaT / (rho_i L A).

This is conditional on the supplied temperature bound, not a confidence interval or a guaranteed recovery forecast. A conversion from this displacement to route-capture time additionally requires a demonstrated geometric response and must not be inferred by dividing by an arbitrary mean melt rate.

## Environmental scale screen

The local relaxation time is tau=C/(K+b); NTU=G/F controls finite-water coupling. At the measured Kennicott stream discharges reported in the primary process study, F is large compared with the conductance of a small wetted patch. For NTU << 1, Gamma approaches K/(K+Gd): much of a hot debris reservoir can be exported as warmer water rather than melted into the local ice. Large initial debris temperature cannot be treated as a free source of novel melt enhancement.

For fixed network forcing, compare only equal initial extensive enthalpy and equal imposed inlet/atmospheric energy. Net absorbed water heat may legitimately differ because outgoing water temperature differs; retain that term instead of forcing equal absorbed heat by changing the experiment unnoticed. A replay with matched deposited heat is a diagnostic counterfactual, not an independently deployable predictor.

## Geometric and evidential cautions

- The legacy face-cover fraction refers to the upper face. Natural toe debris is lower-face cover; merely adding a wet kernel to upper cover will usually leave it unexcited at ordinary stream depths. Use explicit upper/lower mask locations and report the actual wetted overlap.
- Uniform face retreat does not resolve undercutting, collapse or full cliff lifetime. Its face-integrated melt is a reduced output.
- A thin upstream hydraulic obstruction and an attached debris coating are separate model elements. Clearing the former must not be described as removing all face debris.
- Freezing/cold interface physics lies outside this temperate kernel. Negative debris temperature or unsupported refreezing must be refused or handled by a separate verified state model.
- A solver refusal, missing forcing or inactive wet patch is not an observed loss/recovery outcome.

## Closest literature, inspected at primary sources

1. Reid and Brock (2010), [debris energy balance with internal heat conduction](https://www.cambridge.org/core/journals/journal-of-glaciology/article/an-energybalance-model-for-debriscovered-glaciers-including-heat-conduction-through-the-debris-layer/21379FE9ED505F092BCA58276652B7A5). Dynamic debris temperatures and conductive delivery to ice are established. The present identity does not make debris heat storage itself new.
2. Giese et al. (2020), [moisture in debris-covered glacier energy balance](https://tc.copernicus.org/articles/14/1555/2020/). Moisture changes debris thermal behavior and energy partition. A generic first wet-debris or heat-capacity claim is inappropriate.
3. Winter-Billington et al. (2022), [simultaneous heat and water transport](https://www.frontiersin.org/journals/earth-science/articles/10.3389/feart.2022.796877/full). Coupled water/heat and moisture-dependent debris thermal properties predate this project. A lumped wet-contact kernel is a reduced closure, not unprecedented multiphysics.
4. Petersen et al. (2024), [stream hydrology and cliff evolution](https://esurf.copernicus.org/articles/12/727/2024/). Stream undercutting, debris evacuation, avulsion and reburial are documented. The candidate quantitative post-clearance recovery boundary is narrower than discovery of those processes.
5. Evatt et al. (2024), [moving-frame ice-cliff model](https://doi.org/10.1017/aog.2024.34). Cliff geometry/frame effects and endogenous slope selection are established; a planar fixed-frame face must not be sold as full cliff mechanics.

This bounded search does not establish historical priority for the proposed coupled recovery experiment. A consequential robust result remains to be computed and checked against an independent spatial model. The wet kernel and linear identity alone do not establish a Nature-level contribution or competition outcome.
