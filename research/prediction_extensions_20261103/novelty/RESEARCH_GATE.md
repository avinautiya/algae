# Prediction novelty gate and a conservative thermal-path counterexample

## Decision

Do not claim a new general melt law, universal recovery, field superiority or Nature/STSsuitability. The current ideal runoff-pulse campaign is analytically blocked under its common-temperature/fixed-transfer assumptions. More neural-network layers or nuisance parameters do not resolve this.

One defensible candidate for a more informative prediction is **a route-specific thermal-supply correction to a finite-head contact predictor**. This is not a newly discovered mechanism: classical heat advection and stream cooling already exist. A publishable distinction would require showing that a precisely derived, transferable correction predicts a consequential quantity that a strong existing model cannot, with independently supplied physical inputs.

## Established neighbors that must be represented honestly

* Hill & Dow2021, SaDS, https://doi.org/10.1029/2021JF006309 : coupled sheet/channel routing, mass sources from melt and potential-energy dissipation. It is not an omission-free benchmark, but it prevents claiming the first dynamically coupled drainage network. Its channel model simplifies where perimeter melt is distributed.
* Karlstrom, Gajjar & Manga2013, https://doi.org/10.1002/jgrf.20135 : curvature-dependent heat transfer and thermal channel migration already modeled. Merely adding a curvature multiplier is not original.
* Isenko, Naruse & Mavlyudov2005, https://www.sciencedirect.com/science/article/abs/pii/S0165232X04001594 : downstream water-temperature evolution and heat-driven channel melt are prior art. Only the publisher abstract was inspected in this audit; no equation-level priority claim is made from it.
* Petersen etal2024, https://esurf.copernicus.org/articles/12/727/2024/ : stream migration, abandonment, undercutting and cliff reburial are already documented. The qualitative stream-access story is not new.
* Jarosch & Gudmundsson2012, https://tc.copernicus.org/articles/6/493/2012/ : ice mechanics, hydraulic flow and heat transfer have already been coupled for channel evolution. Adding broad multiphysics alone is not novelty.
* Strickland, Covington & Gulley2025 preprint, https://doi.org/10.22541/essoar.175647655.50202412/v1 : the previous pinned novelty audit inspected a closer evolving-topography/debris-transport comparator and identified a finite-head/finite-water-heat distinction. This task inherits that qualified reading, not a claim that all2026follow-ups were exhaustively excluded.

## Exact limited counterexample

For prescribed flowq, heat-capacity rateF=rho cp q, inlet superheatT0, upstream exchange conductanceU and target conductanceH:

    T_target = T0 exp(-U/F)
    P_target = F T0 exp(-U/F)(1-exp(-H/F))
    P_upstream = F T0(1-exp(-U/F))
    P_outlet = F T0 exp(-(U+H)/F)

All energy is allocated once. Target incision is `P_target/(rho_i L A_target)`. For two target routes with unattenuated incision ratesEt0,Eb0, bypass upstream cooling reverses their ranking at

    U_b > F_b log(Eb0/Et0) + (F_b/F_t) U_t,

provided Eb0>Et0 and supplied flows are held fixed. This is elementary algebra from established heat-transfer equations, NOT a new theorem with an originality claim.

The executable example uses the audited ideal pulse flows0.006834 and0.133166m³/s. Its unattenuated bypass incises slightly faster. A bypass upstream conductance exceeding873.99W/K reverses the sign. At twice that conductance the negative gap drift is only about0.0438mm/day. Thus the example shows a missing state variable can change a sign, but **does not establish an environmentally consequential rescue**. It must not be sold as impressive melt improvement.

Prescribed equal target flow is an observation-conditioned comparison, not proof a complete upstream hydraulic network retains the same head under arbitrary path changes. Long paths also change friction and geometry; these must be solved or independently measured. Upstream thermal attenuation represents ice melting there, not energy destruction. Redistribution toward/away from a cliff cannot be called increased total glacier melt without integrating the entire physical support and outlet heat.

## Strong, discriminating experiment, only after data gates

1. Predict local target melt/heat on fixed supports, not cliff lifetime or glacier-scale melt. Supply measured or independently resolved routeq, waterTin, upstream thermal geometry and target geometry. Do not estimateU by fitting the same target melt.
2. Compare (a) finite-head common-temperature reduction, (b) finite-head plus serial thermal paths, (c) an independently discretized advection/heat reference using identical sources, and (d) a strong trained statistical baseline. The reference shares declared closures, so numerical agreement is consistency, not field validation.
3. Hold target initial head/flow partition and thermal support matched, while explicitly accounting for different upstream pressure loss. Report the comparator's initial mismatch rather than compensating with undisclosed calibration.
4. Freeze route-specific temperature and local melt predictions before inspecting test targets. Distinguish upstream loss, target heat, outgoing heat, total melt and local gap drift. Test substantial and negligible NTU, unequal paths and the equal-path limit.
5. Evaluate improvement against the strongest baseline, not only a deliberately temperature-blind variant. Actual field success needs paired upstream/target temperature, targetq, path/contact geometry and same-region ice change. Without these, this remains a synthetic model-structure challenge.

## Kill criteria

* Required thermal-path data absent orUonly identifiable through target melt: no independent field test.
* Attenuation is too weak to exceed numerical/measurement uncertainty, as in the small sign-flip example: no environmental-impact headline.
* A standard stream-cooling model predicts the same result with the same information: no distinct methodological contribution.
* Improved local score disappears under matched physical supports or remains worse than empirical baselines: no superiority claim.
* A negative result is retained. No retuningh,T, pulse duration or test cases to make a rescue happen.

## Delivered execution

`thermal_path.py`, `test_thermal_path.py`, `COUNTEREXAMPLE.json`. Five tests pass: energy closure, serial subdivision, zero flow, exact crossing and local ranking reversal. These are analytic/engineering verification fixtures only. No Claude code, historical result or production campaign was altered.
