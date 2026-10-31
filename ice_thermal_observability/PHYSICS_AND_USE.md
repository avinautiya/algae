# Physical interpretation and integration gate

## What is actually resolved

The transverse temperature profile is resolved while longitudinal velocity is prescribed. An insulated/free-heated top and temperate cold bed define an actual conduction/advection boundary-value problem, not an empirical Nusselt fit. No Navier–Stokes, turbulence or evolving ice interface is solved here. The shear-free parabolic velocity is a conditional laminar-film base state.

For the demonstration H=0.002 m and U=0.025 m/s, Re_h=50 and Pe_h≈376. This is a thin laminar-flow regime. It is not a stand-in for the Re≈1e4–1e5 glacier channels targeted by the active CFD project. With the stated physical properties, the thermal-development length UH²/alpha is about 0.752 m. The probes and target occupy fractions of this length, so unresolved transverse state can affect downstream heat transport.

## Melting feedback has a physically interpretable scale

Let Q=UHW. With initial temperatures bounded by T_cap and additional total supplied surface power P_s,max, energy alone bounds the ratio of generated meltwater mass flux to inlet mass flux by

    (m_dot_melt)/(rho_w Q)
       <= cp T_cap/Lf + P_s,max/(rho_w Q Lf).

This follows because latent power cannot exceed inlet sensible plus extra supplied power in the nonnegative steady heat model. In the no-source demonstration the bound is approximately 0.00253; the specified 0–100 W/m² source over 0.08 m and width 0.1 m adds at most 0.8 W, raising the bound to approximately 0.00300.

This is a mass-flux screening bound, **not a bound on the error of the regional heat forecast**. Small total meltwater addition does not prove that a local thermal or hydraulic perturbation is harmless. The tool separately flags region-mean lowering relative to depth. Even that is a user screen, not a theorem that interface feedback is negligible.

## Why this can help the physical project

The exact observation-placement theorem distinguishes aggregate thermal calibration from support for a local ice-erosion forecast. The LP supplies explicit alternative profiles and an evaluated dual bound rather than inferring uncertainty from the scatter of a few optimizer starts. It identifies which regional predictions remain robust to unresolved thermal state under the declared class.

The strongest conventional comparator here is a two-mode Graetz model fitted to the same observations. It can have smaller average error on the synthetic smooth profile distribution than the interval midpoint. The midpoint's different guarantee is optimal scalar worst-case error over a broader admitted profile class. This package does not demonstrate a superior average glacier predictor.

The target-weighted current-profile moment specifies the thermal information a particular forecast needs. Implementing it would generally require a resolved cross-stream profile or an independently validated state estimator, not one convenient point thermometer. Changes in geometry/flow require a new kernel.

## How to connect Claude's CFD later

For a fixed verified velocity/interface history with a passive linear thermal equation, derive the map from admissible thermal inputs to both observation quantities and **integrated bed heat**. Include known affine boundary offsets. The same LP/duality machinery applies to that map, after linear superposition, units and conservative exchanges are checked by direct replays. It would then give conditional bounds for the CFD-derived operator rather than the straight-film operator.

When thermal buoyancy, melting, geometry or free-surface feedback changes the operator with the unknown state, this linear envelope is not automatically valid. Derive a justified discrepancy enclosure or restrict to a verified linear neighborhood. The current theorem cannot substitute for that work.
