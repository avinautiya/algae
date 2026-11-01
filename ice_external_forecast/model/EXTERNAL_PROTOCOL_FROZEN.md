# Rochemelon 2005 external thermal-component test: model freeze

This fixes the model **before extracting numerical Rochemelon bed-lowering targets**. The source paper has previously been read; this is not a preregistration or a claim that all research choices are blind. Rochemelon targets will not calibrate coefficients, choose a model, choose geometry or select forcing reconstruction. It is external to Swiss model training only if the input/target audit confirms no Rochemelon incision-derived quantities enter the predictors.

## Frozen Swiss transfer models

Primary candidate: Swiss V5 linear-gap nominal local-geometry scaled DB0.4, coefficient **4.3956164544491925**. Retain the five other families: constant normal melt rate **3.4944254768079555e−6m/s**, constant h **6431.4136285363775W/(m²K)**, unfitted DB0.4, unfitted fluid-cooling DB0.3, unfitted smooth-friction Gnielinski. No coefficient is updated from Rochemelon.

The statistical constant-rate reference remains **uncapped**, preserving its original identity. It reports its requested ice removal and heat demand but does not receive a fabricated closed physical heat ledger. Also report a separately named heat-limited constant-rate physical variant; never substitute that variant silently for the statistical baseline. Every h-based family uses the same geometry and finite advective heat supply.

Frozen properties: water density1000kg/m³, cp4210J/(kgK),nu1.8e−6m²/s,Pr13.5,k0.56W/(mK); ice density917kg/m³, latent heat334000J/kg, ice cp2100J/(kgK). These0°C engineering constants are explicit; no site fitting is allowed. The smooth Gnielinski Darcy factor is `(1.82log10(Re)-1.64)^−2`, not an externally incision-fitted roughness coefficient.

## Finite heat and spatial target matching

Each rectangular reach must have independently sourced length L,width w and depth d. `A=w*d`, wetted perimeter `P=w+2d`, hydraulic diameter `D=4A/P`, velocity `v=Q/A`, `Re=vD/nu`, and `h=kNu/D`. A known lake-inlet temperature cannot be applied unchanged at a downstream target. For an imposed uniform reach, `T_out=T_in exp[−h P L/(rho_w cp Q)]`. Segmented reaches pass outlet water/temperature to the next reach. No fresh heat is added between reaches unless an independently sourced additional input is explicitly modelled in a future version.

Removed heat is exactly `rho_w cp Q(T_in−T_out)dt`. Under the declared uniform wetted-surface h assumption, bed receives fraction `w/P` and sidewalls receive `2d/P`. Record both, rather than calling all wetted-interface melt bed lowering. Mean vertical bed lowering is bed-removal volume divided by bed footprint `wL`; section outlet temperature and every watt/joule/volume ledger are retained. Local interface recession is not full glacier mass balance or a solved moving-cliff/free-surface CFD system.

The heat-limited constant-rate variant has constant wetted-interface heat demand until it exhausts available advected heat. Its local cooling is linear until0°C, rather than falsely using an exponential coefficient. Requested and delivered heat are both reported. The uncapped statistical variant remains outside the closed physical ledger and must be flagged if its energy demand exceeds the inlet supply.

Ice temperature must be explicitly supplied.0°C is a declared engineering assumption if no cold-ice observation exists; optional−2/−5°C cases are engineering sensible-heat sensitivities only. The removed-volume cost is `rho_i[L_f+cp_i(0−T_i)]`, split into latent and fresh-ice sensible demand. This is **not** a cold-ice conduction PDE, a growing cold reserve or a validated refreezing treatment. Negative inlet temperature/refreezing is refused.

## Data/applicability gate before any target score

Require supported forcing intervals with an explicit common clock, nonnegative inlet temperature, positive discharge, all relevant section lengths/widths/depths, and the along-flow distance from temperature measurement to each target. Geometry needs source provenance, not an invented default. Every forcing interval needs its own source/basis; point readings cannot be relabelled complete interval means. No temporal extrapolation, silent interpolation, arbitrary diameter or assumed1m field width. A provided source-assumption scenario is explicitly conditional, never relabelled measured geometry.

Initial/boundary geometry remains fixed during each imposed forcing interval. If time-varying Q,d or geometry is unavailable, refuse the claimed time-resolved prediction; a separately supplied sourced piecewise-constant scenario may run but is labelled as such. The executable refuses gaps/overlaps in a requested continuous packet. No observed evolving geometry derived from target incision may enter a supposedly prognostic prediction without declaring that conditioning and rejecting an independence claim.

Common numerical domains require Re≥10000,Re≤5e6 andPr within both standard gates. These limits do not validate transfer of tube correlations to a rough open ice channel. Hydraulically unresolved controls, free-surface effects, atmospheric exchange, streambed motion and changing cross-sections remain scope limits. No target score runs until the parent records the input/observation match gate and target-definition freeze.

Primary source-clock interpretation must follow the recovered source metadata. If the source clock remains ambiguous, report fixed−2/0/+2h alignment scenarios only where input coverage supports each, never choose the best by target error. Width/depth values must come from source-supported measurements or source-assumption bounds; there are no universally invented width/depth perturbations. Every supplied geometry case is retained. Cold-ice cases are likewise engineering scenarios, not confidence intervals.

For scoring after the gate, retain all models and all independently matched target intervals/sections. Compare identical footprint/time integrals. No selecting a station, fitted source-paper coefficient, imputer, favorable date or geometry after inspecting Rochemelon errors. A result can falsify transfer; consistent environmental superiority must not be manufactured by changing this freeze after target extraction.
