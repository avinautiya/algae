# What temperature validation can and cannot establish

For constant homogeneous material, Cv*dT/dt=k*d2T/dz2 with two prescribed temperatures depends on alpha=k/Cv. Multiplying k and Cv by the same positive factor leaves every predicted temperature unchanged but multiplies heat flux by that factor. Excellent held-out interior temperatures therefore do not determine absolute basal heat, melt or conductivity without an independent scale observation. This is an elementary heat-equation identifiability fact, not new mathematics.

The helper `thermal_identifiability.py` supplies an independent harmonic Dirichlet solution. Five tests check boundary values, conductivity/capacity scaling and a flux-moment identity. It does not fit field data.

For a column0<=z<=L and constant material, integration of the heat equation against z/L gives

    qbottom = k(Ttop-Tbottom)/L - Cv*dM/dt,
    M=integral[(z/L)T]dz.

The integrated difference between signed dynamic basal heat and a quasi-steady gradient is an endpoint storage-moment term. Long-season ablation may not distinguish their dynamics even when their hourly interior temperatures differ greatly. Temperature mode convergence also does not automatically establish derivative/flux convergence; basal flux needs an independent precision check. Clipping to positive potential melt changes this signed identity and must be reported separately.

The primary temperature campaign uses two measured probes as boundaries. That is a bounded-subcolumn effective-temperature test, not proof of an unmeasured basal ice temperature. The optional ablation transfer instead assumes a melting-temperature boundary and learns an effective conductivity scale from a different linked stake. Do not describe those two observation operators as the same independently validated experiment.

Original metadata identifies the2020b22cm `Base` as constructed0°C, not a thermometer. It cannot serve as a measured boundary or target. Raw2020sensor coverage improves the missing-temperature problem but leaves uncertain early sensor depths and date-only stake times. The old author's processing-code timestamp differs from the current data's start date; it is not independent evidence resolving that date.

P12's profile/stake thickness difference is spatial representativeness and measurement uncertainty, not proof that its actual deepest thermistor was inside ice. The geometry scenarios describe alternative assumed domains, not retroactive diagnoses of which sensor was physically wrong. General field debris-thickness uncertainty is approximately2cm; the2020excavated profile's quoted1cm uncertainty should not be automatically assigned to all2023probe positions.

Any fitted k and implied Cv=k/alpha must be reported as effective parameters unless independently measured material information supports their physical interpretation. The published profile conductivity values cannot be used as independent validation if they were estimated using the same seasonal mass loss. No glacier-wide, wet-film, cliff-lifetime or publication-level conclusion follows from temperature transfer alone.
