# Independent check of the fixed-interface initial-state bound

The comparison statement matches the actual passive fixed-interface complementarity column in `kennicott_model/column.py` under **identical** material coefficients, capacity, geometry, time discretization and surface-temperature forcing. Let the hotter initial debris and ice state be componentwise no colder, with all initial ice temperatures<=0. The two solutions' temperature order persists. Each backward-Euler conductive matrix has positive diagonal, nonpositive off-diagonal and positive storage, hence nonnegative inverse. Within the warm branch, both temperatures and the net interface supply `Gd*D_last+Gi*I_first` are monotone in the previous state. If the colder solution is warm, the hotter solution is also warm. If the colder solution is cold but hotter is warm, its interface0 is above the cold interface and the discrete maximum principle orders both subcolumns. If both are cold, the joined passive solve preserves order.

Thus instantaneous latent release in the hotter run cannot be smaller. Subtract the two whole-domain energy ledgers:

`ΔL = ΔE_initial + ΔE_in − ΔE_deep − ΔE_final`.

Under the same surface Dirichlet history, `ΔE_in=-Σdt*Gtop*ΔD_top<=0`. Under the same deep0°C boundary, `ΔE_deep=Σdt*Gi*ΔI_bottom>=0`. Positive storage weights and componentwise final order give`ΔE_final>=0`. Consequently

`0 <= ΔL <= ΔE_initial`.

For a2m ice column and uniform initial-temperature difference1K, with source model capacity917*2100J/m³/K and latent heat333700J/kg, the corresponding maximum additional latent release is11.5415kg/m²=11.5415mmwe. The random-pair script tests these signs; it does not establish measured environmental cold content or novelty of comparison principles.

This argument does not permit different fitted debris capacities in the pair, moving melt geometry, ice advection, evolving porosity/moisture, temperature-dependent coefficients that break comparison, above-freezing liquid sensible storage, or externally altered surface forcing. It does not say all initial conditions inside a scalar energy ball are componentwise ordered. A componentwise cold/warm envelope can be used if independently justified; arbitrary equal-energy distributions require another analysis. No field precision bound follows from an engineering choice of initial temperature range.
