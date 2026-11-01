# Fixed-area energy contract

For bare ice, use `BareIceState(cold_reserve_j_m2, melted_kg_m2=0)` and
`step_bare_ice(state, latent_heat_j_kg, dt_s, atmospheric_w_m2, water_w_m2)`.
It applies both fluxes directly to ice enthalpy, without an artificial debris
capacity or temperature. Its `StepResult` uses the same ledger with zero debris
storage. It does not diagnose a surface temperature.

`synthetic_demo.py --output synthetic_ledger.json` compares warm/cold bare ice
and warm/cold debris under three days of explicitly synthetic alternating flux.
It writes per-step states and signed ledgers plus integrated ledgers; optional
`--figure synthetic_comparison.png` requires Matplotlib. Scenario parameters
are illustrative, not observations or calibrated predictions.

`energy.py` has no dependencies beyond Python's standard library. Call
`step_column(ColumnState(...), ColumnParameters(...), dt_s, atmospheric_w_m2,
water_w_m2)`; the returned `StepResult.state` advances the column. All quantities
are per square metre of a constant face area. Water flux must be nonnegative.

Atmospheric heat is signed and enters the debris store. Backward Euler couples
that store to a zero-C reference interface. Interface energy plus water heat pays
the finite cold reserve before latent melt. Negative interface energy increases
cold reserve. Meltwater exits; no freezing, finite ice exhaustion, geometry change,
or cold-ice PDE is represented. In particular, outstanding cold reserve does not
alter the reference temperature used for conduction. This is an explicitly
lumped closure, unsuitable for resolving cold-ice temperature gradients.

The signed ledger is atmospheric + water = debris storage + cold storage +
latent energy; cold storage is the decrease in positive cold reserve, so it is
negative during cooling. `residual_j_m2` exposes numerical closure error. Large
steps remain stable but can reduce temporal accuracy. Supply external fluxes
appropriate to the step; the atmospheric helper does not feed back automatically.

`atmospheric_flux` requires every observation and coefficient explicitly. It uses
face-projected shortwave, gray-body longwave and neutral bulk sensible/vapor
exchange. Specific humidity is kg/kg. See its docstring for equations, constants,
temperature bounds and omissions. No site observations or calibration defaults
are embedded. Tests: `cd forcing && /usr/bin/python3 -m unittest -v`.
