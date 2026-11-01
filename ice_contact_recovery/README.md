# Recovery of stream contact after debris removal

Independent research packet, 1 November 2026. No Claude source, worktree,
checkpoint, frozen campaign or process is modified by this package.

## Research target

Predict the difference between **exporting debris** and **restoring sustained
ordinary-flow contact with a glacier cliff**. The candidate missing state is
the post-abandonment ice topography: a bypass can continue incising while the
toe is abandoned. Complete debris removal need not undo that substrate change.

This is more specific than stream-driven cliff survival, burial, avulsion or
sediment hysteresis, all of which have precedents. Read `prior_art.md` before
writing a headline. A useful new outcome would be a recovery boundary that
improves held-out contact/heat-delivery forecasts beyond a strong model carrying
geometry and water storage. No verified environmental improvement is asserted
by the synthetic examples in this packet.

## Implemented contributions

1. `math/recovery_budget.py`: executable minimum flood-duration and total-water
   forecast, with a sharp proof for a static reservoir/bypass model. A peak
   discharge above the activation threshold is insufficient by itself. Required
   time and total volume grow logarithmically near threshold. Nine tests include
   an independent ODE check and bypass-incision exclusions.
2. `recovery_certificate/bounds.py`: an optimistic finite-width latent-energy
   lower bound on opening a route through a given ice bed. Includes optional
   **discrete whole-cell** shortest-path bound, eight tests, exhaustive small-grid
   comparisons, and an independent adversarial review. The default does not
   falsely charge whole-cell melting for sub-cell channels.
3. `thermal_clock_check.py`: nonlinear two-route matched-dose falsifier. Under
   fixed discharge, fixed debris and homogeneous quasi-steady passive heat,
   temperature history reduces exactly to cumulative thermal dose. Changing
   timing alone must not create a spurious history effect in that regime.
4. `physical.md`: derived forcing-order test with identical joint discharge/
   temperature distributions and total water/heat. Classical flow commutators
   provide an independent small-cycle prediction. This is a candidate physical
   experiment, not a new mathematical technique.

Executed gate: `discharge_order_check.py` finds only 3.29e-10 m of baseline
order-dependent geometry change (0.0000781% of incision). The independent
bracket matches, but this smooth mechanism is environmentally negligible in
that configuration. Do not spend a large CFD campaign on it on this evidence.
The lower-flow sensitivity is a mathematical check, not environmental support.

The methods have deliberately different scopes. The scalar water bound predicts
first zero-depth contact; the spatial energy bound only excludes route opening;
neither guarantees sustained post-flood cliff heat delivery. Do not combine
their outputs into an unconditional survival probability.

## Verification

Run from this folder:

    python3 -m unittest discover -s math -p 'test_*.py' -v
    python3 -m unittest discover -s recovery_certificate -t . -p 'test_*.py' -v
    python3 thermal_clock_check.py --output thermal_clock_results.json
    python3 discharge_order_check.py --output discharge_order_results.json
    MPLCONFIGDIR=/private/tmp/ice_novelty_mpl python3 -m recovery_certificate.demo

For the water-budget figure, run `math/run_examples.py` from its own folder.
Dependencies: Python 3, NumPy, SciPy, Matplotlib. No OpenFOAM, training, new data
download or active solver edit is required. Generated examples are explicitly
synthetic. Classical analytical results are not represented as priority claims.

## Integration order

Read `CLAUDE_HANDOFF.md`. Resolve the ongoing partition/reference gate first.
Then apply one falsifiable question at a time, retaining negative controls and
the current frozen campaign. More compute is useful only after the prediction,
comparison and acceptance rule are fixed.
