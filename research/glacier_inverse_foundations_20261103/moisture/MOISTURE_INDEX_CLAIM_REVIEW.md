# Moisture-index physical and mathematical claim review

The proposed index is a phenomenological moisture/temperature predictor. Its
accuracy must come from the separately frozen target test. It is not an
energy-conserving melt model, and the verified ice-flux module's energy ledger
cannot be attached to it.

## Positive moisture coefficient is a hypothesis

Water replacing pore air can increase effective conductivity under otherwise
comparable packing. Net melt also depends on heat capacity, evaporation,
freezing/melting, subsurface transport and surface forcing. Those competing
effects do not imply that wetting universally increases glacier melt.

The primary [Giese et al.2020 ISBA-DEB study](https://tc.copernicus.org/articles/14/1555/2020/tc-14-1555-2020.html)
already represents moisture-dependent thermal properties and phase changes. Its
dry/saturated comparison produces more melt beneath dry debris: increased water
capacity lowers the effective diffusivity in that configuration. This is a
direct counterexample to treating positive moisture–melt sensitivity as a
thermodynamic law. It also means adding moisture alone is not historically new.

[Rounce et al.2015](https://tc.copernicus.org/articles/9/2295/2015/)
compares latent-flux treatments and finds that omitting latent exchange can
overestimate ablation. [Petersen et al.2022](https://doi.org/10.1029/2021JF006520)
discusses moisture variation together with conductive, convective and latent
transport. These established mechanisms support testing a moisture covariate;
they do not validate a particular fitted exponent or establish causality.
Moisture can itself respond to meltwater, so contemporaneous observed moisture
supports a conditional hindcast/retrieval, not demonstrated prospective forecast
skill or a causal claim that wetting drove the measured melt.

A restriction`beta>=0` can therefore be declared as an empirical
conductivity-dominated hypothesis motivated and selected on Tower1 development
data. It must not be described as required by thermodynamics or expected to
hold under every wet/freezing/evaporative regime. A negative or failed transfer
result would be informative and must not trigger target-based retuning.

## Exact scope of offset invariance

For `D_beta(I)=integral_I Tplus(t)*exp(beta*(theta(t)-theta_ref))dt`, calibrated
prediction is`Y_train*D_beta(test)/D_beta(train)`. A constant additive change
`theta->theta+a` multiplies numerator and denominator by the same factor, which
cancels exactly. The reference theta_ref cancels for the same reason.

This is **additive-moisture-offset invariance**, not invariance to every sensor
calibration error. A multiplicative moisture gain changes the effectivebeta;
nonlinear response, drift, saturation/censoring, mismatched probe depths and
temperature sensitivity generally do not cancel. State that restricted property
plainly. Source readings outside physical moisture bounds still require flags;
invariance does not establish that they are physical water-volume measurements.

The usual characterization follows if a positive continuous weightw satisfies
`w(theta+a)/w(theta)=c(a)` throughout a common domain: continuity reduces the
multiplicative Cauchy relation to`w(theta)=K*exp(beta*theta)`. This is established
functional-equation mathematics under stated assumptions, not a new historical
theorem or an independently derived glacier law.

After TRAIN calibration,

```
d log(M_test)/d beta
  = weighted_mean_theta(test) − weighted_mean_theta(train),
```

where weights are`Tplus*exp(beta*theta)`. Thus even positivebeta need not increase
every calibrated test prediction: relative wetting/drying of the two periods
controls the direction. The index can summarize correlated seasonal state
without isolating causal conductivity.

## Fresh-test scope

Freezebeta, moisture proxy, positivity restriction, missingness policy,
calibration information and every comparator before Tower2 camera release.
Training a parameter on four already released Tower1 intervals makes those
intervals development data. Reusing their calendar periods at a previously
unopened Tower2 stake tests within-campaign spatial transfer with shared weather;
it is not independent-glacier or unseen-weather validation. If a local Tower2
TRAIN-melt gain is allowed, apply the identical information budget to baselines
and label the later comparison as conditional on that calibration.

No Tower2 camera outcome was opened for this review. The index's professional
utility, predictive improvement and historical novelty remain empirical/prior-art
questions; none follows merely from adding an exponential or proving invariance.
