# Gate for a calibrated moisture-modulated temperature index

3 November 2026. Development hypothesis only. Tower 1's four released evaluation periods have already been inspected and are now DEVELOPMENT. No Tower 2 camera outcome was opened in this audit.

## Decision

A compact prospective cross-stake test is defensible if the method is described as an empirical moisture-modulated melt index, not newly discovered moist-debris physics, an energy-conserving melt solver, or guaranteed superiority. The useful question is whether a fixed moisture response learned at a drying stake transfers to a wetting stake after one local training-amplitude calibration. That directional transfer is more informative than another within-season fit. It does not by itself establish a causal conductivity mechanism.

The proposed form is

    M_s(I) = Y_s(TRAIN) * D_s(I;beta)/D_s(TRAIN;beta)
    D_s(I;beta) = integral_I max(Tair,0) exp(beta*theta_s) dt.

The reference theta_ref_s cancels algebraically and need not enter computation. beta>=0 encodes a hypothesized local positive moisture modulation. Moisture can also reduce melt through latent exchange and thermal inertia; the sign is a hypothesis, not a universal physical constraint. Output calibration requires nonzero finite training dose and positive finite training melt. This uses measured air temperature and moisture, not remotely sensed flux or full glacier dynamics.

## Primary neighboring work

- Nakawo and Young (1981), [10.3189/172756481794352432](https://doi.org/10.3189/172756481794352432), already models distinct dry and wet debris layers with moisture-dependent conductivity and series resistance. Moisture-dependent insulation is not new.
- Collier et al. (2014), [10.5194/tc-8-1429-2014](https://tc.copernicus.org/articles/8/1429/2014/), couples reservoir moisture, phase changes and energy exchange in glacier debris. A changing water stock is an established physical state.
- Giese et al. (2020), [10.5194/tc-14-1555-2020](https://tc.copernicus.org/articles/14/1555/2020/), explicitly makes conductivity and heat capacity vary with depth, time, moisture and phase in ISBA-DEB. Their results warn that conductivity-only moisture corrections omit important latent and transport terms. Thus the proposed index is a reduced statistical approximation, not a first moisture-aware glacier model.
- Juen et al. (2013), [10.1111/geoa.12011](https://doi.org/10.1111/geoa.12011), experiments on lithology/grain size/porosity and moisture-controlled thermal properties establish the relevant material dependence and its vertical variation.
- Petersen et al. (2022), [10.1029/2021JF006520](https://doi.org/10.1029/2021JF006520), discusses seasonal/hourly moisture changes, nonconductive flux and their implications for apparent diffusivity. A change in fitted effective k is not unique proof of drying as its cause.
- Melo-Velasco et al. (2025), [10.1029/2025JF008360](https://doi.org/10.1029/2025JF008360), already uses local moisture in CRh/CRi property estimation at these sites. Its source reports moisture-sensor accuracy +/-3%, temperature accuracy +/-0.4 C and sensor-depth uncertainty. Avoid claiming that these observations previously ignored moisture.
- Pellicciotti et al. (2005), [primary article](https://www.cambridge.org/core/journals/journal-of-glaciology/article/an-enhanced-temperatureindex-glacier-melt-model-including-the-shortwave-radiation-balance-development-and-testing-for-haut-glacier-darolla-switzerland/E96A8B8D2903523DE6DBF88E2E06E6D9), develops an enhanced temperature-index model with radiation/albedo. Additional process-aware predictors in a calibrated index are an established framework.

A targeted search for exponential moisture-modulated glacier temperature indices and sensor-bias invariant glacier moisture indices did not locate an exact identical formulation. This is not proof of historical originality. Exponential link functions, calibration ratios and invariance arguments are standard mathematics. The defendable candidate contribution is this specific low-input estimator plus a rigorous prospective transfer comparison, if successful.

## Mathematical claims and their precise limits

### Constant additive sensor offset

For theta_tilde(t)=theta(t)+b_s, exp(beta*b_s) multiplies both doses and cancels exactly. This remains true while fitting beta on Tower 1 if every candidate beta is evaluated with the same calibrated ratio. It does not cancel a multiplicative sensor calibration error: theta_tilde=a_s theta+b_s changes the effective beta to beta*a_s. Transfer between differently calibrated moisture probes remains an important assumption.

### Uniqueness is restricted

Suppose one fixed positive continuous uncentered w obeys

    w(x+b)/w(y+b)=w(x)/w(y)

for all x,y,b in a shift-invariant domain. Set y=0. Then w(x+b)=w(x)w(b)/w(0). Positivity permits logs; continuity turns the additive Cauchy equation into log(w(x)/w(0))=beta*x. Therefore w(x)=A exp(beta*x). This is an elementary characterization of continuous multiplicative characters, not new mathematics.

Crucially, if a reference is itself recomputed from the biased sensor, ANY centered function w(theta-theta_ref) is constant-offset invariant. The exponential is not unique among centered models. Its useful distinction is that its calibrated ratio also does not depend on the chosen reference, along with the simple bounded-drift calculation below.

### Bounded nonconstant drift

For fixed beta and |epsilon(t)|<=delta in both training and evaluation, each positive dose changes by a factor in [exp(-|beta|delta),exp(|beta|delta)]. The calibrated prediction therefore changes by at most [exp(-2|beta|delta),exp(2|beta|delta)]. This is a conditional worst-case multiplicative sensitivity bound, not a confidence interval or measured field accuracy. It does not cover re-estimation of beta under errors, air-temperature errors, camera error, missing forcing, geometry or model discrepancy.

Do not clip shifted moisture values, change sample eligibility by moisture threshold, or fit a new local beta and still assert the same exact invariance. A fixed finite-data mask and unmodified covariate are part of the stated property. Use log-sum-exp for numerical stability; clipping exponents breaks exact guarantees at the clip.

## Why the observed drying is not yet a confirmed mechanism

Tower 1's released SM3 means 0.05019 to 0.01737 differ by 0.03282 volumetric fraction. If the source +/-3% means an absolute +/-0.03 fraction bound per period, arbitrary time-dependent errors permit the true difference to vary by +/-0.06, so the sign is not certified. A single constant sensor offset would cancel the difference; the source does not establish that all error is a constant offset. Repeated five-minute readings do not reduce stated instrument accuracy by sqrt(n). Clarify the manufacturer's error semantics before calling the trend physically quantified.

The apparent k values derived from observed melt and local gradients are useful development diagnostics. They depend on the target, unmeasured basal gap, temperature offsets, ice storage and any nonconductive transport. They must not become independent predictors by recomputing k from Tower 2 evaluation melt.

## Fair prospective test

Freeze one primary model and beta using Tower 1 DEVELOPMENT only. Tower 2 gets one native TRAIN camera amplitude and no other target-derived parameters; retain all eligible Tower 2 evaluation intervals as one test family. No evaluation-driven beta, moisture depth, clock or model selection is allowed.

Required comparators with matched development information and local training gain:

1. Plain air PDD (beta=0) and constant rate.
2. Moisture-only exponential dose (replace max(Tair,0) by 1), to test whether moisture merely carries the seasonal trend.
3. Centered linear-moisture air-PDD dose, constrained positive on a predeclared source/measurement envelope; this tests whether exponential curvature is necessary.
4. Power-moisture air-PDD with a single fitted exponent and a fixed predeclared treatment of nonpositive/missing readings. Label its additive-offset sensitivity; do not rescue it after observing evaluation.
5. Gradient index with identical one-TRAIN amplitude and the existing buried-temperature index.
6. A simple calendar-trend-modulated air PDD with one fitted trend parameter, learned on the same Tower 1 development intervals. This is essential because drying and time trend are confounded at Tower 1. The direction-reversal transfer supplies a more decisive test than within-stake MAE.

All fitted modulation parameters use the same development losses, parameter count, bounds and selection budget. Report the complete model family, ties, source-coverage exclusions, camera endpoint operator and unambiguous primary score. Do not call a comparison against these indices superiority to unreproduced professional coupled energy-balance models. A source-equivalent moisture-aware CRh/CRi or ISBA-DEB comparator remains the stronger professional baseline if faithfully executable with the available inputs.

## Physical storage extension is a separate model

For a basal debris slab, conservation is

    Q_ice = Q_into_slab - change(H_slab) - other outward energy.

If C changes with moisture, writing C(theta)*dT/dt alone is not equivalent to dH/dt. Mass entering/leaving carries enthalpy, and pore freezing/evaporation carries latent energy. Source sensors may not observe the whole slab or the exact interface. Adding an empirical moisture factor does not account for these terms. A conservative cold-storage estimator must use a declared enthalpy definition, phase balance and advective/latent boundary assumptions; otherwise preserve the statistical index label.

## Release gate

Proceed to a small fixed predictive test only if the development beta is numerically stable under held-period deletion, finite source uncertainty sensitivities are reported, and missingness/coverage rules are fixed. If beta is set by one observation or its drift bound is so broad that rankings cannot be distinguished, do not promise a useful operational forecast. Freeze before any Tower 2 outcome access. Success is a measured result, not a prerequisite to manufacture.
