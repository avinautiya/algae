# Research framework after the failed field-superiority tests

**Most specific lead:** use a genuinely observed transition from a melting interface to a cooling ice interface as a natural thermal reference, then test whether temperatures constrain absolute debris properties and basal heat without a melt-fitted gain. This is a plausible new glaciological application, not established historical priority or new inverse-heat mathematics. It is a **GO for a bounded synthetic/temperature-only observability study**, with no new camera outcomes required; it is not yet a production melt-model justification.

The prior model's concrete failure was fitting temperature-response diffusivity while requiring a large empirical melt-amplitude correction. Simply changing temperature weights or adding a fitted contact resistance does not remove the underlying scale ambiguity. The new question must identify what independently supplies the missing absolute heat scale.

## Primary prior-art audit

|Original source|Closest established result|Consequence for this project|
|---|---|---|
|Xu & Wang2014, [*Across-plane thermal characterization…*, DOI10.1063/1.4898330](https://doi.org/10.1063/1.4898330), [original-paper PDF](https://pdfs.semanticscholar.org/1e13/6fabe36893a246166fed06e37ea83f0d1a70.pdf), §§II,IV|Multilayer photothermal responses depend on layer/substrate effusivity ratios and interface resistance. Frequency/amplitude/phase sensitivities and contact-property confounding are explicitly studied.|A known substrate as an absolute thermal reference, multi-frequency separation and contact-resistance inference are established outside glaciers. The glacier proposal cannot claim to invent reference calorimetry.|
|Spantini, Cui, Willcox, Tenorio & Marzouk2017, [DOI10.1137/16M1082123](https://doi.org/10.1137/16M1082123), [author manuscript](https://kiwi.oden.utexas.edu/papers/goal-oriented-bayesian-inverse-problem.pdf)|Goal-oriented Bayesian inversion targets a quantity of interest directly; the paper includes a heat-transfer inverse problem and optimal low-rank posterior approximations.|Adjoint/goal-directed melt inference is not itself a new mathematical framework. Gaussian posterior guarantees cannot be transplanted to unknown camera-error semantics.|
|Djatouti, Waeytens, Chamoin & Chatellier2017, [DOI10.26868/25222708.2017.513](https://publications.ibpsa.org/proceedings/bs/2017/papers/BS2017_513.pdf), §§2–3|A thermal inverse model updates parameters relevant to the sought heat flux rather than reconstructing the whole thermal state. It contrasts global temperature calibration and a goal-oriented functional.|Even the narrower idea “calibrate for heat flux rather than unweighted temperatures” has direct prior art. A glacier-specific validated formulation would need a clear additional contribution.|
|Laha et al.2022, [DOI10.1017/jog.2022.35](https://doi.org/10.1017/jog.2022.35), §4.1|Homogeneous/inhomogeneous regression and Bayesian forward inversion of debris-temperature profiles are compared for ablation estimation.|Forward inversion, Bayesian estimation and layered debris are strong existing glacier baselines, not novel merely because combined with a PINN.|
|Beck & Nicholson2025, [DOI10.5194/tc-19-2715-2025](https://tc.copernicus.org/articles/19/2715/2025/), §§3–5,AppendixA; [public code](https://github.com/calvinbeck/TC-DTD)|Numerical truncation, sensor precision, sampling depth and sensor-position errors can bias diffusion-derived parameters; the authors provide a field-design tool.|A generic observability/sensor-placement study is already close prior art. Distinguish **absolute flux-scale identification through phase/reference changes** from improving derivative diffusivity estimates.|
|Rowan et al.2021, [DOI10.1017/jog.2020.100](https://doi.org/10.1017/jog.2020.100), [author PDF](https://pure.aber.ac.uk/ws/portalfiles/portal/101317580/Article.pdf), §§3.4,4.4–4.5,5.1|Winter basal temperatures are negative, while summer interfaces approach melting. Conductivity is derived from summer diffusivity with assumed rock density, specific heat and porosity. Sensor migration and NG2 boundary anomalies are documented.|The natural phase regimes are observed, but this paper does not perform the proposed ice-reference absolute-property inversion. Observation of a cold interface alone is not a new discovery.|
|Giese et al.2020, [DOI10.5194/tc-14-1555-2020](https://tc.copernicus.org/articles/14/1555/2020/); Collier et al.2014, [DOI10.5194/tc-8-1429-2014](https://tc.copernicus.org/articles/8/1429/2014/)|Debris moisture, pore-ice/water phase changes and seasonally varying thermal properties are explicitly modelled.|Winter identification cannot be carried unchanged into wet summer debris. Adding a pore-phase state is established physics; its identifiability and predictive value must be tested.|
|Tarzia1983, [DOI10.1016/S0017-9310(83)80169-0](https://doi.org/10.1016/S0017-9310(83)80169-0), [author PDF](https://www.austral.edu.ar/investigadores/wp-content/uploads/2020/03/Tarzia-IntJHeatMassTransfer-2619831151-1158.pdf)|Inverse Stefan problems can determine unknown thermal coefficients with additional boundary information.|Phase-change-based parameter identification and additional-data uniqueness results have a long history. A novel theorem must exceed a rederivation of this general principle. The PDF is image-based; bibliographic scope is verified, not every theorem hypothesis.|

This is a bounded primary-source audit, not proof that no earlier glacier paper has used the proposed reference principle. The named observational application remains a hypothesis requiring an explicit priority comparison before publication.

## 1. Phase-switch ice-reference calorimetry

### Exact structural starting point

For a homogeneous debris column with imposed temperatures and an always-melting ice side, the transformation

(k,C,Rc) → (a k,a C,Rc/a),  a>0,

preserves diffusivity k/C, contact length kRc and the entire temperature field, while scaling heat/melt flux by a. Thus a unique best temperature fit need not identify absolute melt. Multiplying residuals by an adjoint-derived weight cannot make this unobserved scale appear.

In a purely cold, linear, initially known semi-infinite ice substrate with independently known effusivity ei=√(kiCi),

Zi(s)=1/[ei√s],  TI(s)=Zi(s)q(s).

With a massless contact and downward debris heat q=−kTd′,

ell(s):=−Td(s)/Td′(s)=kRc+kZi(s).

For two distinct positive real Laplace parameters,

k=[ell(s1)−ell(s2)]/[Zi(s1)−Zi(s2)],
Rc=ell(s1)/k−Zi(s1),  C=k/alpha.

One complex Fourier frequency can algebraically supply two real constraints if phase and the full linear response are known; two frequencies are sufficient, not a new minimality theorem. Multiple frequencies/windows are needed for a meaningful noisy-data consistency check.

For a finite ice column with deep temperature fixed at0°C,

Zi(s)=tanh[Di√(s/alphai)]/[ei√s].

The reference scale comes from independently known ice properties. If ice conductivity/capacity are also scaled freely, the ambiguity returns. `reference_admittance_basis.py` verifies these statements algebraically; five tests pass, including exact synthetic two-frequency recovery and the restored ambiguity when both media are scaled.

### What is genuinely difficult and potentially useful

The prospective contribution would be a **phase-conditioned absolute-property estimator with an explicit observability certificate**, tested against the usual assumed-C conversion and melt-fitted gain. It must establish that natural cooling, rather than imposed laboratory heating, contains enough information to fix the scale in a heterogeneous glacier setting.

The displayed transfer relation does not apply to a Fourier transform formed by deleting all melting samples. Ice retains history: cutting a time window changes initial-state terms, and melting outside the window changes heat delivered to the ice. Use a full causal state model or its exact convolution with initial-profile terms. Do not treat the cold-window sequence as a new periodic steady experiment. Unknown ice initial state, probe offsets, interface side, contact resistance and evolving geometry are nuisance parameters, not automatically known constants.

### Himalayan data gate

Rowan's NG2 record has11 sensors through2m of debris, six-hourly sampling and±0.25°C stated accuracy; winter basal temperatures are clearly below melting. However its uppermost sensor has a documented summer/winter bias relative to the internal profile. Use a verified interior computational boundary with its matching reduced depth or model that discrepancy explicitly. The winter basal gradient is small relative to adjacent-sensor accuracy, so parameter observability must be demonstrated, not inferred from negative temperature alone.

A sensible source window is cooling after a documented summer melting plateau, because it constrains recent near-interface ice history better than arbitrarily starting in midwinter. Even then, profile/deep-boundary sensitivity must be carried. Per-probe systematic offsets should be modelled as persistent nuisance offsets within source accuracy rather than assuming thousands of independent errors that average away.

Pirámide currently fails this natural-reference field gate: in the672 TRAIN half-hours, deepest DT5 reaches only−0.006°C, never below−0.4°C;352 observations lie within±0.4°C. A model-generated cold branch is not an independent observation of the required reference excitation.

### First pass conditions, without opening melt targets

- Demonstrate gauge recovery in synthetic data and failure in the always-melting control.
- Identify a real, source-supported cold interval and interface/probe mapping; propagate offset/depth/initial-state uncertainties.
- Show that conductivity/capacity ranges contract materially relative to their independent prior ranges because of the cold observations; report rejected and unconstrained cases.
- Verify the result across withheld temperature windows and geometry/reference-depth families before touching another melt target.
- Require a genuine future independent melt or heat-flux test for predictive claims; a temperature-only fit cannot validate its own inferred absolute flux.

## 2. Goal-oriented, gauge-aware melt inference

The immediate fix to the objective is to target cumulative basal heat, not temperature RMSE as a surrogate. But the estimator must first check whether that quantity is informed by the measurements.

Let F(theta) be the observed-temperature forward map and Q(theta) the cumulative latent heat. Locally, an admissible direction v with DFv=0 but DQv≠0 makes Q unidentifiable from those data. The warm scaling direction is an explicit example. An adjoint computes DQ efficiently; it does not eliminate that null direction.

A technically concrete method is to build a source-constrained set

S={theta: all temperature residuals lie in declared sensor/model/numerical bounds; energy/passivity and source parameter constraints hold},

then compute Qmin and Qmax over S. The midpoint minimizes worst-case absolute error **conditional on the true model belonging to S**. If the feasible flux interval is too broad for the intended application, issue a quantitative abstention and identify the additional measurement that constrains the harmful direction. If a cold-reference observation removes that direction, quantify how much the flux interval contracts.

The potentially publishable difference is a glacier-specific **melt-information criterion**, not new generic goal-oriented inversion. A fair comparison must include Laha's Bayesian/layered approaches, Beck–Nicholson's sensor-error analysis, a standard inverse heat-flux method and the fixed simple temperature baselines. Do not substitute Gaussian posterior confidence for unknown source-error semantics. This direction can give useful credible flux intervals even when it cannot improve a point forecast.

## 3. Seasonal property transfer constrained by phase physics

A second-season or cross-regime test could combine the cold-reference estimate with a conserved pore-water/ice enthalpy state. Its specific scientific question is whether an absolute scale learned during cooling transfers into melt prediction when pore contents change phase.

Frozen and wet-summer capacity are different combinations of solid, water, pore ice and air. One frozen effective capacity does not uniquely identify those constituents. A spring zero-curtain may supply additional information, but snowmelt infiltration, runoff, latent pore heat and unknown water stock can be confounded. Independent moisture/porosity/lithology constraints are needed; a richer state model is not a license to fit arbitrary seasonal corrections.

The novelty would have to lie in the identifiable cross-regime inference/validation system. The phase-change equations themselves are already in Giese/Collier. This is therefore a **conditional second milestone**, not a reason to open sealed ablation targets now. First prove a real information gain in Direction1; otherwise adding seasonal hydrology only expands the ambiguity.

## Decision

The strongest narrowly differentiated lead is **natural phase-switch reference calorimetry with an observability gate**, followed by goal-oriented flux intervals. Its mathematical ingredients are established, but the proposed glacier use and rigorous failure/eligibility map may be a worthwhile methodological contribution if the Himalayan data support the needed signal.

Proceed only with a bounded synthetic and temperature-only Himalayan observability calculation. Keep all previously scored melt periods as development, retain their failed superiority verdict, and keep unused ablation targets sealed. Do not launch a broad production campaign or claim a new theorem, absolute field validation or a Nature-level result before this gate passes.
