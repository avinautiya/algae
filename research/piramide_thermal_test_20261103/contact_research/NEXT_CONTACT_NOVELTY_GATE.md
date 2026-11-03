# Basal thermal-contact-resistance novelty and plausibility gate

**Decision: NO production contact-resistance campaign now.** A finite debris–ice contact resistance is legitimate established physics, but the present evidence neither establishes its originality in this application nor supports it as the sole explanation of the approximately4× raw melt-amplitude error. A bounded analytical/process falsifier is justified; another sequence of variants tuned against the four already inspected melt periods is not. Tower2 camera targets remain unopened.

This decision does **not** prove that every spatially variable, state-dependent or intermittently wetted contact model is impossible. Those richer hypotheses would need independent evidence and a new unused test; their flexibility is not itself evidence or novelty.

## 1. Closest primary work and what it does

|Primary source|Relevant treatment|What this does and does not establish|
|---|---|---|
|Hammonds & Baker2016, [*Investigating the thermophysical properties of the ice–snow interface…PartII*, DOI10.1016/j.coldregions.2016.01.006](https://doi.org/10.1016/j.coldregions.2016.01.006); [open author ISSW proceedings record](https://arc.lib.montana.edu/snow-science/item/2273)|A porous snow–solid ice interface has limited contact pathways; the authors attribute sharp near-interface thermal gradients primarily to contact resistance using a1D energy analysis.|Direct prior art for an explicit porous-medium/ice contact mechanism. It concerns snow, not supraglacial rock debris. Publisher abstract and author proceedings description were read; the19.5MB proceedings PDF exceeded the browser parser limit. No unread numbered equation is represented as verified.|
|Amschwand et al.2025, [Murtèl rock-glacier heat fluxes, DOI10.5194/esurf-13-365-2025](https://esurf.copernicus.org/articles/13/365/2025/), §3.2.2 and AppendixA|A measured coarse-block cavity exchanges radiative/convective heat. The heat-flux-plate correction explicitly uses a series circuit, Rtot=Rrock+RHFP+Rs, with Rs bounded by an inverse radiative conductance.|Closest inspected coarse-block/cavity analogue. The cited contact resistance is partly an **instrument-contact correction**, not an independently established natural debris–ice Robin coefficient. Its cavity measurements cannot be transferred quantitatively to the Pirámide basal interface without geometry evidence.|
|Garcia & Santamarina2021, [*Heat Flow in Fractured Rocks*, DOI10.1016/j.geothermics.2021.102113](https://doi.org/10.1016/j.geothermics.2021.102113)|Experiments vary stress and moisture in limestone contacts. Apparent temperature discontinuities arise from constriction and low-conductivity gaps; water bridges improve heat transfer and reduce contact resistance.|Direct prior art for roughness, contact pressure and wetting dependence. It does not measure the requested glacier basal coefficient. No rock-fracture stress law should be transplanted to an unconsolidated debris blanket without justification. Publisher primary description/methods were inspected; a complete equation-level numerical law was not accessible.|
|Reid & Brock2010, [DOI10.3189/002214310794457218](https://doi.org/10.3189/002214310794457218)|The debris diffusion model fixes its ice-contact base at melting temperature during ablation, supported by that study's measurements. Heat reaches ice through the debris conductivity/thickness.|A clear Rcontact=0-style reference boundary in a mainstream model. Its assumption is not an unrecognized universal error, and successful local justification cannot be disproved solely by another dataset.|
|Hénot et al.2022, [glacier tables, DOI10.5194/tc-16-2617-2022](https://tc.copernicus.org/articles/16/2617/2022/), Eq7; related [PRL2021 DOI10.1103/PhysRevLett.127.108501](https://doi.org/10.1103/PhysRevLett.127.108501)|Heat transmitted under a rock is represented by Qrock→ice=λ(Trock−Tice)/d1D, with an experimentally constrained effective path length and explicit contact-area geometry.|Existing rock–ice finite-conductance and geometric contact-area physics. It does not separately identify a basal temperature-jump coefficient. Merely adding a resistance to this circuit is not new mathematics.|
|Melo-Velasco et al.2025, [final JGR article DOI10.1029/2025JF008360](https://doi.org/10.1029/2025JF008360), §§2.2.2–2.3.1,4|Conductivity methods, sensor-depth uncertainty and material assumptions produce different melt predictions; optimized effective parameters differ from field-derived values.|The general bulk-versus-field discrepancy is already the paper's central result. It cannot be reintroduced as our new discovery. Warm near-base readings alone do not identify contact resistance.|
|Petersen et al.2022, [DOI10.1029/2021JF006520](https://doi.org/10.1029/2021JF006520), §6.3|Thermistor-profile hysteresis reveals nonconductive contributions; fitting one effective diffusivity need not meaningfully absorb them.|A substantive competing explanation for diffusion-inferred parameters that fail to predict mean melt. Nonconductive transport cannot be declared newly discovered from this project's residual.|
|Giese et al.2020, [DOI10.5194/tc-14-1555-2020](https://tc.copernicus.org/articles/14/1555/2020/)|Moisture/phase changes modify heat capacity, conduction and evaporative/sublimative losses in debris.|Existing porous wet-debris physics. Bulk moisture is not a direct measurement of basal wet-contact fraction or water-film thickness.|

**Novelty result:** the Robin/contact law, resistance-in-series mathematics, rough-contact/water-bridge mechanism and porous ice-interface discontinuity are established. An independently demonstrated debris-specific basal resistance with transferable predictive consequences could still be a contribution, but this bounded search does not establish priority for that narrower claim. Absence from these inspected glacier models is not proof of absence from the field.

## 2. Correct mathematical hypothesis

Separate debris-side basal temperature TDΓ and ice-side interface temperature TIΓ. With area-normalized resistance Rc≥0,

qΓ=(TDΓ−TIΓ)/Rc,

qΓ−qice=L m,  TIΓ≤0,  m≥0,  TIΓm=0.

The Rc→0 limit imposes temperature continuity. A massless contact stores no energy and adds no heat source. Its passive dissipation is qΓ(TDΓ−TIΓ)=Rc qΓ²≥0. In finite volumes the cold-branch conductance is

G=[Δzd/(2kd)+Rc+Δzi/(2ki)]⁻¹.

For the melting branch TIΓ=0, the debris boundary conductance is[Δzd/(2kd)+Rc]⁻¹; **the debris-side basal temperature is not also forced to0**. The existing sensible-plus-latent ledger remains unchanged. This is standard transport/complementarity mathematics, not a new theorem.

For a linear, constant-coefficient, entirely melting column of depthd, the harmonic transfer from surface temperature to basal heat is

H(ω)=kλ/[sinh(λd)+Rc kλ cosh(λd)],  λ=√(iω Cvol/k).

At zero frequency H(0)=1/(d/k+Rc). Unlike a constant output gain, Rc changes frequency response and phase. This supplies a concrete discrimination test under the stated linearization; it does not establish that a heterogeneous, intermittently cold field column satisfies that idealization.

## 3. Identifiability before predictive testing

In local steady downward conduction, a deepest probe at unknown clearanceδ above melting ice measures

Tprobe/q=Rc+δ/k.

Consequently Rc and clearance are confounded. A temperature profile plus a scalar melt measurement does not independently prove an interfacial jump without known geometry/material properties. Allowing k(z), unknown sensor offsets, storage and nonconductive fluxes adds further confounding. A second empirical melt-amplitude gain after fitting Rc would largely reopen the original compensation problem.

**Geometry correction:** the2025 paper says the ablation stake was drilled in the same pit as the thermal/moisture sensors to keep debris thickness consistent (§2.2.3). The3m PVC measurement is stake length. Table1's tower-to-sensor distance is not a demonstrated stake-to-thermistor separation. Therefore a supposed3m thermal/melt sampling mismatch is not an established cause. Sensor placement/contact, ±2cm depths, settling and unresolved lateral variability remain hypotheses rather than confirmed explanations.

## 4. Quantitative TRAIN-only plausibility screen

`contact_screen.py` reads only already released Tower1 TRAIN calibration and temperature data. No Tower2 camera or new held-out melt value is read. Reproducible output is `TRAIN_CONTACT_SCREEN.json`.

For Dec2–16, source-label time-weighted means are DT1=11.4274°C and DT5=0.37380°C. The observed TRAIN latent-equivalent flux is23.5051W/m². The current primary uses k=4.38900W/mK across a0.5m computational column, giving bulk resistanceRd=0.113921m²K/W. Its empirical output factor is0.235027.

A **steady series-resistance explanation of that entire amplitude factor** would require

Rc,required=Rd(1/g−1)=0.370793m²K/W,

and an average contact temperature jump qRc≈**8.716K**. In contrast, the warm deepest probe is only0.374°C. Under the restrictive steady/temperate/colocated assumptions, even taking its temperature0.4°C higher and setting clearance to zero gives Rc≤**0.032921m²K/W**. Positive clearance would reduce the inferred contact resistance further.

This is a strong warning against assigning the entire4× mismatch to one constant basal resistance while retaining the observed near-base profile. It is **not a transient proof against every contact model**: net sensible storage, a cold ice-side temperature, sensor error/geometry, lateral transport and phase changes can invalidate the simple mean-flux inference. Those terms must be quantified rather than assumed away. The empirical factor came from the present model's full transient output, so its steady Rc equivalent is explicitly a diagnostic approximation.

A further **declared sensor-profile proxy**, not a measured total-energy bound, integrates piecewise-linear initial/final TRAIN temperatures with the proxy capacity and a0°C base. It gives debris sensible storage increase1.189MJ/m², or0.983W/m² averaged over TRAIN; the conditional±0.4°C sensor contribution is±0.614W/m². This is small compared with the roughly76W/m² raw-versus-observed latent-flux mismatch. `storage_proxy.py` reproduces `TRAIN_STORAGE_PROXY.json`. Unknown ice storage, pore latent heat, advection and evolving geometry are not thereby bounded or excluded.

A pure-conduction air-gap analogy would translate the required resistance to about9.27mm at an illustrative air conductivity0.025W/mK. That analogy omits radiation, convection and solid/water bridges. For two broad facing surfaces near273K with illustrative effective emissivity0.9, the linearized radiative conductance alone is about4.16W/m²K; that would limit resistance to about0.24m²K/W even before parallel conduction. Actual view factors/contact area are unknown, so this is another plausibility check, not a measured glacier bound. A water-filled gap is generally a better conductor than the same air gap; moisture cannot simply be assigned arbitrarily large insulating contact resistance.

## 5. Alternatives that must be compared first

- **Depth/clearance error:** fit or bound sensor positions only within independent source uncertainty. The near-base mean temperature can arise from an ordinary finite clearance without a discontinuity.
- **Layered material properties:** one fitted diffusivity times one bulk capacity need not represent mean vertical conductivity. The observed TRAIN mean profile has substantially different gradients between depth pairs. Layered inference is already in the cited literature; it is an established comparator, not automatic novelty.
- **Nonconductive transport and wet-pore phase heat:** existing literature predicts deviations from a single diffusion coefficient. A basal resistance should not absorb all unmodelled distributed processes without a local interface observation.
- **Geometry evolution/sensor movement:** gravity-driven settling and changing ice/sensor relations require observation, not an assumed post-score correction.
- **Representative-volume mismatch:** local probe temperatures and area-integrated energy can differ in coarse debris. Same-pit placement mitigates a simple location mismatch but does not establish1D bulk representativity.

## 6. Go/no-go rules

**NO-GO now:** do not launch another production or Tower2-target campaign justified by a claim that contact resistance is unprecedented or already explains the amplitude error. Neither claim passed this gate. Do not tune Rc to the four inspected melt periods and then call their improvement independent validation.

**Permitted small next step, if pursued:** a TRAIN-temperature-only/profile-consistency and identifiability calculation, with independently bounded geometry, explicit storage corrections and the established R=0/layered/nonconductive alternatives. No additional camera target needs to be opened. A candidate Rc must simultaneously satisfy the near-interface temperature constraints, nonnegative passivity, source material/geometry constraints and the full energy budget—not merely match one scalar melt total.

Only if this process gate passes should one freeze a single contact hypothesis and all strong baselines before exposing a genuinely unused Tower2 or external melt target. Prefer deriving Rc from independent interface observations or temperature-frequency response, with no second free output gain. A positive field result would still need uncertainty and priority review; no publication/competition outcome is guaranteed.

**Current strongest conclusion:** the frozen model's field-superiority hypothesis failed. The contact idea is a physically plausible but currently weak and confounded diagnostic hypothesis. Its simple constant-resistance explanation is already in tension with the existing TRAIN near-base temperatures. That analytical exclusion is useful quality control, not a new environmental discovery or a reason to erase earlier results.
