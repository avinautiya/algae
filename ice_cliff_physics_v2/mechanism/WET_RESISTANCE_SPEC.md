# Frozen wet-boundary comparison: hydraulic contact is not thermal delivery

This is a constitutive-model comparison, not a discovery that debris insulates ice. No outcomes from the parent's new resolved-slab or coupled campaign have been read by this agent. The small exact reference formulas below define the expected correctness result before those runs.

## Primary outputs and exact reference

For one wet-covered patch, F=rho cp Q, G=h Awet and K=k Awet/d. Its finite-water conductance is H=F[1-exp(-G/F)]. A prescribed nonnegative net absorbed source qa acts on the debris, the ice interface is0°C and water enters at Tin>=0. Under steady conditions,

    Ts=(qa+H Tin)/(K+H),
    qice=K Ts,
    qwater=H(Tin-Ts),
    Tout=Tin-qwater/F.

Water exchange is signed: hot externally heated debris can warm the outgoing stream. The steady dry-column ice flux under the same imposed source is qa, so

    qice,wet-qice,dry = H/(K+H) [K Tin-qa].

Rewetting reduces local basal melt if qa>K Tin, even when incoming water is above freezing. This follows from a passive resistor/advective heat balance; it is not a new fundamental theorem or a proved first occurrence in glacier literature. It is independent of a warm initial condition.

A wet-treated-bare proxy with the same finite water supply delivers qa+H Tin directly to ice. Its overprediction factor is1+H/K. In the zero-source case, the steady thermal-reconnection ratio relative to bare wet ice is1/(1+H/K). To retain at least fraction r of bare heat delivery, thickness must satisfy

    d <= k Awet (1-r)/(r H).

This is an explicit conditional thickness forecast in the declared steady one-dimensional closure, not a field-validated reconnection deadline.

## Twelve fixed cases

`wet_boundary_cases.json` varies film/solid conductance ratio G/K={0.1,1,10}, wet fraction={0.25,1} and source ratio qa/(K Tin)={0.5,2}, giving12cases. Fixed thickness0.05m, conductivity1W/(m K), volumetric heat capacity1.5e6J/(m3 K), total covered area0.2m2 and inlet superheat0.5K. Water capacity is set so NTU=G/F=0.01 in every case, explicitly an engineering test network, not a measured glacier stream. The dry and wet columns have the same material properties. An independently documented environmental case can be added later; these nondimensional cases must not be labelled a field distribution.

Use zero initial debris/ice enthalpy at0°C and the SAME prescribed absorbed source per area before/after rewetting. Do not remove atmospheric energy from a patch simply because it becomes wet; that would change two causes at once. Conversely, the equal-source policy is a controlled boundary experiment, not an assertion that actual absorbed sunlight stays unchanged under turbid water. A realistic application must supply radiative/advective source partition or show attenuation/zero-source sensitivities. Source energy, latent ice energy, debris storage and outgoing-water energy are all retained.

## Compare four calculations

1. **Reference:** independently discretized1-D diffusion slab(s), separate wet/dry columns, temperate-ice lower boundary and signed finite-water upper boundary. Verify space/time convergence and steady reference before using it as a comparator.
2. **Production approximation:** finite-capacity lumped wet model with otherwise identical inputs. Report transient error separately from steady resistance error; G/K>>0.1 is not automatically a valid lumped-capacity regime.
3. **Wet-treated-bare proxy:** prescribe direct finite-water plus supplied external heat delivery on the same covered footprint, bypassing debris resistance. This is an intentionally simplified proxy, not an assertion that every professional model uses it.
4. **Dry-only port:** same debris material, heat capacity and external source, with water coupling omitted. If the original model refuses wet application, record its abstention. A forced zero-water evaluation is an explicitly misapplied-port diagnostic, not an authentic professional-model prediction.

Dry/wet columns must not share one instantaneous temperature without a tested lateral-mixing approximation. Maintain extensive enthalpy when wet area changes; test wetting/drying remapping separately from the constant-area constitutive campaign.

Run steady approach for a predeclared20diffusiontimes(Cv*d²/k), or use an exact steady solve. Check final ice power against the formula to0.5%, and require a refinement bracket smaller than20% of each reported comparator error. Transient cumulative-melt errors need their own time refinement. No run may be removed based on its direction of error.

## Report, without conflating targets

Primary: mean heat and melt rate over the MATCHED wet-covered footprint Awet. Secondary: cumulative melt over the complete covered area, exposed bare face and all ice separately. Hydraulic connection is supplied in this test and is not melt recovery. Include qa>K Tin, equality and qa<K Tin analytic checks, positive flow, signed water heat, water-temperature rise and every ledger residual.

A >10% proxy error outside the numerical bracket establishes an error in that proxy for that imposed regime; it does not establish environmental prevalence or novelty. Any total-cliff or glacier-mass-balance claim additionally needs actual wet-covered area, material, heat-source and flow distributions. For small patches at field discharges, exported heat may warm the stream by much less than±0.05°C; existing sensor resolution cannot automatically validate this term.

## Priority and baseline boundary

Buri et al.2016 already contains water-driven extra melt and reburial, and SHAW-Glacier2022 already includes advective heat/water transport. Our candidate is a quantified cliff-boundary error/criterion and a future coupled application, not first wet-glacier physics. Inspect `WET_BOUNDARY_PRIOR_ART.md` before writing a headline. No publication or competition outcome follows from this verification campaign.
