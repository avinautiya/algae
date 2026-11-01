# Equal-volume avalanche timing and the rescue window

Focused priority check, 31 October 2026. No Claude work changed.

## Question checked

Can equal total mantle-debris delivery, supplied steadily versus in a short avalanche pulse, produce different subsequent stream–cliff contact under equal water and heat budgets, because deposition shuts the toe route and thermal incision of the bypass subsequently raises the discharge needed to reconnect it?

## Finding

No exact study of this particular supraglacial toe/bypass, externally supplied debris-pulse and thermal-rescue experiment was located in the targeted searches. That is a candidate opening, not proof of historical priority. The new result must isolate a quantitatively consequential **thermal widening of the rescue barrier**. Merely demonstrating same totals/different outcomes or sediment-induced avulsion is already covered by much older river morphodynamics.

### Closest timing objection

[Humphries et al. (2012), Experimental evidence for the effect of hydrographs on sediment pulse dynamics in gravel-bedded rivers](https://doi.org/10.1029/2011WR010419); [author PDF](https://www.sfu.ca/~jvenditt/publications/2012_2011WR010419_Humphries_etal.pdf). Primary publisher full-text results/discussion inspected.

Experiments held released water volume constant and found different sediment-pulse evolution under steady versus differently peaked hydrographs, together with transport hysteresis and reach-topography effects. Therefore equal cumulative inputs do not imply equal response is established general river physics. This study does not have an ablating ice substrate or a cliff-survival output.

### Closest coupled ice/sediment objection

[Beaud et al. (2018), Modeling sediment transport in ice-walled subglacial channels and its implications for esker formation and proglacial sediment yields](https://doi.org/10.1029/2018JF004779). Primary abstract and model framing inspected, not all equations. Its semicircular bedrock-floored subglacial channels evolve through ice melt opening, creep closure, sediment deposition and removal. Thus combining ice melting and sediment conservation is already established; the open-surface toe/bypass geometry and a testable cliff-contact rescue prediction are the potential distinction.

[Ng, Canals under sediment-based ice sheets](https://www.cambridge.org/core/journals/annals-of-glaciology/article/canals-under-sedimentbased-ice-sheets/651DE0F80AD02F6409555108FAE6D66E). Primary abstract inspected. Its canal geometry already depends on coupled sediment and water fluxes with melting/deformation. A broad first-coupled-ice–water–sediment claim would be false. Publication date should be checked from the PDF rather than the website's later online-posting date before citing a year.

Other direct obstacles and inspected sections are in `prior_art.md`: alluvial bifurcation stability and avulsion duration; supraglacial thermal incision; cliff maintenance and reburial after avulsion.

## A useful exact post-loss limit

This is a reference calculation, not a claim of novel mathematics. Let the bypass follow an independently justified free-overfall law

`Q_b=C_b [H−eta_b]_+^(3/2)`.

The toe is dry when upstream stage `H` is below its hydraulic bed `eta_t`. With only the bypass active, an imposed inflow must exceed

`I_rescue(t)=C_b [eta_t(t)−eta_b(t)]^(3/2)`

to rewet the toe. Equality permits a zero-depth limiting state and is not positive toe discharge.

For `g=eta_t−eta_b>0`, the bed budgets give

`gdot=D_t/(1−p_t)−E_t+E_b`

when the dry toe has no sediment export, bypass sediment elevation is unchanged, and the variables use consistent vertical coordinates. `D_t` is solid-volume deposition per plan area, `E_t` the lowering beneath the toe stock, and `E_b` bypass ice incision. Consequently

`I_rescue(t)=C_b [g0 + integral(D_t/(1−p_t)−E_t+E_b)dt]^(3/2)`.

If `E_b>E_t` after the deposition pulse ends, the rescue discharge can keep growing even when the added sediment volume stops growing. A bounded flood class `I<=I_max` cannot rescue after `g>(I_max/C_b)^(2/3)`, provided the assumed hydraulic regime persists. This creates a finite operational rescue window. It does not prove permanent abandonment under arbitrarily large future floods or changing geometry.

For a time-varying coefficient `C_b`, changing geometry, frictional backwater or downstream control, recompute the hydraulic threshold instead of retaining this weir formula. Reservoir storage affects when the stage reaches the threshold; use the volume budget to test finite flood duration too.

## Decisive counterfactual

Keep identical complete forcing histories and debris histories. Disable only the **difference in substrate incision** after toe loss, then compare `I_rescue(t)` and actual finite-flood reattachment. If the effect remains equally large on a frozen substrate, it is predominantly ordinary sediment avulsion, not the proposed thermally growing rescue barrier.

Use a second comparison with the same debris total spread steadily, but preserve its source composition and heat capacity. A rapid deposit can generate an artificial impulse if its placement displaces water without conserving volume. Count displaced water, pore volume, work and heat explicitly.

The endpoint should be a predeclared reattachment criterion over a finite horizon, time-to-loss and cumulative cliff-maintenance heat. A model without an evolving exposed face must report the last quantity as a proxy, not claim a validated cliff lifetime or glacier mass-loss improvement.
