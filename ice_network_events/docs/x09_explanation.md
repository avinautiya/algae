# Why averaged forcing shortens outlet-2 drying in X09 (stored outputs only; X09 is development evidence)

X09 (sequence study): episode A inflow 0.0206 m³/s, episode B 0.1166 m³/s, 12 h each, 12 h continuation at the mean 0.0686 m³/s; sill gap 0.20 m; lake 1125 m². Numbers from `x09_diagnostics.json` (network r17, `atol_energy`=1e-3 J, resolved reference 200 cells).

| run | outlet-2 dry h (AB / BA / AVG) | final share outlet 2 (AB / BA / AVG) | sill 1 / sill 2 lowering, mm (AB / BA / AVG) |
|---|---|---|---|
| network r17 thermal | 12.4 / 11.0 / 0.7 | 0.125 / 0.161 / 0.279 | 124/64, 111/66, 94/121 |
| network r17 fixed geometry | 12.2 / 10.6 / 0.7 | 0.273 / 0.273 / 0.273 | 0 |
| resolved ref200 thermal | 12.4 / 11.3 / 0.8 | 0.088 / 0.116 / 0.216 | 122/52, 109/51, 93/98 |

1. **Water-storage timing (explains the wet/dry mismatch).** I*₂ = Q₁(gap) = 0.0269 m³/s; episode A (0.0206) is below it for 12 h ≫ τ_c = 1.65 h, so the stage cannot reach sill 2; the mean (0.0686) is 2.5× I*₂, so the averaged run is wet after the initial fill (0.7 h dry). The fixed-geometry run reproduces the dry durations (12.2/10.6/0.7 h) — heat and sill lowering are not needed for the mismatch. Statement [R]: comparison theorem (`event_theory.md` §1.2): inflow below I*₂ throughout episode A with the stage initially below sill 2 keeps outlet 2 dry for the whole 12 h; the dry duration is then set by the episode length, not by the mean.
2. **Thermal storage and evolving sills (explain the 12–15 pp final-share discrepancy).** With fixed geometry the long-time share is identical for all orders (0.273) — the ordering/averaging difference there is only transient. With thermal lowering, outlet 2 receives heat only while wet: AB spends the first 12 h dry, so sill 2 lowers 64 mm (AVG 121 mm) while sill 1 lowers 124 mm (AVG 94 mm); the gap then widens and the final share is 0.125 vs 0.279 (15 pp; resolved reference 0.088 vs 0.216, 13 pp). This is the wetness–temperature covariance: the mean enthalpy is delivered to a different outlet.
3. **Numerical event handling.** Wet/dry transitions are 1–3 per run (no chatter at the 0.5 % of mean-inflow threshold with 15 min debounce); the resolved reference reproduces the dry durations to 0.1–0.3 h and the share discrepancy within 2.4 pp; water/mass/energy residuals ≤ 1e-15 relative. The stall class (first overtopping with legacy tolerances) is not involved (tolerance corrected).
4. **Conclusion.** Averaging fails here because the mean inflow sits on the wrong side of I*₂ for half the horizon (water storage + threshold) and, secondarily, because heat is delivered only to the wet outlet (thermal storage + evolving sills). It is not a numerical artefact.
