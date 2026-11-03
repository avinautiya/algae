# Initial-state comparison for the passive contact column

For identical forcing, material parameters and geometry, consider two solutions of the fixed-interface passive contact model with componentwise ordered initial temperatures. The backward-Euler diffusion matrices are M-matrices; the zero-temperature unilateral contact is an order-preserving obstacle. Thus the initially warmer state remains componentwise warmer, and its latent-sink rate cannot be smaller (the warm-branch candidate rate is a positive linear functional of previous temperatures, thresholded atzero).

Let delta denote hotter minus colder. Summing the discrete energy balances gives

    deltaQ_latent = deltaE_initial − deltaE_final
                    + deltaQ_surface_input − deltaQ_deep_output.

State ordering gives deltaE_final>=0. The shared prescribed surface temperature gives deltaQ_surface_input<=0, while the shared deep boundary gives deltaQ_deep_output>=0. Consequently

    0 <= deltaQ_latent <= deltaE_initial.

For arbitrary initial error boxes, compare the lower and upper admissible initial states; the largest latent difference is bounded by their sensible-energy difference. This remains a MODEL-class bound, not a proof that chosen temperature boxes contain the glacier's actual initial state. It requires common capacities/coefficients and boundary forcing. Separately recalibrating capacity changes the comparison problem, so this bound does not automatically enclose all calibrated forecasts. Moving-interface advection, liquid/refreezing mass, and unmodeled sources are outside this contact model.

For the two-metre ice slab alone, the difference between uniform−1C and0C initial states is 2*917*2100 J/m². Dividing by333700J/kg gives11.54mm water-equivalent potential latent difference. Debris initial-state differences add their own capacity-weighted bound. This can explain why an explicitly modeled initial-state sensitivity need not invalidate an entire long-season conditional comparison.

These are established comparison/passivity ideas, not a historical new theorem. The independent randomized check uses100different forcing/initial-state pairs, checks final temperature ordering, and compares actual computed latent heat against initial energy differences. A numerical check does not replace the assumptions above or a formal independent proof review.
