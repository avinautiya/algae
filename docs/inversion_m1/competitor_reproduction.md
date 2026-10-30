# Competitor reproduction (Baseline E): what was reproduced, and what the number means

**Method labelled accurately:** MAP (maximum a posteriori, multi-start L-BFGS-B, 6 starts) of the **archived** unconstrained target density, in NumPy. It is **not** NUTS and **does not reproduce posterior uncertainty**. The published posteriors (4 chains x 1000 draws) were not regenerated.

| Check | Result | Evidence |
|---|---|---|
| Scaler constants | match the raw pickle bytes (pickle not executed); internal identities hold | `equivalence_check.json` A |
| Emulator, Keras (TF 2.21) vs NumPy, 2000 inputs | max abs difference 6.0e-7 | B |
| Target density, original TFP 0.25 block (verbatim) vs NumPy, 91 jointly finite points | max abs diff 1.2e-3 (max relative 7e-7) | C |
| Spectral error over the 105 archive spectra, 350-1200 nm, fixed sigma | mean MAE 0.0025; **mean relative error 1.79 %, max 11.9 %** vs published 1.8 % and 11.8 % (n = 104; their error definition is not stated, ours is mean(\|resid\|/obs) per spectrum) | `darkzone_repro/summary.json` |
| Retrieved ranges | no spectrum at an algal, cryoconite or dust prior cap | summary |
| Optimiser agreement across starts | for 42 of 105 spectra the best and worst of 6 starts differ by more than 1 nat in log-target. **This is optimiser variability on a rough landscape, not an uncertainty estimate**; NUTS with 4 random inits reportedly converges to one mode for 95 % of their spectra | `per_spectrum.json` |

Not verified: the NUTS kernel and gradients; the training configuration of the emulator; the exact CSV the pipeline reads (resampling 1 to 5 nm unspecified); the 2.5 nm grid hypothesis (section D of the equivalence check, rms 0.002 HCRF vs sigma 0.003); the paper's 104 vs 105 spectra.

## Does the algal ambiguity persist under the competitor's observation model?
Profile of the archived target over algae (12 spectra, every 9th; 11 log-spaced values 1e2 to 4.9e5 cells/mL, spacing 0.37 dex; other six coordinates re-optimised from the MAP plus 3 random starts; 1.92-nat drop):
* with the priors: median interval 0.18 dex (narrower than the grid for most), 1 spectrum bound-limited at the lower and 1 at the upper end, maximum width 1.5 dex;
* likelihood only (no priors, no Jacobian): median below the grid step, maximum 2.6 dex, 2 lower- and 1 upper-bound-limited.
So under sigma = 0.003 and a perfect model, algae is mostly **closed by the data**, with a minority of spectra (up to about 25 %) where it is not. This conditional narrowness says nothing about correctness: no concentration measurement exists for these spectra, and Baseline A shows that equally narrow intervals can sit 0.7 dex from counts. Subset only (12 of 105), approximate threshold (a Student-t likelihood is not chi-squared).
