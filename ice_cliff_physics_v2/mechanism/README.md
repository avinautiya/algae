# Independent wet-debris mechanism audit

This packet does not modify a solver. It provides an independent passive-network identity, a scale screen and a frozen coupled-experiment specification.

Run:

```sh
python -m unittest discover -s . -p 'test_*.py'
python run_screen.py
```

Fifteen tests pass, including independent RK4 integration, arbitrary piecewise forcing, the dry limit, the storage-to-ice factor bound and disappearance of an equal-total forcing-order effect after common thermal relaxation.

The48 illustrative scale-screen combinations give maximum whole-face equivalent lowering22.3micrometres and maximum coated-footprint lowering0.446mm for the declared terminal storage excursion bounds. These are assumed fixed-network calculations, not a bound on nonlinear evolving geometry and not field observations. They argue against treating the wet-storage kernel alone as a consequential environmental discovery. Any much larger coupled effect must demonstrate actual geometry amplification beyond numerical uncertainty.

`THEORY_AND_NOVELTY.md` derives the identity, scope and literature boundaries. `CAMPAIGN_SPEC.md`, `campaign_cases.json` and `CAMPAIGN_FREEZE.json` precede the proposed coupled outcomes. No coupled outcome was computed by this agent. Thermal storage, moisture and stream-maintained cliffs have clear precedents; no new fundamental theorem or historical priority is claimed.

`WET_RESISTANCE_SPEC.md` and its separate freeze define12 steady wet-boundary cases; `wet_resistance.py` supplies exact finite-water/resistance and signed heat-export predictions, plus a conditional thickness limit. `write_steady_references.py` generates the analytic reference CSV. These are expected constitutive limits, not independent observations or a new-physics discovery. Parent resolved/coupled results remain separate.
