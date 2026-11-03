# Independent pre-execution model audit

Reviewed before actual field-temperature fitting, camera outcome access or target
scoring. Current portable model freeze SHA256:
`f7c7ff19a320060b8cdcd7a7dcbd74b970956903606ccc73d7286db254b2970a`.
All six frozen model/protocol/spec hashes and all four input hashes match.
Execute from `model/`; the four data paths are relative to that working directory.

## Operator, indexing and units

The full source selection includes both endpoints of the28-day span. It contains
1345 half-hour labels and1344 intervals. Split is label index672. Thus
`heat[:672]` integrates exactly the14-day training camera window and
`heat[672:]` exactly the14-day evaluation window. Thermal calibration excludes
the split sample: its last observation is02:30, before the03:00 evaluation start.
The state propagates through the split; no evaluation profile reinitialization
occurs. Right-hold forcing uses the shared split boundary temperature during the
last training half-hour; this is declared conditional forcing, not later melt.

For normalized storage `∑ Δz T`, flux is `α ΔT / Δz`, in K·m/s. Its time integral
is K·m. Multiplying by bulk volumetric capacity gives J/m²; dividing by latent
heat gives kg/m², numerically mm water equivalent. A fitted melt gain has units
`mmwe/(K m)`. Its implied `gain*Lf` is an **effective** capacity, not an independently
measured material property. Negative upper/debris temperatures are signed, not
zero-clipped. Both interval heat signs and every internal basal-flux substep are
checked before no-refreezing melt conversion.

DT1 is an observed buried upper computational boundary. DT5 is interior to each
declared geometry; it is not forced to0C. The ice boundary lies below it. This is
a conditional geometry/phase model, not independent verification that the actual
ice interface is at those exact coordinates. Source marginal depth bounds do not
establish a probabilistic joint geometry distribution.

## Timing and numerical qualification

A synthetic constant-boundary14-day solve using the same32cells,300s step and
four interior output depths took0.093s on the available machine. This involved
4032 steps or129024cell-steps. Six scenarios times25grid points imply about14s
of training solves; twelve full28-day method/scenario solves about2.2s more.
Allow20–60s including callback fitting, serialization and machine variability.
This is a synthetic workload estimate, not a measured field-campaign runtime.

The actual frozen grid is25points from`1e-8` to`3e-6 m²/s`, not`1e-5`.
Backward Euler is stable, but stability does not establish accuracy. At the
smallest diffusivity, diurnal skin depth is approximately0.0166m, around one
32-cell width. At the largest diffusivity, cell Fourier number is about3.7 for
the300s step. A target-independent mesh/time refinement using frozen selected
parameters is still needed before a quantitative numerical-accuracy claim.
Do not retune after camera evaluation to make refinement improve the score.

## Analysis boundary and remaining checks

The analysis uses one independent TRAIN-camera gain for each candidate or
baseline and preserves signed camera endpoint differences. Its primary evaluation
is one14-day interval from one tower. All geometry/hold scenarios must be retained.
The derivative-OLS comparator is not the paper's normalized-TLS method or a full
professional glacier model. The air-temperature baseline explicitly reconstructs
gaps up to18h; this is not fully native observed forcing. Observed evaluation
upper temperatures make the experiment a conditional hindcast, not a forecast
of future weather.

Two pre-execution linkage checks were requested of the parent analysis owner:
bind the CLI normalized-prediction argument to the exact receipt prediction SHA;
verify the DT1/air CSV hashes against the decoded process manifest before dosing.
The source/calibration model itself does not read camera outcomes. No favorable
point prediction, physical discovery or professional-model superiority is implied
by this audit.
