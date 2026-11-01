# Draft methods insert — environmental extension

Use after integration and adjust verb tense to the actual completed runs.
The numerical examples below are **synthetic verification**, not environmental
results. This text does not replace a full coupled validation section.

## Distinct geometric observables

Let a vertical cliff at horizontal position `x` rise from channel-floor elevation
`f` to crest `c`, with extrusion width `W`, left boundary `a`, right ice boundary
`b` and basal elevation `zb`. Its nonconvex stepped cross-section has ice volume

`V = W [(b-a)(f-zb) + (b-x)(c-f)]`.

Face-normal retreat and floor lowering are independent removals. Applying net
latent energies `Ef` and `Eb` successively gives

`dx = Ef/[rho_i L W(c-f)]`,

`-df = Eb/[rho_i L W(x+dx-a)]`.

The two removed regions are disjoint and their total volume is exactly
`(Ef+Eb)/(rho_i L)`. The independent polygon-area check verifies this equality.
The synthetic 3 MJ example removes 0.0097950228 m³; its polygon/energy volume
residual is approximately 1e-15 m³. This is an algebraic/component check, not
an accuracy estimate for a glacier. The model homogenizes face retreat and
holds the crest fixed; it does not resolve undercutting, slumping or fracture.

Wet bare face, dry bare face and debris-covered face are partitioned without
overlap. Heat on these faces is kept distinct from channel-floor heat. A mask
partition by itself does not establish a spatially varying face-evolution law.

## Environmental energy ports

The atmosphere helper evaluates projected/shaded absorbed shortwave, net
longwave and explicitly parameterized neutral bulk sensible and latent exchange.
All meteorology, surface coefficients and physical constants are required
inputs. Air temperature is never substituted for water temperature. Terrain
radiation/shading must be resolved upstream rather than inferred by this helper.

A fixed-area debris column obeys

`Cd (Td_new-Td_old) = dt [qatm - G Td_new]`.

The reference ice interface is at the melting point. Conductive heat plus
separately delivered water heat first repays a nonnegative cold reserve `Cc`;
only the remainder is assigned to latent melting. Cooling can increase `Cc`.
The code reports

`Eatm + Ewater = DeltaEdebris + DeltaEice_sensible + Elatent`.

A direct bare-ice branch omits debris storage. These are lumped energy states,
not a resolved subfreezing ice-temperature field. Departed meltwater cannot
refreeze in this model. Evolving area needs an explicit remapping of storage;
fixed-area states cannot simply be multiplied by changing area without an
additional energy term.

Serial water heat exchange prevents assigning the same advected heat to the
floor and cliff face. For water heat-capacity rate `Cw=rho_w cp Q`, the floor
receives `Cw Tw [1-exp(-UA_floor/Cw)]`, and the face exchanger receives heat from
the remaining downstream sensible energy. This is a plug-flow closure at a
melting-temperature boundary. It is not a resolved turbulent heat-transfer law.

## Transmission and applicability

Parallel ratings divide the physical corridor width once. Each route has an
actual ice floor, a debris crest and an explicitly parameterized aperture or
effective Darcy conductance. No path admits water below intact ice. Hydrostatic
pressure is integrated over the vertical aperture; a separate fully wetted
narrow-slot closure includes viscous and entrance losses through

`DeltaH = [12 mu ell/(rho g w a^3)] Q + Q^2/[2g(Cd w a)^2]`.

The slot model requires both ends above its roof and smooth laminar narrow-gap
conditions. Thin-opening results are not transferable to long microscopic
under-rock gaps without assessing friction, development and roughness. These
closures define conditional alternatives, not a measured site transmission.

Depth/length, Froude number and curvature acceleration `u²/(gR)` are reported
without converting them to an automatic hydrostatic-validity pass. Unknown
curvature remains unknown. A resolved crest comparison is still required.

## Observed geometry and the decisive prediction

The [Lirung/Langtang archive](https://doi.org/10.5281/zenodo.13073856) supplies
83 profile epochs across 9 cliff identifiers. The analysis includes every
profile, with fixed binning and explicit survey-error handling. These data
provide observed geometry, not co-located hydraulic or thermal forcing.
Profile points are not independent cliff tests, and profile change is not
automatically attributed to melt.

The testable extension is improved prediction of common-support **face melt**
from dynamic hydraulic recovery, relative to prescribed contact, under matched
environmental energy and transmission assumptions. Previous work already
connects stream hydrology with cliff survival and reburial:
[Petersen et al. (2024)](https://doi.org/10.5194/esurf-12-727-2024), while
[Kneib et al. (2022)](https://doi.org/10.5194/tc-16-4701-2022) documents temporal
cliff evolution. The proposed specific approximation error must be tested;
these new components alone do not establish its magnitude, priority or field
predictive benefit.
