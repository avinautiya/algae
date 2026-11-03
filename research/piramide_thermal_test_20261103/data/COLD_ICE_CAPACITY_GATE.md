# Cold ice, wet debris and independent heat-capacity gate

No camera outcomes were opened in this audit. This is a feature-only physical-consistency check for the predeclared training interval, not a fit to melt or a new superiority result.

## Primary literature

[Giese et al.2020](https://doi.org/10.5194/tc-14-1555-2020) already represents moisture transport, water/ice phase changes and subfreezing glacier ice in ISBA-DEB. Accounting for moisture only by changing conduction is insufficient to reproduce all its energy effects. [SHAW-Glacier2022](https://doi.org/10.3389/feart.2022.796877) also treats coupled heat/water transport and phase change. These are strong physical comparators, not proof that our simpler model includes their processes.

The final [Pirámide2025 paper](https://doi.org/10.1029/2025JF008360), section2.3.3, describes transient debris diffusion plus conduction in a2m ice layer with an isothermal lower boundary. Its source conductivity inference is method-dependent; optimized conductivity uses ablation. Thus adding cold-ice storage/admittance repairs our physical scope but is not historically novel. Native in-debris measurements do not independently measure an initial glacier-ice temperature profile below the interface.

## Training-only moisture and sensible heat capacity

`training_capacity_audit.py` reads only source-bound process CSVs. It uses [2December03:00,16December03:00), the672 half-hourly training samples. The physical feature rule is declared in the script: flag theta<0 or theta>assumedphi and exclude such samples from a proxy mean, never clip. All672samples at all3probes pass for this interval. This does not make all observations elsewhere valid.

| Probe | Depth | Mean volumetric water-content reading |
|---|---:|---:|
| SM1 | 0.30m | 0.208972917 |
| SM2 | 0.40m | 0.212370685 |
| SM3 | 0.50m | 0.050192560 |

The equal-probe arithmetic mean0.157178720 is a **proxy**, not a measured volume average: the upper30cm has no moisture probe, depths are uneven and probe calibration/representativeness remains uncertain. Source porosity0.44 is assumed in calibration code, not independently measured here.

For fixed-phase liquid water/air pores, the volume-additive sensible capacity is

C=(1−phi)rho_r*c_r + theta*rho_w*c_w +(phi−theta)rho_a*c_a.

Using independently decoded lithology rho_r=2684.933333kg/m³, c_r=796.866667J/kg/K, assumedphi=.44, approximate water rho*c=4.18MJ/m³/K and air rho*c=1206J/m³/K gives:

| Conditional scene | theta | C MJ/m³/K |
|---|---:|---:|
| Dry, air-filled pores | 0 | 1.198670 |
| Training equal-probe proxy | 0.157178720 | 1.855487 |
| Fully liquid-water-saturated pores | .44 | 3.037339 |

These are three fixed finite sensible-capacity scenes, not confidence bounds on locally measured bulkC. The dry case conflicts with observed deeper moisture if interpreted literally for the whole layer; it is a sensitivity endpoint. The proxy scene is an explicitly conditional nominal assumption. Do not claim that these values solve independent heat-capacity identification.

## Why phase change cannot be hidden inside a finiteC

Liquid water and pore ice carry latent energy. A bulk enthalpy law must include a phase fraction, e.g. H(T,w)=H_sensible(T,w)+rho_w*Lf*theta_liquid, with water-mass/volume bookkeeping consistent when ice density differs. For unknown pore phase, sensibleC alone does not bound dH/dT at freezing. Atphi=.44 the potential latent term for water volumephi is146.96MJ/m³, far larger than a one-kelvin sensible adjustment. It cannot be represented by arbitrarily enlarging constantC or silently clipping negative temperatures/fluxes.

Native training temperatures contain negative values in37/18/1/11/5samples atDT1..DT5. Deep small negatives may lie within the stated±0.4°C sensor accuracy, so those values do not alone prove actual pore freezing. Nevertheless liquid-water-only fixedC is an assumption to report, not a demonstrated phase state. Cold-glacier-ice storage fixes a different energy reservoir from freezing within wet debris.

## Conservative ice-boundary framework

For a semi-infinite stationary ice column with diffusivityalpha_i, initial uniform temperatureTi0 and imposed interfaceTb(t), the conductive flux into ice is

q_i(t)=k_i/sqrt(pi*alpha_i) * [(Tb(0)−Ti0)/sqrt(t) + integral_0^t Tb'(tau)/sqrt(t−tau) dtau].

This is the classical heat-equation/Duhamel boundary admittance, not a new theorem. With q_d positive downward from debris, a melt-only interface requires

Tb<=0; m_dot>=0; Tb*m_dot=0; q_d−q_i=rho_i*Lf*m_dot.

When Tb<0, m_dot=0 and the flux cools or warms ice. When Tb=0, the residual after ice sensible storage may become latent melt. This gives a legitimate path for negative flux without inventing refreezing of absent liquid water. It does **not** independently establish the initialTi0, finite/deep boundary, pore-water phase, or moving-interface/advection effects at Pirámide. These are conditional model inputs/gates.

A finite-volume enthalpy alternative may be easier to conserve and verify. It should check a zero-forcing rest state, cooling/reheating with zero net melt before cold content is repaid, a known semi-infinite heat solution, and an integrated sensible-plus-latent ledger. Passing these tests verifies implementation only; melt accuracy still needs untouched camera evaluation.

## Decision

Use the source-independent proxyC only as a labeled conditional nominal case and dry/saturated scenes as declared sensitivity endpoints. Freeze all scenes before reading camera outcomes; do not select the winningC after evaluation. Require the cold-ice solver to preserve cooling energy, and carry wet-debris phase uncertainty explicitly. If the model cannot pass its own physical gate, preserve that failure rather than unsealing camera targets to tune around it.
