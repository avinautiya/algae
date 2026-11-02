State stores total joules over positive area: signed debris thermal energy
(zero at 0 C), nonnegative cold debt, and nonnegative cumulative latent energy.
The conserved accounted quantity is debris energy minus cold debt plus latent.
Latent heat is explicit for converting newly melted energy to kilograms;
the state itself uses joules. There is no finite ice mass/exhaustion model.

`step(state, parameters, dt_s, atmospheric_w_m2, water_w_m2)` holds area and
parameters fixed. Both ports are signed positive inward; atmospheric heat enters
debris, water heat enters the ice reservoir. Backward Euler matches the existing
fixed-area closure for its supported nonnegative water flux. Conductive and
water energy repay cold debt before generating latent energy. Cooling increases
debt and does not refreeze departed meltwater. This is a lumped zero-C reference
interface, not a cold-ice temperature PDE.

`remap_area` requires explicit `AreaContent` for added or removed material,
including its debris energy, cold debt, and material-associated latent history
per square metre. Zero content must be declared explicitly. Existing totals
are never multiplied by a new area. Returned transfers are positive inward;
their signed energy is debris transfer minus debt transfer plus latent transfer.
Contraction exports declared material, including latent history; if whole-system
lifetime melt is needed, accumulate `step.latent_j` outside this moving domain.
Remap and time-step ledgers are separate. No undeclared geometry work or mixing
source is modeled. Zero area, nonfinite arithmetic, invalid debt, and unsupported
transfer combinations raise errors. Floating point residuals should be checked
using a scale-aware tolerance by callers.

`step_bare(state, dt_s, atmospheric_w_m2, water_w_m2)` sends both signed
ports directly into ice and requires zero debris energy. It takes no fictitious
debris capacity. An externally computed conductive exchange may be included
in either supplied signed port, with that partition tracked by the caller.

Run: `python3 -m unittest research.environment_coupled_2026_11_01.energy_ports.test_energy`
