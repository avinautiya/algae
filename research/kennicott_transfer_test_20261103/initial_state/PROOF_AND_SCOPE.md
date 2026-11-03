# A computable bound on initial thermal-state memory

For a fixed passive thermal chain with node capacities c_i>0 and edge conductances g_e>0, compare two solutions with identical Dirichlet boundary histories. Their temperature difference u has homogeneous boundary conditions. Let p solve the discrete harmonic problem p_left=0,p_right=1. In a chain, p_i equals cumulative resistance to nodei divided by total resistance. Discrete summation by parts gives

    d/dt sum_i c_i*p_i*u_i = −q_right_difference.

Hence cumulative right-boundary heat difference equals sum_i c_i*p_i*(u_i(0)−u_i(T)). Positivity of the passive semigroup implies the impulse-to-right absorption fraction by any finite time is between0andp_i. For |u_i(0)|<=epsilon_i,

    |DeltaQ_right(0,T)| <= sum_i c_i*p_i*epsilon_i

for every finiteT. Divide by latent heat only on an independently justified heat-to-melt branch. This bound is tighter than assigning all initial heat to the ice, since some exits the upper boundary. The code's independent augmented matrix-exponential test integrates physical right flux and checks the identity directly.

This is established adjoint/absorption-probability mathematics. Its proposed role is to quantify a real installation/initial-state uncertainty in the glacier comparison rather than silently deleting initial days. It does not certify the chosen temperature-error radii, missing boundary forcing, unknown field capacity, time-varying conductance, nonlinear phase, moving geometry or the actual observed target. Time-varyingp introduces an extra derivative term, so applying this fixed-graph bound to an adaptive graph without further derivation is invalid.
