# Prescribed-flux moving-profile benchmark

This exact geometric solution is provided as a verification reference, not a new
mathematical theorem or a solution of the conjugate thermal/flow problem.

Let width W, initial front x0, upstream edge U, initial floor z0 and crest Z be fixed.
A supplied stage zw lies strictly between z0 and Z. The ice is initially a step.
Prescribe constant **net latent** vertical-face flux qw below zw,qa above zw and
horizontal-floor flux qb. New horizontal ledges have prescribed zero normal flux.
Set `aw=qw/(rho*L)`, `aa=qa/(rho*L)`, `b=qb/(rho*L)`, `Hw=zw-z0`, `Ha=Z-zw`.

The front kinematics are

    x_w(t)=x0+aw*t,   x_a(t)=x0+aa*t,   z_floor(t)=z0-b*t.

As long as the front and floor remain inside the finite domain, direct area
integration gives total geometric removal

    DeltaV/W = (aw*Hw + aa*Ha + b*(x0-U))*t + aw*b*t².

The cross term is not an extra source. The face newly exposed by floor lowering
adds `aw*b*t²/2` to face melt, and growth of the basal footprint adds the other
half to floor melt. Integrating the prescribed boundary fluxes gives exactly the
same total. A per-step face-first/floor-second splitting converges to this
continuum geometry at first order in coupling timestep when beginning-step
floor-flux area is used. Spatial bands must resolve the jump at zw; a separate
linear-flux benchmark tests first-order approximation in height.

For comparison, the reduced model enforces a single horizontal front at every
height while receiving the same height-integrated vertical flux. Its exact ODE is

    dx/dt = aw + (aa-aw)*Ha/(Hw+Ha+b*t).

For b>0 the solution is

    x(t)=x0+aw*t+(aa-aw)*Ha/b * log(1+b*t/(Hw+Ha)).

Its volume removal is

    DeltaV_uniform/W = (Hw+Ha)*(x(t)-x0) + b*t*(x(t)-U).

This reduced model is mathematically consistent but cannot produce the
nonuniform front or its undercut. Checking both exact limits prevents mistaking
a representational limitation for a failure of energy bookkeeping.

The implementation also checks actual polygon area independently. At nonuniform
bands it reports ceiling/terrace area so a later thermal solver cannot
accidentally omit these boundaries while claiming complete interface physics.
