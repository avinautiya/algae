# Conditional recovery exclusion bound

This is a useful elementary comparison argument, **not claimed new mathematics**.
It applies to an explicit uniform-toe reduced model. It is not a certificate for
an arbitrary spatially evolving cliff or an unobserved englacial pathway.

Write `A(eta) eta_dot = I - Q_b - Q_t`, with `A(eta)>0`, `Q_t>=0`,
`I<=Imax`, and a bypass rating bounded throughout the horizon by
`Q_b >= kmin max(eta-bmax,0)^(3/2)`. All sources of reservoir water, including
any meltwater returning to it, must be in `Imax`. Then

`eta <= etaU = max(eta(0), bmax+(Imax/kmin)^(2/3))`.

Proof: at the upper boundary the net volume rate is nonpositive. The scalar
barrier cannot be crossed. This does not require a constant storage area.

Let the uniformly lowering toe floor have area `Ai`, initial elevation `zi0`,
and a contact-path elevation at least `zi+Dmin`, with `Dmin>=0`. If the total
positive heat accessible to that floor is at most `Eplus`, its spatially uniform
cold reserve is at least `Cc`, and no other mechanism removes ice, then

`zi >= zi0 - max(Eplus-Cc,0)/(rho_i L Ai)`.

Thus contact is excluded if this floor bound plus `Dmin` exceeds `etaU`.
Otherwise the result is **undetermined**: available heat does not establish
delivery to ice, and a stage upper bound need never actually be attained.
The code deliberately does not convert this necessary condition into a forecast
of recovery or a statistical probability.

The heat bound must include advected heat, atmospheric floor heat, debris heat
release and any mechanical dissipation assigned to that floor. Heat on a
vertical face cannot be counted as floor lowering. Nonuniform incision can
open a narrow path with much less energy than uniform lowering; in that case
this uniform-toe bound is inapplicable. Collapse or an under-floor conduit also
invalidates the proof. A measured multi-epoch cliff profile does not by itself
constrain the required heat or discharge bounds.
