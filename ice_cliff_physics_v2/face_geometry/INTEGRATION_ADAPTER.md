# Coupled latent-energy integration adapter

`integration_adapter.py` connects the parent model's already-conserved net latent
energy partition to the height-band geometry. It does not introduce an additional
water-heat budget, split stored debris heat twice, or use the routing proxy as the
physical ice volume.

```python
profile = initialize_profile(step_geometry_dict, n_bands=32)
profile, ledger, proxy = advance_with_partition(
    profile,
    {'bare_air': E_air, 'bare_water': E_water,
     'covered_dry': E_covered_dry, 'covered_water': E_covered_water},
    floor_latent_J=E_floor,
    water_stage=stage,
    cover_z_intervals=[(covered_lo, covered_hi)],
)
```

The exact seven initial geometry keys match SteppedCliff. The four class-energy
keys are required, even for zero-valued inactive classes. Positive energy for a
zero-area old-profile class fails; it is not reassigned to a convenient surface.
The function captures the old geometry's four exclusive class-area partitions,
allocates each class's energy proportionally to its per-band area, and then calls
the existing disjoint face-first/floor-second operator. Allocation uses a
nonnegative remainder scheme and reports its residual.

A partial-height class inside a band produces a band-average displacement. This
is an explicit height-discretization approximation; the adapter does not silently
create growing sets of adaptive bands. The numerical test refines that spatial
partition, independently of total-energy conservation.

The returned proxy has exactly the seven SteppedCliff fields and uses the lowest
band's front as face_x_m. It is a necessary reduced routing representation, not a
resolved water profile and **not the actual remaining ice volume**. Use
`profile.ice_volume_m3` or `ledger['after_profile_ice_volume_m3']` for physical
volume. The ledger also reports the proxy volume and signed difference so misuse
is visible.

The parent must require temperate ice, no cold-debt/refreezing path, nonnegative
latent class energies, and an explicit adiabatic-horizontal-ledge policy for this
bounded integration mode. The helper does not infer this policy or claim full
interface physics. Height masks are Eulerian, not transported debris; stage is a
supplied hydraulic approximation. New face area, material enthalpy and reservoir
remapping remain parent responsibilities.

Seven added tests verify all four classes plus floor, the planar uniform-flux
limit, spatial partition refinement, the difference between true and proxy
volume, snapshot timing before floor-area growth, and invalid energy refusals.
Together with the earlier fifteen geometry tests, all twenty-two profile tests
pass. This is coupled implementation verification, not environmental accuracy or
historical originality.
