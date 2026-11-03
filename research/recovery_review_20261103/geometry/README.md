# Fixed physical divider/pool geometry for bare v5

`FixedGeometrySim` constructs a new supported bare-ice v5 state with a fixed0.25m lane on each side of a fixed0.25m divider. The pool is0.75m wide and the lengths0.5/0.5/0.25m are constant. The divider occupies1/2/4cells atdx0.25/0.125/0.0625m instead of shrinking physically withdx. The hydraulic kernel is the exact vendored pinned kernel, unchanged. This is a new geometry version, not a correction to old saved results.

Defaults preserve the declaredQ0.07m³/s,Tw0.5K,h200Wm⁻²K⁻¹. No larger heating is silently introduced. Defaults are ideal inputs, not measured field forcings. The default minimum is4cells per lane;1/2-cell geometry checks require an explicit`min_lane_cells=1` override. Four cells is only a minimum applicability gate, not evidence of convergence. The next supported refinement family isdx0.0625/0.03125/0.015625m (4/8/16lane cells), with identical physical geometry and inputs.

Construction validates alignment, then replaces the temporary parent initial geometry/state BEFORE any evolution. The parent phase solver's documented initial-condition reset updates the reference and records`initial_condition_definition`; no evolved mass, heat or time history is erased. A post-evolution reset retains its explicit correction and removes the physical-only certificate. Native physical-time phase and thermal guards are inherited. Deposits are refused by this builder.

The builder reconstructs all padded arrays, masks, routes, inlet span, tail and ghost bed values. `byp_cols` is overridden for the wider divider. Physical lane/reach/wall areas and total inlet width remain identical across meshes. The actual inlet boundary's momentum integral isQ for all three meshes. Wall cells are nonice and cannot melt.

Eight tests cover dimensions, areas, fixed input width, actual boundary values, route slices, ghost exclusion, initial and post-evolution ledgers, source step, no wall melt and alignment/refusal. A single0.01s physical step is checked only as an integration smoke test. No full20s campaign, event-timing convergence, contact accuracy, recovery, independent validation or field superiority is established.

Usage:

```python
from fixed_geometry import FixedGeometrySim
s=FixedGeometrySim(dx=0.0625)
result=s.step_physical(0.01)
```

Run`python3 -m unittest -v test_fixed_geometry`. NumPy and Numba are required. Vendored files are reproducibility dependencies from the published repair ZIP, not modified production code. Claude's worktree, jobs, original solver and frozen outputs are untouched.
