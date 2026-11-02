# Follow-up: sharp-class alignment resolves the audited basal alias

This is a versioned reanalysis of the same frozen synthetic6h input used in
geometry_review.md, after the parent enabled lossless insertion of the current
water-stage and cover interfaces before thermal allocation. Old outputs and the
reviewed source snapshot remain preserved; this is not an independent new case.

| Initial bands | Actual final/max bands | Basal front m | Floor z m |
|---|---:|---:|---:|
|4|78|0.0070294452171853665|-0.08704064977938475|
|16|90|0.0070294452171853665|-0.08704064977938475|
|64|138|0.0070294452171853665|-0.08704064977938475|
|128|202|0.0070294452171853665|-0.08704064977938475|
|256|330|0.0070294452171853670|-0.08704064977938475|

Basal-front spread is8.67e-19m across initial grids. This resolves findingG2 for
this sharp-class alias: coarse initial bands now preserve the exposed2cm basal
patch instead of averaging its recession with other classes. The actual adaptive
band count—not just the initial request—is reported. This does not remove the
need to refine timestep, initial curved geometry or within-class thermal forcing.

The parent also added surface_partition_sample_time_s and state_sample_time_s,
explicitly identifying old-step surface masks and end-step geometry/reservoirs.
That resolves the specific output-time ambiguity in findingG1 without replacing
old energy-domain fields by new geometry values.

Executable reproduction: reproduce_aligned_geometry_audit.py with the package
path. Results and source/input hashes are in
geometry_height_alignment_confirmation.json. The original source audit and
unaligned refinement files remain distinct. No field accuracy or originality
claim follows from this numerical correction.
