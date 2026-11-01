# Reproduced bed-update defects

Pinned source f3602e35e1595bcdbee4208a354a34944f1ae974, swe2d.py line175.
The global maximum(zb-dz,pool_floor) changes non-ice tail geometry even at
zero elapsed time, and returned melt uses requested rather than actual lowering.
At dx=.025, the zero-time call raises some non-ice entries by1.07375m;
interior non-ice bed-volume change is0.79m3. Exact evidence is saved in
source_reproduction.json. Uniform100W/m2 heat fixtures are synthetic and do
not quantify the completed campaign error.

The paired3.5-day dx=.1 development diagnostic compares the pinned operator
with ice-only unbounded lowering. Loss times are3.164769045 and3.188464521d.
Counted-minus-geometric melt is0.0166157975 and0m3. These are diagnostic
results, not revised medium/fine production scores or proof about day12 recovery.
Full records and a reproducible figure are included.

operator.py supplies a pure ice-only lowering operator with explicit refusal
at an optional domain floor. It is not automatically installed, and is not a
complete liquid mass/momentum or multirate energy correction. Removing the
old floor requires boundary/domain revalidation; a physical bedrock interface
requires its own heat treatment instead of clipping latent melt silently.

Do not run the expensive finest case on the defective pinned operator.
Repair, verify and refreeze first. Recompute the grid sequence for the revised
model; the old GCI cannot certify it. See ../numerics and CLAUDE_HANDOFF.md.
