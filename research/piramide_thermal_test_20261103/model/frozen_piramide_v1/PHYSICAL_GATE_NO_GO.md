# V1 physical gate: no-go before camera access

The initial protocol and sources were published at `1b6830fb5edd90238bb0332e4287854c8627986d`. All twelve normalized thermal method/scenario predictions were published at `5ea3781815eb04b9e10b598b03d712d61a2e21a1`, before any training or evaluation camera outcomes were read.

All six forward-inference cases produce negative basal heat on some intervals. The frozen melt-only conversion explicitly refuses this. The predeclared primary therefore cannot proceed to melt calibration or scoring under V1. No candidate or geometry is dropped to manufacture a valid family, and the derivative comparator is not substituted as a new primary.

For the middle-gap/left primary, the inferred diffusivity is 2.3654178463540555e-6 m²/s. Fifty native half-hour intervals have negative integrated basal heat. Their total negative magnitude is about 0.588874 K m, approximately 0.523% of the net 28-day heat. Small aggregate magnitude does not change the failed physical scope condition.

Independent refinement preserves the negative flux. The evaluation/train heat ratio shifts from 0.72713922 to about 0.72708061, −0.00806%, while the derivative comparator shifts about +0.000779%. Numerical heat-ratio stability is not a repair for missing cold-content physics.

The physically justified next version adds cold-ice energy storage and a complementary melting boundary. Negative heat cools ice; later heat first warms it before melting. It does not interpret negative latent-equivalent energy as refreezing without a water supply. This is established phase-boundary physics and a new implementation milestone, not a historical novelty claim. Its code/specification must be frozen separately before outcome access. Both training and evaluation camera outcomes remain sealed at this decision.
