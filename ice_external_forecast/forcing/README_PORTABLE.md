# Rochemelon forcing packet

Data and original code only; measured bed-change targets are absent. Read FORCING_INPUT_GATE.md before any model evaluation.

FORCING_GRAPH_SAMPLES.csv is factual graph-digitized forcing, with per-row source pixel coordinates. FORCING_PROVENANCE.json records source URL, hash, anchors and conservative coordinate error. No raw sensor CSV was located. Do not describe samples of a plotted line as55 independent observations.

The publisher PDF, source JPEG, page render and selection-proof image stay in the private working directory. Their reuse license was not verified, so they are deliberately omitted from this publication-ready packet. Source links remain in the manifest. Own plots of the extracted numerical values are included separately.

Reproduce sources into a private directory with prepare_sources.py (requiresPyMuPDF1.24.14). Then run digitize_forcing.py with environment ROCHE_SOURCE_DIR pointing there. This requiresNumPy andPillow. The extraction script reads only forcing-panel pixel rows; the withheld middle bed panel is not sampled. Graphical proof is generated locally for quality control and should not be republished without source rights review.

The script output is a2-hour sampling of the published curves, with a preferred12–108h regular-breach window. Any finer model forcing should be a labelled interpolation, with no extrapolation. Source absolute-clock timezone, control section and downstream geometric matching remain unresolved.
