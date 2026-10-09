#!/usr/bin/env python
"""Separately named k-prior sensitivity of the held-out experiment (`heldout_v2_kfree`); the frozen
`records/heldout_v2` result is not modified (docs/data_use_register.md, deviation D3).

The k prior (HCRF/albedo) of heldout_v2 is built from S6 2017 ARF and overlaps the reverse-fold test plots.
Here it is replaced IN MEMORY (priors._defaults patched; no source edit) by:
  independent : N(1.0, 0.25)  - instrument assumption already used for Stibal 2017 (no S6 information)
  range       : N(0.8|0.9|1.0, 0.175) and N(0.9, 0.0875 | 0.35)
Physics tables are rebuilt in a NEW directory; the frozen emulators are linked read-only (content-tagged).

    python common/run_budgeted.py --background --mem-mb 2000 --name heldout_kfree -- \
        python3 phase4/heldout_kfree.py --variant independent
"""
import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "phase2"))

VARIANTS = {"independent": (1.0, 0.25), "mean0.8": (0.8, 0.175), "mean0.9": (0.9, 0.175), "mean1.0": (1.0, 0.175),
            "sd_half": (0.9, 0.0875), "sd_double": (0.9, 0.35)}
OPTICS = ("tierA_empirical", "tddft_D", "measured_mac_C")


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--variant", choices=list(VARIANTS), required=True)
    a = p.parse_args(argv)
    mk, sk = VARIANTS[a.variant]
    import priors
    orig = priors._defaults

    def patched(rho_bottom=None):
        d = orig(rho_bottom)
        d.update(mu_k=mk, sd_k=sk)
        return d
    priors._defaults = patched
    out = os.path.join(HERE, "results", "heldout_v2_kfree", a.variant)
    cache = os.path.join(out, "cache")
    os.makedirs(cache, exist_ok=True)
    src = os.path.join(HERE, "results", "heldout_v2", "cache")
    for f in os.listdir(src):
        if f.startswith("heldout_") and f.endswith(".npz") and any(f"heldout_{o}_sza" in f for o in OPTICS):
            dst = os.path.join(cache, f)
            if not os.path.exists(dst):
                os.symlink(os.path.join(src, f), dst)
    import heldout
    heldout.main(["--outdir", out, "--optics", *OPTICS, "--workers", "1"])
    open(os.path.join(out, "K_PRIOR.txt"), "w").write(f"k ~ N({mk}, {sk}) [{a.variant}]\n")


if __name__ == "__main__":
    main()
