#!/usr/bin/env python
"""
Surrogate (look-up table) interface for coupling the pigment-aware optics to other models (e.g. an SEB or
land-ice model): broadband albedo, algal albedo reduction and instantaneous forcing as functions of

    log10 abundance (cells mL^-1), community fraction f_n, bubble radius r (um), dust (ppb), solar zenith

for a frozen optical model ('tddft_D', 'tddft_C', 'measured_mac_C' or 'tierA_empirical', as in heldout.py).
Interpolation: multilinear within each SZA node's emulator, linear in SZA between nodes.

Every call returns domain flags; values are NOT extrapolated silently:
    in_domain      all inputs inside the tabulated ranges
    flags          which inputs were outside (and clipped to the boundary)
Accuracy against direct BioSNICAR runs is measured by `benchmark()` (off-grid states, off-node SZA) and
recorded in records/surrogate_benchmark.json; the API refuses to load a table without that record unless
allow_unbenchmarked=True.

    python phase4/surrogate_api.py --optics tddft_D --sza 40 50 60 --outdir phase4/results/surrogate
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "phase2"))

OUTPUTS = ("bba", "rf_algae")


class Surrogate:
    def __init__(self, emulators: dict, optics: str, benchmark: dict | None = None):
        self.optics = optics
        self.sza = np.array(sorted(emulators), float)
        self.ems = [emulators[z] for z in sorted(emulators)]
        self.benchmark = benchmark
        e0 = self.ems[0]
        self.names = e0.names
        self.ranges = {n: (float(e0.axes[n][0]), float(e0.axes[n][-1])) for n in e0.names}
        self.ranges["sza"] = (float(self.sza[0]), float(self.sza[-1]))
        self._I = [{k: em.interpolator(k) for k in OUTPUTS} for em in self.ems]
        clean = {}
        for k, em in enumerate(self.ems):           # clean-ice reference: lowest abundance node
            clean[k] = em
        self._clean = clean

    def _coords(self, em, state):
        return np.array([[em.coord(n)[0] if len(em.axes[n]) == 1 else
                          (np.log10(state[n] + 100.0) if n == "dust_ppb" else state[n]) for n in em.active]])

    def __call__(self, log_b, sza, r_um=3000.0, f_n=0.5, dust_ppb=0.0):
        state = dict(log_b=float(log_b), r_um=float(r_um), f_n=float(f_n), dust_ppb=float(dust_ppb), sza=float(sza))
        flags = []
        for n, (lo, hi) in self.ranges.items():
            if n in state and not lo - 1e-9 <= state[n] <= hi + 1e-9:
                flags.append(f"{n}={state[n]:g} outside [{lo:g}, {hi:g}] (clipped)")
                state[n] = float(np.clip(state[n], lo, hi))
        j = int(np.clip(np.searchsorted(self.sza, state["sza"]) - 1, 0, len(self.sza) - 2)) if len(self.sza) > 1 else 0
        if len(self.sza) == 1:
            ws = [(0, 1.0)]
        else:
            t = (state["sza"] - self.sza[j]) / (self.sza[j + 1] - self.sza[j])
            ws = [(j, 1 - t), (j + 1, t)]
        out = {k: 0.0 for k in OUTPUTS}
        clean_bba = 0.0
        for k, w in ws:
            em = self.ems[k]
            x = self._coords(em, state)
            for o in OUTPUTS:
                out[o] += w * float(self._I[k][o](x)[0])
            xc = self._coords(em, dict(state, log_b=self.ranges["log_b"][0]))
            clean_bba += w * float(self._I[k]["bba"](xc)[0])
        out["dalpha_algae_vs_lowest_abundance"] = clean_bba - out["bba"]
        out.update(in_domain=not flags, flags=flags, optics=self.optics,
                   benchmarked=self.benchmark is not None)
        return out


def build(optics, sza_nodes, phase1_l2=None, biosnicar=None, workers=2, cache_dir=".", rho_bottom=690.0):
    import emulator as E
    from heldout import OPTICS
    ems = {}
    for z in sza_nodes:
        cfg = E.EmulatorConfig(sza=int(z), spacecraft="S2A", rho_bottom=rho_bottom, dust_ppb=E.DUST_NODES_PPB,
                               photosynthetic=True, **OPTICS[optics])
        ems[float(z)] = E.build_emulator(cfg, phase1_l2=phase1_l2, biosnicar=biosnicar, workers=workers,
                                         cache=os.path.join(cache_dir, f"surrogate_{optics}_sza{int(z)}.npz"))
    return ems


def benchmark(sur: Surrogate, optics, n=40, seed=0, phase1_l2=None, biosnicar=None, rho_bottom=690.0):
    """Surrogate vs direct BioSNICAR at random off-grid states and an off-node SZA (mid-way between
    nodes). Direct runs use the SZA by itself (BioSNICAR illumination at the integer SZA)."""
    import emulator as E
    from heldout import OPTICS
    rng = np.random.default_rng(seed)
    z = float(np.round(0.5 * (sur.sza[0] + sur.sza[1]))) if len(sur.sza) > 1 else float(sur.sza[0])
    cfg = E.EmulatorConfig(sza=int(z), spacecraft="S2A", rho_bottom=rho_bottom, dust_ppb=E.DUST_NODES_PPB,
                           photosynthetic=True, **OPTICS[optics])
    lo = sur.ranges
    st = [dict(log_b=rng.uniform(lo["log_b"][0] + 0.5, lo["log_b"][1] - 0.5),
               f_n=rng.uniform(0, 1) if lo["f_n"][1] > lo["f_n"][0] else lo["f_n"][0],
               r_um=float(np.exp(rng.uniform(np.log(800), np.log(12000)))),
               dust_ppb=float(np.exp(rng.uniform(np.log(5e4), np.log(1e6))))) for _ in range(n)]
    d = E.direct_forward(cfg, st, phase1_l2=phase1_l2, biosnicar=biosnicar)
    import biosnicar_bridge as bb
    runner = bb.BioSNICARRunner(bb.locate_biosnicar(biosnicar))
    s = np.array([sur(x["log_b"], z, runner.snap_radius(x["r_um"], cfg.ice_mode), x["f_n"], x["dust_ppb"])["bba"]
                  for x in st])
    err = s - d[:, 4]
    return dict(optics=optics, sza_test=z, sza_nodes=sur.sza.tolist(), n=n, bba_mae=float(np.mean(np.abs(err))),
                bba_max_abs=float(np.max(np.abs(err))), bba_bias=float(err.mean()))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--optics", default="tddft_D")
    p.add_argument("--sza", type=float, nargs="+", default=[40.0, 50.0, 60.0])
    p.add_argument("--phase1-l2", default=os.path.join(HERE, "..", "phase1", "results", "level2"))
    p.add_argument("--biosnicar", default=None)
    p.add_argument("--workers", type=int, default=2)
    p.add_argument("--outdir", default=os.path.join(HERE, "results", "surrogate"))
    a = p.parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)
    ems = build(a.optics, a.sza, a.phase1_l2, a.biosnicar, a.workers, a.outdir)
    sur = Surrogate(ems, a.optics)
    bm = benchmark(sur, a.optics, phase1_l2=a.phase1_l2, biosnicar=a.biosnicar)
    print(json.dumps(bm, indent=1))
    import provenance as PV
    PV.atomic_write_text(os.path.join(a.outdir, f"benchmark_{a.optics}.json"), json.dumps(bm, indent=1))


if __name__ == "__main__":
    main()
