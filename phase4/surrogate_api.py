#!/usr/bin/env python
"""
Surrogate (look-up table) interface for coupling the pigment-aware optics to other models (e.g. an SEB or
land-ice model): broadband albedo, algal albedo reduction and instantaneous algal forcing as functions of

    log10 abundance (cells mL^-1), community fraction f_n, bubble radius r (um), dust (ppb), solar zenith

for a frozen optical model ('tddft_D', 'tddft_C', 'measured_mac_C' or 'tierA_empirical', as in heldout.py).
Interpolation: multilinear within each SZA node's emulator, linear in SZA between nodes.

Zero-algae reference: the emulator's rf_algae is the forcing relative to a run WITHOUT algae (same ice,
dust and SZA), so the algal albedo reduction is d_alpha = rf_algae / SW_model, with SW_model the irradiance
the emulator used - a genuine algae-free reference, not the lowest-abundance node.

Qualification: a surrogate is usable only with a benchmark record (benchmark()) that
  * was computed for exactly these emulators (their content tags),
  * samples the whole supported domain (all axes incl. their ends, SZA at and between nodes),
  * meets the tolerances BBA max |error| <= TOL_BBA and rf_algae max |error| <= TOL_RF_ABS + TOL_RF_REL x |rf|.
  * is bound to the surrogate's own code, tolerances and benchmark design (`benchmark_fingerprint()`); a
    record made by other code or tolerances is STALE and rejected (the emulator tags already bind the
    configuration, Phase 1 content, measured data, physics code and BioSNICAR revision).
Surrogate(...) raises without such a record unless allow_unqualified=True, and every call reports
`qualified`. Requests outside the tabulated ranges RAISE `DomainError` by default; strict=False clips AND
flags them (`in_domain`, `flags`) for exploratory use only.

Contract (units and conventions):
  * abundance: log10 cells per mL of MELTWATER, 1.0-6.0 (10 to 1e6); zero algae is NOT a domain point:
    the algae-free reference enters only through dalpha/rf (relative to a run without algae);
  * f_n: fraction A. nordenskioeldii (0-1); r_um: bubble optical radius (um) of bubbly ice, crust
    450 kg m-3 over 690 kg m-3, 2 cm crust; dust_ppb: mineral dust (ng g-1) uniform in the crust;
  * illumination: BioSNICAR clear-sky sub-Arctic summer spectrum (incoming=3), direct beam, SZA between
    the emulator nodes (linear in SZA; no extrapolation);
  * spectral support: 300-2500 nm; bba = irradiance-weighted broadband albedo on that band;
    rf_algae = instantaneous algal forcing (W m-2) with the clear-sky model SW at that SZA, SIGNED;
    dalpha_algae = rf_algae / SW_model (signed, > 0 darkening).
Tolerances (set from the error budget BEFORE any benchmark result: half the protocol's minimum meaningful
albedo change 0.01): |bba error| <= 0.005, |dalpha error| <= 0.005, |rf error| <= 0.5 W m-2 + 5 %.

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
TOL_BBA = 0.005
TOL_RF_ABS, TOL_RF_REL = 0.5, 0.05          # W m^-2, fraction
TOL_DALPHA = 0.005
DESIGN_VERSION = "lhs+corners+mid-sza/2"


class QualificationError(RuntimeError):
    pass


class DomainError(ValueError):
    pass


def benchmark_fingerprint() -> str:
    """Binds a benchmark record to this module's code, the tolerances and the design version."""
    import provenance as PV
    return PV.fingerprint(dict(code=PV.file_sha256(os.path.abspath(__file__)),
                               tol=[TOL_BBA, TOL_DALPHA, TOL_RF_ABS, TOL_RF_REL], design=DESIGN_VERSION))


def _tags(ems):
    return {str(float(z)): str(em.meta.get("tag", "")) for z, em in sorted(ems.items())}


class Surrogate:
    def __init__(self, emulators: dict, optics: str, benchmark: dict | None = None, allow_unqualified=False):
        self.optics = optics
        self.sza = np.array(sorted(emulators), float)
        self.ems = [emulators[z] for z in sorted(emulators)]
        self.benchmark = benchmark
        self.qualified = qualifies(benchmark, emulators)
        if not self.qualified and not allow_unqualified:
            raise QualificationError("no passing benchmark record for these emulators (run benchmark()); "
                                     "pass allow_unqualified=True only for benchmarking itself")
        e0 = self.ems[0]
        self.names = e0.names
        self.ranges = {n: (float(e0.axes[n][0]), float(e0.axes[n][-1])) for n in e0.names}
        self.ranges["sza"] = (float(self.sza[0]), float(self.sza[-1]))
        self._I = [{k: em.interpolator(k) for k in OUTPUTS} for em in self.ems]
        self._sw = [float(em.meta.get("sw_down", np.nan)) for em in self.ems]

    def _coords(self, em, state):
        return np.array([[(np.log10(state[n] + 100.0) if n == "dust_ppb" else state[n]) for n in em.active]])

    def __call__(self, log_b, sza, r_um=3000.0, f_n=0.5, dust_ppb=None, strict=True):
        state = dict(log_b=float(log_b), r_um=float(r_um), f_n=float(f_n), sza=float(sza),
                     dust_ppb=float(self.ranges["dust_ppb"][0] if dust_ppb is None and "dust_ppb" in self.ranges
                                    else (dust_ppb or 0.0)))
        flags = []
        for n, (lo, hi) in self.ranges.items():
            if n in state and not lo - 1e-9 <= state[n] <= hi + 1e-9:
                flags.append(f"{n}={state[n]:g} outside [{lo:g}, {hi:g}] (clipped)")
                state[n] = float(np.clip(state[n], lo, hi))
        if flags and strict:
            raise DomainError("; ".join(f.replace(" (clipped)", "") for f in flags) + " - out of the supported domain")
        if len(self.sza) == 1:
            ws = [(0, 1.0)]
        else:
            j = int(np.clip(np.searchsorted(self.sza, state["sza"]) - 1, 0, len(self.sza) - 2))
            t = (state["sza"] - self.sza[j]) / (self.sza[j + 1] - self.sza[j])
            ws = [(j, 1 - t), (j + 1, t)]
        out = {k: 0.0 for k in OUTPUTS}
        dalpha = 0.0
        for k, w in ws:
            x = self._coords(self.ems[k], state)
            vals = {o: float(self._I[k][o](x)[0]) for o in OUTPUTS}
            for o in OUTPUTS:
                out[o] += w * vals[o]
            dalpha += w * vals["rf_algae"] / self._sw[k]          # relative to the algae-free run
        out["dalpha_algae"] = dalpha
        out.update(in_domain=not flags, flags=flags, optics=self.optics, qualified=self.qualified)
        return out


def qualifies(benchmark, emulators) -> bool:
    if not benchmark:
        return False
    return (bool(benchmark.get("passed")) and benchmark.get("emulator_tags") == _tags(emulators)
            and benchmark.get("benchmark_fingerprint") == benchmark_fingerprint())


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


def domain_design(ranges, sza_nodes, n_per_sza=24, seed=0, radii=None):
    """States covering the supported domain: a Latin hypercube over (log B, f_n, ln r, ln dust) plus all
    corner values of log B and dust, at every SZA node and mid-way between neighbouring nodes."""
    from scipy.stats import qmc
    rng = np.random.default_rng(seed)
    zs = list(sza_nodes) + [0.5 * (a + b) for a, b in zip(sza_nodes[:-1], sza_nodes[1:])]
    U = qmc.LatinHypercube(d=4, seed=seed).random(n_per_sza)
    lo = {k: v[0] for k, v in ranges.items()}
    hi = {k: v[1] for k, v in ranges.items()}
    st = []
    for z in zs:
        rows = [dict(log_b=lo["log_b"] + u[0] * (hi["log_b"] - lo["log_b"]),
                     f_n=lo["f_n"] + u[1] * (hi["f_n"] - lo["f_n"]),
                     r_um=float(np.exp(np.log(lo["r_um"]) + u[2] * (np.log(hi["r_um"]) - np.log(lo["r_um"])))),
                     dust_ppb=float(np.exp(np.log(lo["dust_ppb"]) + u[3] * (np.log(hi["dust_ppb"]) - np.log(lo["dust_ppb"])))))
                for u in U]
        for lb in (lo["log_b"], hi["log_b"]):
            for d in (lo["dust_ppb"], hi["dust_ppb"]):
                rows.append(dict(log_b=lb, f_n=float(rng.uniform(lo["f_n"], hi["f_n"])), dust_ppb=d,
                                 r_um=float(np.exp(rng.uniform(np.log(lo["r_um"]), np.log(hi["r_um"]))))))
        for r in rows:
            r["sza"] = float(round(z))
            if radii is not None:                       # direct BioSNICAR runs use its look-up-table radii
                r["r_um"] = float(radii[np.argmin(np.abs(np.log(radii) - np.log(r["r_um"])))])
        st += rows
    return st


def benchmark(emulators, optics, n_per_sza=24, seed=0, phase1_l2=None, biosnicar=None, rho_bottom=690.0):
    """Surrogate vs direct BioSNICAR (BBA and algal forcing) over domain_design(); returns the record
    used for qualification."""
    import emulator as E
    import biosnicar_bridge as bb
    from heldout import OPTICS
    sur = Surrogate(emulators, optics, allow_unqualified=True)
    runner = bb.BioSNICARRunner(bb.locate_biosnicar(biosnicar))
    e0 = next(iter(emulators.values()))
    radii = np.asarray(e0.axes["r_um"], float)
    states = domain_design(sur.ranges, list(sur.sza), n_per_sza, seed, radii=radii)
    rows = []
    for z in sorted({s["sza"] for s in states}):
        cfg = E.EmulatorConfig(sza=int(z), spacecraft="S2A", rho_bottom=rho_bottom, dust_ppb=E.DUST_NODES_PPB,
                               photosynthetic=True, **OPTICS[optics])
        sub = [s for s in states if s["sza"] == z]
        d = E.direct_forward(cfg, sub, phase1_l2=phase1_l2, biosnicar=biosnicar)
        for s, dv in zip(sub, d):
            p = sur(s["log_b"], z, s["r_um"], s["f_n"], s["dust_ppb"])
            sw = bb.sw_down_clear_sky(z, None, cfg.day_of_year)
            rows.append(dict(**s, at_node=bool(np.any(np.isclose(z, sur.sza))), bba_direct=float(dv[4]),
                             bba_sur=p["bba"], rf_direct=float(dv[5]), rf_sur=p["rf_algae"],
                             dalpha_direct=float(dv[5]) / sw, dalpha_sur=p["dalpha_algae"]))
    import pandas as pd
    t = pd.DataFrame(rows)
    eb, er, ed = t.bba_sur - t.bba_direct, t.rf_sur - t.rf_direct, t.dalpha_sur - t.dalpha_direct
    ok_rf = np.abs(er) <= TOL_RF_ABS + TOL_RF_REL * np.abs(t.rf_direct)
    rec = dict(optics=optics, sza_nodes=sur.sza.tolist(), n=int(len(t)), emulator_tags=_tags(emulators),
               tolerances=dict(bba_max_abs=TOL_BBA, dalpha_max_abs=TOL_DALPHA, rf_abs=TOL_RF_ABS, rf_rel=TOL_RF_REL),
               bba_max_abs=float(np.abs(eb).max()), bba_mae=float(np.abs(eb).mean()),
               bba_max_abs_off_node=float(np.abs(eb[~t.at_node]).max()) if (~t.at_node).any() else None,
               rf_max_abs=float(np.abs(er).max()), rf_mae=float(np.abs(er).mean()),
               rf_fraction_within_tol=float(ok_rf.mean()),
               dalpha_max_abs=float(np.abs(ed).max()),
               dalpha_max_abs_lowest_abundance=float(np.abs(ed[t.log_b == t.log_b.min()]).max()),
               bba_max_abs_by_sza={str(z): float(np.abs(eb[t.sza == z]).max()) for z in sorted(t.sza.unique())},
               benchmark_fingerprint=benchmark_fingerprint(), design_version=DESIGN_VERSION,
               passed=bool(np.abs(eb).max() <= TOL_BBA and ok_rf.all() and np.abs(ed).max() <= TOL_DALPHA))
    return rec, t


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--optics", default="tddft_D")
    p.add_argument("--sza", type=float, nargs="+", default=[40.0, 50.0, 60.0])
    p.add_argument("--phase1-l2", default=os.path.join(HERE, "..", "phase1", "results", "level2"))
    p.add_argument("--biosnicar", default=None)
    p.add_argument("--workers", type=int, default=2)
    p.add_argument("--n-per-sza", type=int, default=24)
    p.add_argument("--outdir", default=os.path.join(HERE, "results", "surrogate"))
    a = p.parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)
    ems = build(a.optics, a.sza, a.phase1_l2, a.biosnicar, a.workers, a.outdir)
    rec, tab = benchmark(ems, a.optics, a.n_per_sza, phase1_l2=a.phase1_l2, biosnicar=a.biosnicar)
    print(json.dumps(rec, indent=1))
    import provenance as PV
    PV.atomic_write_text(os.path.join(a.outdir, f"benchmark_{a.optics}.json"), json.dumps(rec, indent=1))
    PV.atomic_to_csv(tab, os.path.join(a.outdir, f"benchmark_{a.optics}_points.csv"), index=False)


if __name__ == "__main__":
    main()
