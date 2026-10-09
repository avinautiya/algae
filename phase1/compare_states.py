"""
Approximation checks for the Level 2 excited states (repair item 3), with states matched by overlaps,
not by index.

Comparisons (each is a measurement; acceptance targets apply only to numerical approximations):
  A   B3LYP full TD (30 roots)  vs  B3LYP TDA (15 roots)          TDA approximation (physical model change)
  B   B3LYP TDA (15)            vs  CAM-B3LYP TDA (15)             functional effect at equal approximation
  R   B3LYP TDA (15)            vs  B3LYP TDA (25)                 root-count convergence of the low roots
  T   per run: stage 1e-3 -> 1e-4 -> 1e-5                          solver-tolerance convergence
  G   B3LYP TDA (15) at the start geometry vs at the relaxed one   geometry sensitivity

Numerical acceptance (R, T): every matched bright state (f >= BRIGHT_F) changes by < 0.01 eV and the
broadened absorption maximum by < 2 nm. A and B are reported, not accepted or rejected.

State matching:
  * both runs carry transition vectors (td_ckpt/*/latest.npz from tdcheckpoint): overlap of AO-basis
    transition densities T = C_occ X C_vir^T, |<T_i|S|T_j>| / (|T_i||T_j|), with the AO overlap S
    recovered from the orbitals (S = (C C^T)^-1; C is square for these runs). For different geometries
    (G) the AO overlap of the first run is used: an approximation for small displacements, reported.
  * otherwise (the legacy full-TD run kept no vectors): cosine similarity of the dominant-transition
    weights in the states table ("HOMO-1->LUMO (0.20)"). Components below 0.1 are not stored, so this
    matching is coarser and is labelled as such.
  Assignment: Hungarian algorithm on 1 - overlap; a pair with overlap < MIN_OVERLAP is reported as
  unmatched rather than forced.

Usage: python compare_states.py [--out results/v2/approximation_checks.json]
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import re
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import completion  # noqa: E402
import spectra  # noqa: E402

BRIGHT_F = 0.05
MIN_OVERLAP = 0.5
TOL_DE_EV, TOL_PEAK_NM = 0.01, 2.0
V2 = os.path.join(HERE, "results", "v2")


class Run:
    def __init__(self, outdir: str, func: str, level: str = "level2", allow_provisional: bool = True):
        self.outdir, self.func = outdir, func
        ok, st, why = completion.validate_run(outdir, allow_provisional=allow_provisional)
        self.valid, self.status, self.reasons = ok, st, why
        stem = os.path.join(outdir, f"{level}_{func}")
        df = pd.read_csv(stem + "_states.csv")
        self.e, self.f = df.Energy_eV.to_numpy(float), df.Oscillator_Strength.to_numpy(float)
        self.labels = df.get("Transitions", pd.Series([""] * len(df))).fillna("").tolist()
        summ = json.load(open(os.path.join(outdir, "summary.json")))
        self.molar_mass = float(summ["molar_mass_g_mol"])
        self.tda = bool(summ.get("tda", False))
        hist = os.path.join(outdir, f"{level}_{func}_stage_history.json")
        self.stages = json.load(open(hist)).get("stages", []) if os.path.isfile(hist) else []
        self.vec = self.C = self.nocc = None
        tag = f"{func}_{'TDA' if self.tda else 'RPA'}_{len(self.e)}"
        for p in glob.glob(os.path.join(outdir, "td_ckpt", tag, "latest.npz")):
            with np.load(p) as z:
                if int(z["nspin"]) == 1 and z["vectors"].shape[0] == len(self.e):
                    self.vec, self.C = np.asarray(z["vectors"], float), np.asarray(z["mo_coeff_0"], float)
                    self.nocc = int(np.sum(np.asarray(z["mo_occ_0"]) > 0))

    def transition_densities(self):
        C = self.C
        nocc, nvir = self.nocc, C.shape[1] - self.nocc
        out = []
        for v in self.vec:
            x = v[:nocc * nvir].reshape(nocc, nvir)
            if not self.tda:                                  # RPA vectors are [X, Y]: use X + Y
                x = x + v[nocc * nvir:2 * nocc * nvir].reshape(nocc, nvir)
            out.append(C[:, :nocc] @ x @ C[:, nocc:].T)
        return out


def _label_vector(label: str) -> dict:
    out = {}
    for h, p, w in re.findall(r"(HOMO[ab]?(?:-\d+)?)->(LUMO[ab]?(?:\+\d+)?) \(([-0-9.]+)\)", label):
        out[(h, p)] = out.get((h, p), 0.0) + float(w)
    return out


def overlap_matrix(a: Run, b: Run):
    """(|overlap| matrix, method)."""
    if a.vec is not None and b.vec is not None and a.C.shape == b.C.shape:
        S = np.linalg.inv(a.C @ a.C.T)
        Ta, Tb = a.transition_densities(), b.transition_densities()
        na = [np.sqrt(abs(np.sum((S @ t @ S) * t))) for t in Ta]
        nb = [np.sqrt(abs(np.sum((S @ t @ S) * t))) for t in Tb]
        O = np.array([[abs(np.sum((S @ ta @ S) * tb)) / (x * y) for tb, y in zip(Tb, nb)] for ta, x in zip(Ta, na)])
        return O, "transition-density overlap"
    La, Lb = [_label_vector(s) for s in a.labels], [_label_vector(s) for s in b.labels]

    def cos(u, v):
        k = set(u) | set(v)
        x, y = np.array([u.get(i, 0.0) for i in k]), np.array([v.get(i, 0.0) for i in k])
        d = np.linalg.norm(x) * np.linalg.norm(y)
        return float(x @ y / d) if d > 0 else 0.0
    return np.array([[cos(u, v) for v in Lb] for u in La]), "dominant-transition label cosine (coarse)"


def match(a: Run, b: Run):
    from scipy.optimize import linear_sum_assignment
    O, how = overlap_matrix(a, b)
    r, c = linear_sum_assignment(1.0 - O)
    pairs = [(int(i), int(j), float(O[i, j])) for i, j in zip(r, c)]
    return pairs, how


def peak_nm(e, f, mm, fwhm, lo=300.0, hi=750.0):
    lam = np.arange(lo, hi + 0.05, 0.1)
    m = spectra.epsilon_to_mac(spectra.gaussian_broaden(spectra.nm_to_ev(lam), e, f, fwhm), mm)
    return float(lam[int(np.argmax(m))])


def compare(a: Run, b: Run, name: str, kind: str, n_compare: int | None = None) -> dict:
    """kind: 'numerical' (acceptance targets apply) or 'model' (measured difference only)."""
    n = min(len(a.e), len(b.e)) if n_compare is None else n_compare
    pairs, how = match(a, b)
    rows = []
    for i, j, o in sorted(pairs):
        if i >= n:
            continue
        rows.append(dict(state_a=i + 1, state_b=j + 1, overlap=round(o, 4), e_a=a.e[i], e_b=b.e[j],
                         dE_ev=b.e[j] - a.e[i], f_a=a.f[i], f_b=b.f[j],
                         bright=bool(max(a.f[i], b.f[j]) >= BRIGHT_F), matched=bool(o >= MIN_OVERLAP)))
    bright = [r for r in rows if r["bright"]]
    # peak over the low roots common to both (so a larger root set does not move the peak by itself)
    k = min(len(a.e), len(b.e))
    peaks = {str(w): (peak_nm(a.e[:k], a.f[:k], a.molar_mass, w), peak_nm(b.e[:k], b.f[:k], b.molar_mass, w))
             for w in (0.3, 0.6)}
    out = dict(name=name, kind=kind, a=a.outdir, b=b.outdir, a_status=a.status, b_status=b.status,
               matching=how, states=rows,
               max_abs_dE_bright_ev=max((abs(r["dE_ev"]) for r in bright), default=None),
               unmatched_bright=[r["state_a"] for r in bright if not r["matched"]],
               peak_nm={w: dict(a=p[0], b=p[1], shift=p[1] - p[0]) for w, p in peaks.items()})
    if kind == "numerical":
        out["accept"] = bool(out["max_abs_dE_bright_ev"] is not None and out["max_abs_dE_bright_ev"] < TOL_DE_EV
                             and not out["unmatched_bright"]
                             and all(abs(v["shift"]) < TOL_PEAK_NM for v in out["peak_nm"].values()))
    return out


def tolerance_check(run: Run) -> dict | None:
    """T: energies/f recorded when every root had converged at 1e-3, 1e-4, 1e-5 (index order within one
    run; reorderings are flagged by a change in the f pattern)."""
    st = {s["tol"]: s for s in run.stages if s.get("energies_au") is not None}
    if not st:
        return None
    tols = sorted(st, reverse=True)
    final = np.asarray(st[tols[-1]]["energies_au"]) * 27.211386245988
    ff = np.asarray(st[tols[-1]]["osc"]) if st[tols[-1]].get("osc") is not None else run.f
    br = ff >= BRIGHT_F
    rows = []
    for t in tols:
        e = np.asarray(st[t]["energies_au"]) * 27.211386245988
        f = np.asarray(st[t]["osc"]) if st[t].get("osc") is not None else None
        pk = (peak_nm(e, f, run.molar_mass, 0.3), peak_nm(final, ff, run.molar_mass, 0.3)) if f is not None else None
        rows.append(dict(tol=t, max_abs_dE_bright_vs_final_ev=float(np.max(np.abs(e - final)[br])) if br.any() else None,
                         max_abs_dE_all_vs_final_ev=float(np.max(np.abs(e - final))),
                         peak_shift_vs_final_nm=None if pk is None else pk[0] - pk[1]))
    return dict(run=run.outdir, stages=rows, note="relative to the final (tightest) tolerance; the final tolerance's "
                "own error is not bounded by this comparison")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(V2, "approximation_checks.json"))
    a = ap.parse_args(argv)
    runs = {}
    for key, d, func in [("FULL30", os.path.join(HERE, "results", "level2"), "B3LYP"),
                         ("TDA15", os.path.join(V2, "L2_B3LYP_TDA15"), "B3LYP"),
                         ("TDA25", os.path.join(V2, "L2_B3LYP_TDA25"), "B3LYP"),
                         ("CAM15", os.path.join(V2, "L2_CAM_TDA15"), "CAM-B3LYP"),
                         ("RELAX15", os.path.join(V2, "L2_B3LYP_TDA15_RELAXED"), "B3LYP"),
                         ("COO15", os.path.join(V2, "L2_COO_OPT_TDA15"), "B3LYP")]:
        try:
            r = Run(d, func, level="level2_carboxylate" if key == "COO15" else "level2")
        except (OSError, KeyError) as e:
            print(f"{key}: not available ({e.__class__.__name__})")
            continue
        if r.status not in ("production", "provisional_geometry"):
            print(f"{key}: status {r.status} - not compared ({'; '.join(r.reasons)})")
            continue
        runs[key] = r
    res = dict(bright_f=BRIGHT_F, min_overlap=MIN_OVERLAP, targets=dict(dE_ev=TOL_DE_EV, peak_nm=TOL_PEAK_NM),
               comparisons=[], tolerance=[])
    for name, x, y, kind in [("A: B3LYP full vs TDA", "FULL30", "TDA15", "model"),
                             ("B: B3LYP TDA vs CAM-B3LYP TDA", "TDA15", "CAM15", "model"),
                             ("R: B3LYP TDA 15 vs 25 roots", "TDA15", "TDA25", "numerical"),
                             ("G: start vs relaxed geometry (B3LYP TDA 15)", "TDA15", "RELAX15", "model"),
                             ("P: neutral acid vs carboxylate (B3LYP TDA 15)", "TDA15", "COO15", "model")]:
        if x in runs and y in runs:
            c = compare(runs[x], runs[y], name, kind, n_compare=15)
            res["comparisons"].append(c)
            acc = f", accept={c['accept']}" if "accept" in c else ""
            print(f"{name}: max |dE| bright {c['max_abs_dE_bright_ev']}, peak shift "
                  f"{ {w: round(v['shift'], 1) for w, v in c['peak_nm'].items()} } nm, unmatched {c['unmatched_bright']}"
                  f" [{c['matching']}]{acc}")
        else:
            print(f"{name}: pending (needs {x} and {y})")
    for k, r in runs.items():
        t = tolerance_check(r)
        if t:
            res["tolerance"].append(t)
            print(f"T {k}: " + ", ".join(f"{s['tol']:g}: {s['max_abs_dE_bright_vs_final_ev']}" for s in t["stages"]))
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    tmp = a.out + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(res, fh, indent=1, default=float)
    os.replace(tmp, a.out)
    return res


if __name__ == "__main__":
    main()
