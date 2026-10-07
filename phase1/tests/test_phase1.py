"""
Validation tests for Phase 1 that do not need a DFT run.   Run:  python phase1/tests/test_phase1.py  (or pytest)

1. Gaussian broadening conserves oscillator strength: integral of eps dE = 2.8707e4 * sum(f).
2. eps -> MAC conversion (Napierian and decadic) and the nm <-> eV relations.
3. Both molecules' built-in geometries have the target formula, charge and multiplicity (no SCF).
4. A tiny real PySCF calculation (H2O, sto-3g) runs through make_mf + run_tddft (PCM off).
"""

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))

import spectra as S  # noqa: E402


def test_broadening_conserves_oscillator_strength():
    e = np.linspace(0.5, 10.0, 20001)
    lines, f = np.array([3.0, 4.2, 5.1]), np.array([0.2, 0.5, 0.05])
    eps = S.gaussian_broaden(e, lines, f, fwhm_ev=0.3)
    assert abs(np.trapezoid(eps, e) / (S.EPS_INTEGRAL_K * f.sum()) - 1) < 1e-6
    assert abs(S.EPS_INTEGRAL_K - 2.8707e4) < 5.0                      # Hilborn (1982)
    # the peak position is the line position
    one = S.gaussian_broaden(e, np.array([3.0]), np.array([1.0]), 0.3)
    assert abs(e[np.argmax(one)] - 3.0) < 1e-3


def test_units():
    assert abs(S.nm_to_ev(500.0) - 2.4796839) < 1e-6 and abs(S.ev_to_nm(S.nm_to_ev(413.0)) - 413.0) < 1e-9
    eps = np.array([10000.0])                                          # L mol^-1 cm^-1
    dec = S.epsilon_to_mac(eps, 426.33, napierian=False)
    nap = S.epsilon_to_mac(eps, 426.33, napierian=True)
    # decadic MAC [m^2 kg^-1] = eps * 0.1 / M[kg/mol]^-1 ... = eps * 1000 * 0.1 / M
    assert abs(dec[0] - 10000.0 * 100.0 / 426.33) < 1e-6
    assert abs(nap[0] / dec[0] - np.log(10.0)) < 1e-12


def test_molecule_definitions():
    import molecules as M
    from rdkit import Chem
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula
    for key, formula in (("level1", "C11H8O5"), ("level2", "C18H18O12")):
        spec = M.MOLECULES[key]
        assert spec["formula"] == formula
        assert CalcMolFormula(Chem.MolFromSmiles(spec["smiles"])) == formula
        mol = M.get_mole(key, basis="sto-3g")                          # builds the Mole, checks the formula
        assert mol.charge == 0 and mol.spin == 0


def test_tiny_tddft_runs():
    from pyscf import gto
    import qc
    mol = gto.M(atom="O 0 0 0; H 0 0.757 0.587; H 0 -0.757 0.587", basis="sto-3g", verbose=0)
    mf = qc.make_mf(mol, "B3LYP", solvent=None)
    mf.kernel()
    assert mf.converged
    out = qc.run_tddft(mf, nstates=3)
    e = np.asarray(out["energies_ev"])
    assert e.size >= 3 and np.all(np.diff(e) >= -1e-9) and e[0] > 5.0        # water's first band is UV


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print(f"PASS {name}")
