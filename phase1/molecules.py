"""
Molecular definitions for Phase 1 (Level 1 and Level 2) glacier-algal pigment models.

Level 1  : Purpurogallin                                     C11H8O5   (M = 220.18 g/mol)
           2,3,4,6-tetrahydroxy-5H-benzo[7]annulen-5-one (benzotropolone core)
Level 2  : Purpurogallin carboxylic acid-6-O-beta-D-glucopyranoside   C18H18O12 (M = 426.33 g/mol)
           (Remias et al. 2012, FEMS Microbiol. Ecol. 79:638; Ancylonema / Mesotaenium)

Benzo[7]annulene numbering used throughout:

            OH(4)   O(5)
             |      ||
     HO(3)-C4--C4a--C5--C6-OR(6)        R = H (Level 1) / beta-D-Glcp (Level 2)
           |    |         \\
     HO(2)-C3   |          C7-H
            \\  |          |
         C2==C1-C9a==C9---C8-X          X = H (Level 1) / COOH (Level 2)

Both systems are closed-shell neutral singlets:
    charge = 0, multiplicity = 1   (PySCF: mol.charge = 0, mol.spin = 2S = 0)

STRUCTURAL ASSUMPTION (Level 2): the original identification places the carboxyl
group on the tropolone ring; we use C8 (the regiochemistry produced by the
gallic-acid + pyrogallol oxidative coupling route). If your NMR/MS assignment
differs, change LEVEL2_SMILES only - every downstream step is built from it.

Coordinates are given two ways:
  1. build_xyz(): programmatic RDKit build (ETKDGv3 conformer search + MMFF94s),
     the recommended route so the conformer search can be re-run / extended.
  2. LEVEL1_XYZ / LEVEL2_XYZ: hard-coded Cartesian coordinates (Angstrom) of the
     lowest-MMFF94s conformer produced by build_xyz() with the default seed, so
     the pipeline runs even without RDKit.
"""

from __future__ import annotations

import numpy as np

# --------------------------------------------------------------------------- #
# SMILES (stereochemistry verified: the glucosyl fragment reproduces the      #
# InChIKey of arbutin, a known 4-hydroxyphenyl beta-D-glucopyranoside)         #
# --------------------------------------------------------------------------- #
_BETA_D_GLCP = "[C@@H]9O[C@H](CO)[C@@H](O)[C@H](O)[C@H]9O"   # attached via anomeric C1'

LEVEL1_SMILES = "C1=CC2=CC(O)=C(O)C(O)=C2C(=O)C(O)=C1"
LEVEL2_SMILES = "OC(=O)C1=CC2=CC(O)=C(O)C(O)=C2C(=O)C(O" + _BETA_D_GLCP + ")=C1"

# Hard-coded starting geometries (Angstrom). Filled in below from build_xyz().
LEVEL1_XYZ = """
C      2.501024    -1.543779    -0.119046
C      1.210068    -1.891919    -0.122920
C     -0.073667    -1.009734     0.031322
C     -1.233225    -1.688975     0.459108
C     -2.429185    -1.006196     0.626133
O     -3.537959    -1.682315     1.041429
C     -2.494921     0.349940     0.367920
O     -3.685780     1.003867     0.530209
C     -1.362395     1.034658    -0.056767
O     -1.531897     2.375145    -0.310656
C     -0.134717     0.367249    -0.218244
C      1.106983     1.149832    -0.632482
O      0.929248     2.216308    -1.208723
C      2.573996     0.862105    -0.298029
O      3.422993     1.965693    -0.304970
C      3.101956    -0.350594    -0.130919
H      3.193542    -2.391394    -0.090657
H      1.029630    -2.969150    -0.149905
H     -1.211719    -2.756019     0.668215
H     -4.242150    -1.003289     1.081771
H     -3.483021     1.932829     0.291434
H     -0.723770     2.699864    -0.766936
H      2.885179     2.705005    -0.659093
H      4.189786    -0.369132    -0.028195
"""

LEVEL2_XYZ = """
O      1.398697    -4.782659    -2.028392
C      2.466064    -3.967342    -1.926555
O      3.560684    -4.338431    -2.314652
C      2.133147    -2.618606    -1.354125
C      3.012681    -1.660749    -1.711696
C      3.215673    -0.241152    -1.048622
C      4.535018     0.220306    -0.915705
C      4.777909     1.435756    -0.293873
O      6.060979     1.878301    -0.159542
C      3.725275     2.191528     0.199245
O      3.998555     3.379837     0.824304
C      2.410630     1.752423     0.068347
O      1.460911     2.566842     0.618825
C      2.152906     0.534162    -0.559381
C      0.767086    -0.063320    -0.678685
O     -0.026307     0.479110    -1.432469
C      0.411148    -1.455330    -0.108862
O     -0.620833    -1.610326     0.813651
C     -1.768958    -0.785460     0.571579
O     -1.484284     0.583430     0.847618
C     -2.586019     1.421493     0.464960
C     -2.104076     2.872728     0.495878
O     -0.960438     3.025660    -0.362684
C     -3.764918     1.145934     1.400305
O     -4.903313     1.925581     1.036614
C     -4.124433    -0.346047     1.322230
O     -5.122205    -0.639445     2.312374
C     -2.895332    -1.242695     1.520456
O     -3.283496    -2.604745     1.283988
C      1.060109    -2.548343    -0.541874
H      1.767366    -5.592468    -2.440745
H      3.788678    -1.925188    -2.434225
H      5.375548    -0.364716    -1.280003
H      5.978424     2.736934     0.303836
H      3.115462     3.714079     1.090763
H      0.589276     2.454944     0.173228
H     -2.082489    -0.931244    -0.472254
H     -2.860948     1.195500    -0.574937
H     -2.880990     3.563893     0.156653
H     -1.784889     3.160818     1.502225
H     -1.127257     2.499646    -1.172020
H     -3.504576     1.404338     2.434259
H     -5.641188     1.572614     1.575170
H     -4.598790    -0.564739     0.356799
H     -5.216726    -1.613603     2.302607
H     -2.562857    -1.189579     2.564620
H     -2.470146    -3.139065     1.369095
H      0.613242    -3.490604    -0.218273
"""

MOLECULES = {
    "level1": dict(
        name="Purpurogallin",
        smiles=LEVEL1_SMILES,
        formula="C11H8O5",
        molar_mass=220.180,      # g/mol
        charge=0,
        multiplicity=1,
        xyz=LEVEL1_XYZ,
        n_confs=60,
    ),
    "level2": dict(
        name="Purpurogallin carboxylic acid-6-O-beta-D-glucopyranoside",
        smiles=LEVEL2_SMILES,
        formula="C18H18O12",
        molar_mass=426.330,      # g/mol
        charge=0,
        multiplicity=1,
        xyz=LEVEL2_XYZ,
        n_confs=300,             # flexible glycoside: sample more conformers
    ),
    # Level 3: Fe(III)-purpurogallin model complexes (fe_complex.py); geometry built on demand
    "level3_catecholate": dict(
        name="[Fe(III)(purpurogallin-catecholate)(H2O)4]+",
        smiles=None, formula="C11H14FeO9", molar_mass=380.08, charge=1, multiplicity=6, xyz=None,
        ligand_molar_mass=220.180, n_confs=1,
    ),
    "level3_tropolonate": dict(
        name="[Fe(III)(purpurogallin-tropolonate)(H2O)4]2+",
        smiles=None, formula="C11H15FeO9", molar_mass=381.09, charge=2, multiplicity=6, xyz=None,
        ligand_molar_mass=220.180, n_confs=1,
    ),
}


def build_xyz(smiles: str, n_confs: int = 100, seed: int = 20240611,
              expected_formula: str | None = None) -> tuple[str, float]:
    """Build 3D coordinates from SMILES with RDKit.

    ETKDGv3 embeds `n_confs` conformers, each is relaxed with MMFF94s, and the
    lowest-energy one is returned as a PySCF-style XYZ block ("El x y z" lines).

    Returns (xyz_block, mmff_energy_kcal_per_mol).
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula

    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if expected_formula is not None:
        got = CalcMolFormula(mol)
        if got != expected_formula:
            raise ValueError(f"SMILES gives {got}, expected {expected_formula}")

    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.pruneRmsThresh = 0.25
    cids = list(AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, params=params))
    if not cids:
        raise RuntimeError("RDKit failed to embed any conformer")

    mp = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
    results = AllChem.MMFFOptimizeMoleculeConfs(mol, mmffVariant="MMFF94s", maxIters=5000)
    energies = np.array([e if conv == 0 else np.inf for conv, e in results])
    best = cids[int(np.argmin(energies))]
    del mp

    conf = mol.GetConformer(best)
    lines = []
    for atom in mol.GetAtoms():
        p = conf.GetAtomPosition(atom.GetIdx())
        lines.append(f"{atom.GetSymbol():2s} {p.x:12.6f} {p.y:12.6f} {p.z:12.6f}")
    return "\n".join(lines), float(energies.min())


def get_xyz(key: str, rebuild: bool = False) -> str:
    """Return the starting XYZ block for 'level1' or 'level2'."""
    spec = MOLECULES[key]
    if key.startswith("level3_"):
        import fe_complex
        return fe_complex.build_xyz(key.split("_", 1)[1])
    if rebuild:
        xyz, _ = build_xyz(spec["smiles"], spec["n_confs"], expected_formula=spec["formula"])
        return xyz
    return spec["xyz"].strip()


def get_mole(key: str, basis: str = "6-31g*", xyz: str | None = None,
             cart: bool = True, verbose: int = 4, max_memory: int = 10000):
    """Build a pyscf.gto.Mole for 'level1' / 'level2'.

    cart=True reproduces the Gaussian convention for 6-31G(d) (six Cartesian d
    functions), so energies/spectra are directly comparable with Gaussian results.
    """
    from pyscf import gto

    spec = MOLECULES[key]
    mol = gto.Mole()
    mol.atom = xyz if xyz is not None else get_xyz(key)
    mol.unit = "Angstrom"
    mol.basis = basis
    mol.cart = cart
    mol.charge = spec["charge"]
    mol.spin = spec["multiplicity"] - 1          # PySCF spin = N_alpha - N_beta = 2S
    mol.verbose = verbose
    mol.max_memory = max_memory                  # MB
    mol.build()

    # Sanity check: elemental composition must match the target formula.
    from collections import Counter
    counts = Counter(mol.atom_pure_symbol(i) for i in range(mol.natm))
    formula = "".join(f"{el}{counts[el] if counts[el] > 1 else ''}" for el in ("C", "H", "Fe", "O") if counts[el])
    if formula != spec["formula"]:
        raise ValueError(f"{key}: geometry has formula {formula}, expected {spec['formula']}")
    return mol


if __name__ == "__main__":
    # Regenerate and print starting geometries (used to fill LEVEL*_XYZ above).
    for key, spec in MOLECULES.items():
        xyz, e = build_xyz(spec["smiles"], spec["n_confs"], expected_formula=spec["formula"])
        n = len(xyz.splitlines())
        print(f"# {key}: {spec['name']}  ({spec['formula']}, {n} atoms, "
              f"charge {spec['charge']}, multiplicity {spec['multiplicity']}, E_MMFF94s = {e:.3f} kcal/mol)")
        print(xyz)
        print()
