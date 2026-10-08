"""
Level 3 model of the Fe(III)-purpurogallin pigment complex (computed counterpart of the measured
Fe-purpurogallin spectrum of Prochazkova et al. 2025, Fig. 4, which was taken on the aglycone).

Model: one purpurogallin (Level 1 core) bound bidentately to high-spin Fe(III), the remaining four
octahedral sites filled by water, in IEF-PCM water:

  catecholate  [Fe(PG-2H)(H2O)4]+   two adjacent phenolates of the trihydroxy-benzo ring (the galloyl
               site that binds Fe(III) in gallates/pyrogallols); charge +1, multiplicity 6 (S = 5/2)
  tropolonate  [Fe(PG-H)(H2O)4]2+   carbonyl O + adjacent tropolone O(-); charge +2, multiplicity 6

Which mode the algal pigment adopts is not known (Prochazkova et al. 2025 show complexation by Raman,
not the binding site), so the catecholate model is the main (DFT-optimised) structure and the
tropolonate model a sensitivity case. The ground spin state is checked by single points at
multiplicity 6, 4 and 2 (spin_check).

Starting structures are built here (RDKit purpurogallin + idealised octahedral Fe-O 2.0 A, Fe-OH2
2.1 A) and relaxed with GFN2-xTB/ALPB(water) in run_phase1.py before DFT.
"""

from __future__ import annotations

import numpy as np

FE_O = 2.00          # A, Fe(III)-phenolate
FE_OW = 2.10         # A, Fe(III)-OH2
OH = 0.97            # A
HOH = np.radians(104.5)

MODES = {
    "catecholate": dict(charge=1, multiplicity=6, formula="C11H14FeO9"),
    "tropolonate": dict(charge=2, multiplicity=6, formula="C11H15FeO9"),
}


def _pg_rdkit(seed=20240611):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from molecules import LEVEL1_SMILES
    m = Chem.AddHs(Chem.MolFromSmiles(LEVEL1_SMILES))
    p = AllChem.ETKDGv3()
    p.randomSeed = seed
    AllChem.EmbedMolecule(m, p)
    AllChem.MMFFOptimizeMolecule(m, mmffVariant="MMFF94s", maxIters=5000)
    return m


def _binding_oxygens(m, mode):
    """Indices (O1, O2) of the two donor oxygens and the hydrogens to remove."""
    def h_of(o):
        return [n.GetIdx() for n in m.GetAtomWithIdx(o).GetNeighbors() if n.GetSymbol() == "H"]

    hydroxyl = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == "O" and h_of(a.GetIdx())]
    carbonyl = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == "O" and not h_of(a.GetIdx())]
    c_of = {o: next(n.GetIdx() for n in m.GetAtomWithIdx(o).GetNeighbors() if n.GetSymbol() == "C")
            for o in hydroxyl + carbonyl}
    bonded = lambda a, b: m.GetBondBetweenAtoms(a, b) is not None  # noqa: E731
    if mode == "catecholate":
        # adjacent phenolic OH pairs on the aromatic (benzo) ring; take the pair farthest from the
        # carbonyl, whose neighbouring OH is engaged in the O-H...O=C hydrogen bond
        pos = m.GetConformer().GetPositions()
        pairs = [(a, b) for i, a in enumerate(hydroxyl) for b in hydroxyl[i + 1:]
                 if m.GetAtomWithIdx(c_of[a]).GetIsAromatic() and m.GetAtomWithIdx(c_of[b]).GetIsAromatic()
                 and bonded(c_of[a], c_of[b])]
        co = carbonyl[0]
        pair = max(pairs, key=lambda p: min(np.linalg.norm(pos[p[0]] - pos[co]), np.linalg.norm(pos[p[1]] - pos[co])))
        return pair, h_of(pair[0]) + h_of(pair[1])
    if mode == "tropolonate":
        co = carbonyl[0]
        o = next(o for o in hydroxyl if bonded(c_of[o], c_of[co]))
        return (co, o), h_of(o)
    raise ValueError(mode)


def build_xyz(mode: str = "catecholate") -> str:
    """Idealised starting geometry (Angstrom, 'El x y z' lines): Fe first, then ligand, then waters."""
    m = _pg_rdkit()
    (o1, o2), drop = _binding_oxygens(m, mode)
    pos = m.GetConformer().GetPositions()
    sym = [a.GetSymbol() for a in m.GetAtoms()]
    keep = [i for i in range(m.GetNumAtoms()) if i not in drop]
    p1, p2 = pos[o1], pos[o2]
    mid = 0.5 * (p1 + p2)
    # in-plane direction away from the ligand: from the ring centroid through the O-O midpoint
    ring_c = pos[[i for i in range(len(sym)) if sym[i] == "C"]].mean(axis=0)
    u = mid - ring_c
    v = (p1 - p2) / np.linalg.norm(p1 - p2)
    u -= u.dot(v) * v
    u /= np.linalg.norm(u)
    half = 0.5 * np.linalg.norm(p1 - p2)
    fe = mid + u * np.sqrt(max(FE_O ** 2 - half ** 2, 0.5))
    d1 = (p1 - fe) / np.linalg.norm(p1 - fe)
    d2 = (p2 - fe) / np.linalg.norm(p2 - fe)
    w = np.cross(d1, d2)
    w /= np.linalg.norm(w)
    lines = [f"Fe {fe[0]:12.6f} {fe[1]:12.6f} {fe[2]:12.6f}"]
    lines += [f"{sym[i]:2s} {pos[i][0]:12.6f} {pos[i][1]:12.6f} {pos[i][2]:12.6f}" for i in keep]
    for d in (-d1, -d2, w, -w):
        o = fe + FE_OW * d
        # water hydrogens pointing away from Fe, in a plane containing the Fe-O axis
        perp = np.cross(d, w if abs(d.dot(w)) < 0.9 else d1)
        perp /= np.linalg.norm(perp)
        for s in (+1, -1):
            h = o + OH * (np.cos(HOH / 2) * d + s * np.sin(HOH / 2) * perp)
            lines.append(f"H  {h[0]:12.6f} {h[1]:12.6f} {h[2]:12.6f}")
        lines.append(f"O  {o[0]:12.6f} {o[1]:12.6f} {o[2]:12.6f}")
    return "\n".join(lines)


def spin_check(xyz: str, mode: str = "catecholate", functional="B3LYP", basis="6-31g*", max_memory=7000):
    """UKS/PCM single points at multiplicity 6, 4 and 2 on one geometry. Returns {mult: (E_Eh, <S^2>)}."""
    from pyscf import gto
    import qc
    out = {}
    for mult in (6, 4, 2):
        mol = gto.M(atom=xyz, basis=basis, cart=True, charge=MODES[mode]["charge"], spin=mult - 1,
                    max_memory=max_memory, verbose=3)
        mf = qc.make_mf(mol, functional, "pcm")
        e = mf.kernel()
        s2, _ = mf.spin_square()
        out[mult] = (float(e), float(s2), bool(mf.converged))
        print(f"[spin] {mode} M={mult}: E = {e:.6f} Eh, <S^2> = {s2:.3f} (ideal {(mult - 1) / 2 * ((mult - 1) / 2 + 1):.2f}), "
              f"converged {mf.converged}", flush=True)
    return out


if __name__ == "__main__":
    import sys
    print(build_xyz(sys.argv[1] if len(sys.argv) > 1 else "catecholate"))
