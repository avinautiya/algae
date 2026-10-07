# Phase 1: DFT/TD-DFT modelling of glacier-algal pigments (Levels 1 and 2)

| | Level 1 | Level 2 |
|---|---|---|
| Compound | Purpurogallin (benzotropolone core) | Purpurogallin carboxylic acid-6-O-β-D-glucopyranoside |
| Formula / M | C11H8O5 / 220.18 g mol⁻¹ | C18H18O12 / 426.33 g mol⁻¹ |
| SMILES | `C1=CC2=CC(O)=C(O)C(O)=C2C(=O)C(O)=C1` | `OC(=O)C1=CC2=CC(O)=C(O)C(O)=C2C(=O)C(O[C@@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)=C1` |
| Atoms / electrons | 24 / 114 | 48 / 222 |
| Basis functions (6-31G(d), 6d) | 256 | 486 |
| Charge, multiplicity | 0, singlet (PySCF `spin=0`) | 0, singlet (PySCF `spin=0`) |

## Files

| File | Purpose |
|---|---|
| `molecules.py` | SMILES, charge/multiplicity, RDKit builder (ETKDGv3 + MMFF94s conformer search), hard-coded starting Cartesians, formula checks |
| `qc.py` | RKS (B3LYP / CAM-B3LYP, 6-31G(d)), IEF-PCM or SMD water, xTB pre-optimisation, geomeTRIC optimisation with per-step checkpoints, non-equilibrium TD-DFT |
| `spectra.py` | Gaussian broadening in energy, ε(λ) → MAC(λ), DataFrames, plotting (NumPy/SciPy/pandas only) |
| `run_phase1.py` | End-to-end command-line driver |
| `phase1_colab.ipynb` | Colab notebook (install → Level 1 → Level 2 → plots) |

## Method

1. **Start geometry:** RDKit ETKDGv3 conformer search (60 conformers for Level 1, 300 for Level 2) with MMFF94s. The lowest conformer is stored in `molecules.py`.
2. **Pre-optimisation:** GFN2-xTB with ALPB(water) through `tblite`. This takes seconds and cuts the number of DFT steps.
3. **Ground state:** RKS B3LYP/6-31G(d) with IEF-PCM (ε = 78.4), optimised by geomeTRIC to Gaussian-default thresholds.
   - `B3LYP` maps to PySCF `B3LYPG` (VWN-RPA), which matches Gaussian's B3LYP.
   - Cartesian 6d functions also match Gaussian's 6-31G(d).
   - Density fitting (def2-universal-jkfit) is on by default; use `--no-df` to turn it off.
   - `--solvent smd` switches to SMD(water).
4. **Excited states:** full TD-DFT (RPA), 30 singlet roots, with **non-equilibrium** linear-response PCM. The electronic response uses the optical dielectric constant ε∞ = 1.78, which is correct for vertical absorption.
   - Pass `--tddft-functionals B3LYP CAM-B3LYP` to get both spectra at the same B3LYP geometry.
   - That gives the CAM-B3LYP benchmark without a second optimisation. If you want a CAM-B3LYP geometry instead, use `--opt-functional CAM-B3LYP`.
5. **Spectrum and MAC:**
   - Each line is broadened with an area-normalised Gaussian in energy: ε(E) = 2.8707×10⁴ Σᵢ fᵢ gᵢ(E) L mol⁻¹ cm⁻¹ (Hilborn 1982).
   - FWHM is 0.3 eV, which is about 41 nm at 413 nm.
   - MAC [m² kg⁻¹] = 1000 · ln(10) · 0.1 · ε / M (Napierian, the BioSNICAR convention; `--decadic` drops the ln 10).
   - The spectrum is evaluated on a 1 nm grid from 300 to 800 nm.

## Outputs (`results/<level>/`)

- `<level>_<F>_states.csv`: one row per state, with columns `State, Wavelength_nm, Energy_eV, Oscillator_Strength, MAC_estimated, Transitions`. `MAC_estimated` is the broadened MAC at that state's λ.
- `<level>_<F>_spectrum.csv`: the 300–800 nm grid, with columns `Wavelength_nm, Energy_eV, Oscillator_Strength, Epsilon_L_mol-1_cm-1, MAC_estimated`. Here `Oscillator_Strength` is the broadened f-density (eV⁻¹).
- `<level>_<F>_spectrum.npz` / `.png`: plot-ready arrays and a figure.
- `opt_B3LYP_traj.xyz`, `_last.xyz`, `_final.xyz`: the optimisation trajectory, a resume point, and the final structure.
- `summary.json`

## Runtime

Measured on a 4-vCPU Xeon (2.8 GHz) for Level 1 at B3LYP/6-31G(d)/IEF-PCM with density fitting:

| Step | Measured cost |
|---|---|
| SCF from scratch | 140–235 s |
| SCF during optimisation (restarts from the previous density) | about 60–80 s |
| Analytic gradient | about 45 s |
| TD-DFT Davidson iteration (10 roots) | about 47 s; split roughly evenly between XC kernel, PCM response and DF-JK |
| `--grid-level 2 --lebedev 17` | SCF 74 s, gradient 27 s, TD iteration 34 s (about 1.5–2× faster) |

The full-job estimates below are extrapolations from these numbers:

| Job | 4-vCPU CPU (estimate) | Colab free CPU, 2 vCPU (estimate) |
|---|---|---|
| Level 1 optimisation (about 10–20 steps after xTB) | about 30–60 min | about 1–2 h |
| Level 1 TD-DFT, 30 roots, per functional | about 0.5–1 h | about 1–2 h |
| Level 2 (486 basis functions), per optimisation step | about 6–10 min | about 15–25 min |
| Level 2 TD-DFT, 30 roots, per functional | several hours | several hours |

Recommendations:

- **Memory:** PySCF honours `--max-memory` only loosely. With the default 10 000 MB, a Level 2 B3LYP/6-31G(d)/PCM job reached 13.6 GB resident and was killed on a 16 GB machine. Use `--max-memory 4000` there, which ran at about 4–5 GB.

- **CPU:** run Level 1 first. Run Level 2 as a separate session, resuming with `--start-xyz results/level2/opt_B3LYP_last.xyz` whenever a session ends. Save `--outdir` to Google Drive.
- **Speed:** `--tda` roughly halves the TD-DFT cost. Oscillator strengths are slightly less reliable with TDA, so use full TD-DFT for final numbers.
- **GPU:** `--gpu` (gpu4pyscf, any Colab GPU runtime) is expected to be about an order of magnitude faster. The GPU path follows gpu4pyscf's documented construction (`density_fit → to_gpu → PCM`) but was **not executed** in development, because no GPU was available. Check the first run's SCF energy against a CPU single point (`--skip-opt`); they should agree to about 1e-6 Eh.

## Tests (`tests/test_phase1.py`)

These checks run without a production calculation:
- broadening conserves ∫ε dE = 2.8707×10⁴ Σf;
- ε → MAC units (Napierian/decadic) and nm ↔ eV;
- both molecules' formulas, charges and multiplicities;
- a tiny real PySCF RKS + TD-DFT run.

## Scientific caveats to keep in view

- **Regiochemistry of the COOH (Level 2)** is set to C8 of the benzo[7]annulene skeleton. If your assignment differs, edit `LEVEL2_SMILES` in `molecules.py` and rerun with `--rebuild`.
- **Protonation state:** the carboxylic acid (pKa ≈ 3–4) is likely deprotonated at vacuolar pH. Both models are run neutral, as specified. A charge −1 carboxylate run is a recommended sensitivity test: set `charge=-1` in `MOLECULES['level2']`.
- **Tautomers and conformers:** benzotropolones have tropolone OH/C=O tautomers and OH rotamers, and the glycoside is flexible. Spectra should ideally be Boltzmann-averaged over the lowest few DFT conformers. The single MMFF-lowest conformer used here is a starting point.
- **Functional:** B3LYP tends to red-shift and over-stabilise charge-transfer states. CAM-B3LYP usually blue-shifts π→π* bands of polyphenols by about 0.2–0.4 eV. Compare both with measured HPLC-DAD/UV-vis spectra before trusting absolute λmax.
- **Vertical TD-DFT gives no vibronic structure.** The 0.3 eV Gaussian width is phenomenological, so test FWHM 0.25–0.40 eV with `spectra.build_spectrum`, which re-broadens saved results in seconds.
- `--freq` runs an analytic PCM Hessian on CPU to confirm a true minimum. It is expensive, so run it on Level 1 at least once.
