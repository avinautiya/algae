"""Fill the DarkZone_2026_data manifest entry from the inspected Zenodo archive (10.5281/zenodo.18826013)."""
import hashlib, json, os, sys
import h5py, numpy as np, pandas as pd
AR = os.environ["DARKZONE_ARCHIVE"]; ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
nc = os.path.join(AR, "data/spectra/HCRF_ASDFieldspec4_spectra_Chevrollier_et_al_2026.nc")
f = h5py.File(nc, "r"); tags = [t.decode() if isinstance(t, bytes) else t for t in f["measurement_tag"][:]]
dates = pd.Series([pd.to_datetime(t[:6], format="%d%m%y").date().isoformat() for t in tags])
repo = pd.read_csv(os.path.join(ROOT, "data/empirical/chevrollier2023_hcrf.csv"))
wl = f["wavelength"][:]; R = f["reflectance"][:]
rtags = [c for c in repo.columns if c != "wavelength_nm"]
match = []
for c in rtags:
    y = repo[c].to_numpy(float); w0 = repo["wavelength_nm"].to_numpy(float)
    best = min(((np.nanmax(np.abs(y[(w0 >= 400) & (w0 <= 1300)] - np.interp(w0, wl, R[:, j])[(w0 >= 400) & (w0 <= 1300)])), t) for j, t in enumerate(tags)))
    match.append((c, best[1], float(best[0])))
sha = hashlib.sha256(open(nc, "rb").read()).hexdigest()
entry = dict(
    id="DarkZone_2026_data", title="Dark Zone paper data/software archive (Chevrollier et al. 2026)", source="Zenodo 10.5281/zenodo.18826013 (CC BY 4.0), data_software.zip md5 0955d9c3ddf93ff9b69229d4ac28aa17",
    files=[dict(role="spectra", path="<archive>/data/spectra/HCRF_ASDFieldspec4_spectra_Chevrollier_et_al_2026.nc", sha256=sha),
           dict(role="emulator", path="<archive>/data/emulator/mlp_icelaps_101025.h5"), dict(role="code", path="<archive>/src/mcmc_pipeline.py")],
    spectral_quantity="HCRF, nadir, ASD FieldSpec4, fore-optic 10 degree FOV (directional; NOT albedo)", wavelengths_nm=[int(wl[0]), int(wl[-1])],
    n_spectra=len(tags), n_per_date=dates.value_counts().sort_index().to_dict(), sites="S21 (61.10 N, 2021) and S22A (69.43 N, 2022)",
    measured_components="none in the archive. Paper: dust MAC n=10 samples, mineralogy n=17 (Table S1), cryoconite spectrum derived from 2 reflectance measurements at S21; NO per-spectrum concentrations, counts or pigments",
    sky_conditions="recorded per spectrum in the paper's supplement tables S2/S4 (clear to fully cloudy) for 21 tabulated spectra only; not in the archive file",
    overlap_with_repo=dict(repo_dataset="SGreenland_2021_hcrf_counts", repo_spectra=len(rtags), identical_in_archive=int(sum(m[2] < 1e-4 for m in match)), max_abs_diff=max(m[2] for m in match)),
    independent_glaciers=2, role_in_this_study="reference for the competitor reproduction and for reconstruction-only checks; NOT a concentration-validation set; the 20 spectra already in the repo are reused data",
    ancestry=["the archive's algae and cryoconite optics were measured/derived at S21, the same location as 2021 spectra", "the 20 repo transfer-check spectra are a subset"],
    leakage_notes="S21 is both an optics source and a retrieval target in the competitor's workflow",
)
p = os.path.join(ROOT, "records/inversion_m1/data_manifest.json"); m = json.load(open(p))
m["datasets"] = [entry if d["id"] == "DarkZone_2026_data" else d for d in m["datasets"]]
json.dump(m, open(p, "w"), indent=1, default=str)
print(json.dumps({k: entry[k] for k in ("n_spectra", "n_per_date", "overlap_with_repo")}, indent=1, default=str))
