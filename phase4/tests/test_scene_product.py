"""Scene product: screening classes, physics interpretation on a small synthetic emulator, export
round trip (COG, NetCDF, sidecar checksums), immutability. No network, no real scene."""
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(HERE, "..", "..", "phase2"))

import scene_product as SP  # noqa: E402
from emulator import Emulator  # noqa: E402

C = np.array([0.62, 0.60, 0.56, 0.46])
A = np.array([0.9, 0.75, 0.35, 0.03])          # algae absorb blue-green >> red > NIR
D = np.array([0.30, 0.28, 0.26, 0.22])         # dust: flat
SW = 700.0


def fake_emulator(model="ours"):
    axes = dict(log_b=np.arange(1.0, 6.01, 0.25), f_n=np.linspace(0, 1, 3) if model == "ours" else np.array([0.5]),
                r_um=np.geomspace(300, 20000, 6), dust_ppb=np.array([3e4, 1.5e5, 6e5, 1.2e6]))
    lb, fn, r, d = np.meshgrid(*axes.values(), indexing="ij")
    rf = 1.0 - 0.15 * (np.log(r) - np.log(300)) / (np.log(20000) - np.log(300))
    B = 10 ** lb
    clean_d = C * rf[..., None] * np.exp(-D * d[..., None] / 1e6)
    bands = clean_d * np.exp(-A * B[..., None] / 2e4)
    w = np.array([0.3, 0.3, 0.25, 0.15])
    bba = bands @ w
    rf_alg = SW * (clean_d @ w - bba)
    return Emulator(axes, dict(bands=bands, bba=bba, rf_algae=rf_alg, rf_total=rf_alg, pigment_ug_l=B * 1e-3),
                    dict(sw_down=SW, model=model))


@pytest.fixture
def phys():
    cal = dict(sigma=0.01, radius_median_um=2000.0, radius_ln_sd=0.9, tau_dex=0.7)
    return {SP.PRIMARY: dict(em=fake_emulator(), cal=cal, path="x", sha256="x"),
            SP.ALTERNATIVE: dict(em=fake_emulator("tierA"), cal=dict(cal, tau_dex=1.0), path="y", sha256="y")}


def test_screening_classes():
    R = np.array([[0.5, 0.48, 0.45, 0.35],        # ice, SCL 11
                  [0.5, 0.48, 0.45, 0.35],        # cloud
                  [0.30, 0.15, 0.08, 0.04],       # pond (NDWI_ice > 0.25)
                  [0.95, 0.93, 0.9, 0.8],         # snow
                  [np.nan, 0.1, 0.1, 0.1],        # nodata
                  [0.3, 0.3, 0.3, 0.3]])          # SCL 5 bare soil
    scl = np.array([11, 9, 11, 11, 11, 5])
    cls, flags, _ = SP.screen(R, scl)
    assert list(cls) == [255, 1, 2, 3, 0, 9]
    assert flags[2] & (1 << 2) and flags[3] & (1 << 3) and flags[1] & (1 << 1)


def test_interpretation_separates_algae_clean_and_excluded(phys):
    em = fake_emulator()
    # pixel 0: heavy algae, low dust; pixel 1: clean ice, low dust; pixel 2: dark SCL 7, heavy algae; 3: cloud
    def node(lb, d, r=2):
        i = int(np.argmin(abs(em.axes["log_b"] - lb)))
        return em.data["bands"][i, 1, r, d]
    R = np.vstack([node(4.75, 0), node(1.0, 0), node(4.75, 0), node(4.0, 0)])
    out, diag = SP.interpret(R, np.array([11, 11, 7, 8]), phys)
    cls = out["interpretation_class"]
    assert cls[0] == 5 and cls[1] == 4 and cls[2] == 5 and cls[3] == 1
    assert out["quality_flags"][2] & (1 << 4)                     # dark SCL kept but flagged
    assert out["dalpha_algae_mean"][0] > SP.DALPHA_MIN and abs(out["dalpha_algae_mean"][1]) < SP.DALPHA_MIN
    assert np.isnan(out["log10_cells_mean"][3])
    assert out["log_bayes_factor_algae"][0] > SP.LOG_BF_MIN
    assert "f_n_mean" not in out                                  # no species fractions


def test_export_roundtrip_and_immutability(phys, tmp_path):
    from rasterio.crs import CRS
    from rasterio.transform import from_origin
    em = fake_emulator()
    H, W = 4, 5
    rng = np.random.default_rng(0)
    idx = rng.integers(0, em.shape[0], H * W)
    R = em.data["bands"][idx, 1, 2, 0] + rng.normal(0, 0.005, (H * W, 4))
    scl = np.full(H * W, 11)
    scl[:3] = 9
    R[5] = np.nan
    vals, diag = SP.interpret(R, scl, phys)
    maps = {k: v.reshape(H, W) for k, v in vals.items()}
    tr, crs = from_origin(500000.0, 7450000.0, 20.0, 20.0), CRS.from_epsg(32622)
    prov = dict(scene_id="TEST", validation_status="test fixture")
    paths = SP.export(maps, tr, crs, str(tmp_path), "TEST", prov)
    checks = SP.validate_outputs(paths, maps, tr, crs)
    assert all(checks.values())
    with pytest.raises(FileExistsError):
        SP.export(maps, tr, crs, str(tmp_path), "TEST", prov)
    # a tampered output fails the checksum check
    with open(paths["netcdf"], "ab") as fh:
        fh.write(b"0")
    with pytest.raises(AssertionError):
        SP.validate_outputs(paths, maps, tr, crs)
    rep = SP.report(maps, diag, prov, 20.0)
    assert "not independently validated" in rep


def test_out_of_domain_sza_refused():
    with pytest.raises(ValueError, match="out of domain"):
        SP.load_physics(70.2)
