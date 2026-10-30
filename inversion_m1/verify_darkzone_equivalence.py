"""Numerical equivalence check: NumPy re-implementation (reproduce_darkzone.py) versus the ORIGINAL archived code path.

Runs in an environment with TensorFlow 2.21 + TensorFlow Probability 0.25 (the paper reports TFP 0.25). Compares, on identical
inputs: (A) scaler constants read from the raw pickle bytes (the pickle is NOT executed) against those used in the NumPy
code, with the scaler's own internal identities; (B) emulator outputs (Keras vs NumPy); (C) the full unconstrained target
log-density, with the TFP block copied verbatim from src/mcmc_pipeline.py l.49-158; (D) sensitivity of the observation
slice to the (inferred, unverified) 2.5 nm band-centre offset.
Scaler stand-in: X*scale_+min_ with scale_/min_ read from the file (sklearn MinMaxScaler.transform semantics).
"""
import json, os, struct, sys
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import reproduce_darkzone as RD
AR = os.environ["DARKZONE_ARCHIVE"]
out = {}

# ---- (A) constants from raw bytes -----------------------------------------------------------------------------
b = open(os.path.join(AR, "data/emulator/minmax_scaler_icelaps_101025.save"), "rb").read()
blocks = {}
for i in range(len(b) - 56):
    v = np.array(struct.unpack("<7d", b[i:i + 56]))
    if np.all(np.isfinite(v)) and np.all(np.abs(v) < 1e12) and np.all((v == 0) | (np.abs(v) > 1e-9)) and np.any(v != 0):
        blocks[i] = v
offs = sorted(blocks)
cand = {o: blocks[o] for o in offs}
def find(pred): return [o for o in cand if pred(cand[o])]
o_scale = find(lambda v: np.allclose(v, 1 / RD.RANGE, rtol=1e-9) )
o_dmin = find(lambda v: np.allclose(v, RD.DATA_MIN)); o_dmax = find(lambda v: np.allclose(v, RD.DATA_MAX)); o_rng = find(lambda v: np.allclose(v, RD.RANGE))
o_min_ = find(lambda v: np.allclose(v, -RD.DATA_MIN / RD.RANGE, rtol=1e-9, atol=1e-12))
out["A_scaler"] = dict(scale__matches_1_over_range=bool(o_scale), data_min__matches=bool(o_dmin), data_max__matches=bool(o_dmax),
                       data_range__matches=bool(o_rng), min__matches_minus_dmin_times_scale=bool(o_min_),
                       offsets=dict(scale_=o_scale, min_=o_min_, data_min_=o_dmin, data_max_=o_dmax, data_range_=o_rng),
                       method="raw-byte scan of float64 blocks; arrays identified by value, so a block is only 'found' if it equals the NumPy-code constant")

# ---- (B) emulator outputs -------------------------------------------------------------------------------------
import tensorflow as tf

model = tf.keras.models.load_model(os.path.join(AR, "data/emulator/mlp_icelaps_101025.h5"), compile=False)
layers = RD.load_mlp(os.path.join(AR, "data/emulator/mlp_icelaps_101025.h5"))
rng = np.random.default_rng(1)
X = rng.uniform(0, 1, size=(2000, 7))
y_tf = model(X.astype("float32")).numpy().astype("float64"); y_np = RD.mlp(layers, X)
d = np.abs(y_tf - y_np)
out["B_emulator"] = dict(n_inputs=2000, output_shape=list(y_tf.shape), max_abs_diff=float(d.max()), mean_abs_diff=float(d.mean()),
                         max_abs_diff_in_fit_slice_10_181=float(d[:, 10:181].max()), note="float32 TF vs float64 NumPy; unit-cube inputs, uniform random")
# ---- (C) target log density, original TFP block ---------------------------------------------------------------
import tensorflow_probability as tfp
tfd, tfb = tfp.distributions, tfp.bijectors
class Scaler:                                         # sklearn MinMaxScaler.transform semantics with the file's constants
    def transform(self, X): return np.asarray(X) / RD.RANGE + (-RD.DATA_MIN / RD.RANGE)
class Emu: pass
my_emulator = Emu(); my_emulator.scaler = Scaler(); my_emulator.emulator = model
n_params = 7; sigmoid_bijector = tfb.Sigmoid(); inverse_sigmoid_bijector = tfb.Invert(sigmoid_bijector)
beta_prior_dz1 = tfd.TransformedDistribution(distribution=tfd.Beta(2.0, 3.5), bijector=tfb.Scale(0.4666))
unconstrained_beta_prior_dz1 = tfd.TransformedDistribution(beta_prior_dz1, inverse_sigmoid_bijector)
dz2_beta = tfd.TransformedDistribution(distribution=tfd.Beta(2.0, 3.5), bijector=tfb.Scale(0.7999992))
prior_uniform_dz2 = tfd.TransformedDistribution(dz2_beta, inverse_sigmoid_bijector)
dd = tfd.Beta(0.95, 3.0)
beta_prior_laps = [tfd.TransformedDistribution(distribution=dd, bijector=tfb.Scale(tf.constant(my_emulator.scaler.transform(np.tile(scale, (7, 1)).T)[:, idx][0], dtype=tf.float32)))
                   for scale, idx in zip([5e6, 1.5e6, 5e5], [-3, -2, -1])]
unconstrained_beta_priors_laps = [tfd.TransformedDistribution(beta_prior_laps[i], inverse_sigmoid_bijector) for i in range(3)]
def log_prior_fun(z):
    z = tf.reshape(z, [-1, n_params]); dz1_z = z[..., 2]; dz2_z = z[..., 3]; laps_z = z[..., 4:]
    return (unconstrained_beta_prior_dz1.log_prob(dz1_z) + prior_uniform_dz2.log_prob(dz2_z)
            + tf.add_n([p.log_prob(laps_z[..., i]) for i, p in enumerate(unconstrained_beta_priors_laps)]))
def poly_tensor_vec(x, coeffs):
    powers = tf.range(len(coeffs), dtype=tf.float32)
    return tf.reduce_sum(coeffs * tf.pow(tf.expand_dims(x, -1), powers), axis=-1)
POLY_COEFFS = [0.02812912, -1.39078195, 5.12318554, -5.49135067, 2.72952866]
def rho2_min_fn(rho1):
    return tf.where(rho1 < 0.479, 0.07410701, poly_tensor_vec(rho1, POLY_COEFFS))
def rho_forward(z1, z2):
    rho1 = sigmoid_bijector.forward(z1); rho2_min = rho2_min_fn(rho1)
    return tf.stack([rho1, rho2_min + (1.0 - rho2_min) * sigmoid_bijector.forward(z2)], axis=-1)
def log_det_jacobian(z1, z2):
    rho1 = sigmoid_bijector.forward(z1); drho1_dz1 = rho1 * (1.0 - rho1); rho2_min = rho2_min_fn(rho1)
    sz2 = sigmoid_bijector.forward(z2); return tf.math.log(drho1_dz1 * (1.0 - rho2_min) * sz2 * (1.0 - sz2) + 1e-12)
import h5py
f = h5py.File(os.path.join(AR, "data/spectra/HCRF_ASDFieldspec4_spectra_Chevrollier_et_al_2026.nc"), "r")
obs = f["reflectance"][:][0:851:5, :]
sigma = 0.003
def target_original(z, y_obs):
    z = tf.reshape(z, [-1, n_params]); rhos = rho_forward(z[..., 0], z[..., 1])
    dz1 = sigmoid_bijector.forward(z[..., 2]); dz2 = sigmoid_bijector.forward(z[..., 3]); laps = sigmoid_bijector.forward(z[..., 4:])
    theta = tf.concat([rhos, dz1[..., None], dz2[..., None], laps], axis=-1)
    y_pred = my_emulator.emulator(theta)[:, 10:10 + 171]; resid = (y_pred - y_obs) / sigma; nu = 3.0
    ll = -0.5 * (nu + 1.0) * tf.reduce_sum(tf.math.log1p((resid ** 2) / nu), axis=-1)
    return ll + log_prior_fun(z) + log_det_jacobian(z[..., 0], z[..., 1])
diffs, rel, n_cmp = [], [], 0
for trial in range(400):
    j = int(rng.integers(0, obs.shape[1]))
    y = obs[:, j]
    if not np.all(np.isfinite(y)): continue
    z = np.concatenate([rng.normal(0, 1.5, 4), rng.normal(-1.5, 1.5, 3)])
    a = float(target_original(tf.constant(z[None], tf.float32), tf.constant(y[None], tf.float32)).numpy()[0])
    bnp, _ = RD.log_target(z, layers, y)
    if np.isfinite(a) and np.isfinite(bnp):
        diffs.append(abs(a - bnp)); rel.append(abs(a - bnp) / max(abs(a), 1.0)); n_cmp += 1
    else:
        out.setdefault("C_nonfinite_pairs", []).append([float(a), float(bnp)])
out["C_target_density"] = dict(n_compared=n_cmp, max_abs_diff=float(np.max(diffs)), median_abs_diff=float(np.median(diffs)), max_rel_diff=float(np.max(rel)),
                               note="original TFP block (verbatim) vs NumPy; float32 vs float64; z drawn from broad normals incl. regions far from the data")
if "C_nonfinite_pairs" in out: out["C_nonfinite_pairs"] = out["C_nonfinite_pairs"][:10]
# ---- (D) 2.5 nm observation-offset sensitivity ----------------------------------------------------------------
wl = f["wavelength"][:]; R = f["reflectance"][:]
full = R[(wl >= 350) & (wl <= 1210), :]; w = wl[(wl >= 350) & (wl <= 1210)]
pt = np.stack([np.interp(np.arange(350, 1201, 5), w, full[:, j]) for j in range(R.shape[1])], 1)
sh = np.stack([np.interp(np.arange(352.5, 1202.6, 5), w, full[:, j]) for j in range(R.shape[1])], 1)
dd_ = (sh - pt)[:, np.all(np.isfinite(pt), 0)]
out["D_band_centre_offset"] = dict(hypothesis="emulator outputs are band centres 302.5+5k (440 outputs from [300,2500,5] band mode, session.py l.126-137), observations sampled at 350+5k",
                                   verified=False, rms_diff=float(np.sqrt((dd_ ** 2).mean())), p95_abs_diff=float(np.percentile(np.abs(dd_), 95)), max_abs_diff=float(np.abs(dd_).max()),
                                   likelihood_sigma=sigma)
json.dump(out, open(os.path.join(HERE, "..", "records", "inversion_m1", "darkzone_repro", "equivalence_check.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
