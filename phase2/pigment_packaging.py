"""
Intracellular pigment packaging (Duysens "flattening" / self-shading) corrections.

Physics
-------
A pigment dissolved homogeneously in a medium absorbs with the solution
(Napierian) absorption coefficient

    a_sol(lambda) = MAC(lambda) * c          [m^-1]      (MAC in m^2 kg^-1, c in kg m^-3)

When the same mass of pigment is packed into discrete cells at intracellular
concentration c_i, each cell is optically thick at strongly absorbed
wavelengths: molecules at the back of the cell sit in the "shadow" of those
at the front. The absorption cross-section of one cell is then (geometric /
anomalous-diffraction optics, refractive index of the cell close to that of
its surroundings; Duysens 1956, BBA 19:1; Morel & Bricaud 1981, Deep-Sea Res. 28:1375):

    sigma_abs = A_proj * < 1 - exp(-a_i * l) >                         (1)

where a_i = MAC * c_i is the intracellular absorption coefficient, l is the
geometric chord length of a ray through the cell, <.> averages over all rays
hitting the randomly oriented cell, and A_proj = S/4 is its mean projected
area (Cauchy). Because <l> = 4V/S = V/A_proj (Cauchy's mean-chord theorem),
the unpackaged ("solution") cross-section of the same pigment is

    sigma_sol = a_i * V = a_i * <l> * A_proj                           (2)

and the packaging (flattening) factor is

    Q*(lambda) = sigma_abs / sigma_sol = <1 - exp(-a_i l)> / (a_i <l>),   0 < Q* <= 1   (3)

so that the in vivo mass absorption coefficient is

    MAC_vivo(lambda) = Q*(lambda) * MAC(lambda).                        (4)

Sphere of diameter d (closed form, Morel & Bricaud 1981, rho' = a_i d):

    Q_a(rho') = 1 + 2 exp(-rho')/rho' + 2 (exp(-rho') - 1)/rho'^2
    Q*(rho')  = (3/2) Q_a(rho') / rho'                                  (5)

Limits: Q* -> 1 - 3 rho'/8 (rho' -> 0), Q* -> 3/(2 rho') (rho' -> inf,
black sphere: sigma_abs -> pi d^2/4).

Ancylonema cells are short cylinders (rod-shaped, ~10 um wide, 20-40 um long),
so for non-spherical cells we evaluate Eq. (3) exactly with a Monte-Carlo
sample of mu-randomness chords (isotropic, uniform random lines), which
reproduces Eq. (5) for spheres to <0.3 % (see tests/test_packaging.py).

Assumptions (state them in the paper): ray optics without refraction (relative
refractive index of cell vs. meltwater/ice ~1.02-1.08), pigment homogeneously
distributed in the pigment-bearing compartment, monodisperse cells.
Set `vacuole_fraction` < 1 to confine the pigment to a vacuole occupying that
fraction of the cell volume (treated as an equivalent-volume compartment of
the same shape) - this increases c_i and hence the packaging effect.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


# --------------------------------------------------------------------------- #
# Analytic sphere (Duysens / Morel & Bricaud)                                  #
# --------------------------------------------------------------------------- #
def q_star_sphere(rho):
    """Packaging factor Q*(rho') for a homogeneous sphere, rho' = a_i * d.

    Uses the Taylor series 1 - 3x/8 + x^2/10 - x^3/48 + x^4/280 for x < 1e-2
    to avoid catastrophic cancellation in Eq. (5).
    """
    x = np.asarray(rho, dtype=float)
    out = np.ones_like(x)
    small = x < 1e-2
    xs = x[small]
    out[small] = 1.0 - 3.0 * xs / 8.0 + xs**2 / 10.0 - xs**3 / 48.0 + xs**4 / 280.0
    xl = x[~small]
    ex = np.exp(-xl)
    q_a = 1.0 + 2.0 * ex / xl + 2.0 * (ex - 1.0) / xl**2
    out[~small] = 1.5 * q_a / xl
    return out


# --------------------------------------------------------------------------- #
# General convex cells: chord-length distribution                              #
# --------------------------------------------------------------------------- #
def _orthonormal_basis(u):
    """Two unit vectors perpendicular to each row of u (n,3)."""
    helper = np.where(np.abs(u[:, [2]]) < 0.9, np.array([[0.0, 0.0, 1.0]]), np.array([[1.0, 0.0, 0.0]]))
    e1 = np.cross(u, helper)
    e1 /= np.linalg.norm(e1, axis=1, keepdims=True)
    e2 = np.cross(u, e1)
    return e1, e2


def _random_lines(n, bound_radius, rng):
    """mu-random lines through a bounding sphere: isotropic direction u and an
    impact point uniform over the disc of radius `bound_radius` normal to u."""
    u = rng.normal(size=(n, 3))
    u /= np.linalg.norm(u, axis=1, keepdims=True)
    e1, e2 = _orthonormal_basis(u)
    r = bound_radius * np.sqrt(rng.random(n))
    phi = 2.0 * np.pi * rng.random(n)
    p = (r * np.cos(phi))[:, None] * e1 + (r * np.sin(phi))[:, None] * e2
    return p, u


def chords_sphere(radius, n=200_000, seed=1):
    rng = np.random.default_rng(seed)
    b = radius * np.sqrt(rng.random(n))        # impact parameter uniform over the disc
    return 2.0 * np.sqrt(radius**2 - b**2)


def chords_cylinder(radius, length, n=200_000, seed=1, max_iter=50):
    """Chord lengths of mu-random lines through a right circular cylinder.

    Cylinder axis = z, |z| <= length/2, x^2 + y^2 <= radius^2. Lines that miss
    are rejected; hits are therefore uniform over the projected area for each
    orientation, and orientations are weighted by projected area as required
    for mu-randomness.
    """
    rng = np.random.default_rng(seed)
    h = 0.5 * length
    bound = np.sqrt(radius**2 + h**2)
    chords = []
    n_have = 0
    for _ in range(max_iter):
        p, u = _random_lines(2 * (n - n_have) + 1000, bound, rng)
        # infinite cylinder: (px + t ux)^2 + (py + t uy)^2 <= R^2
        A = u[:, 0] ** 2 + u[:, 1] ** 2
        B = 2.0 * (p[:, 0] * u[:, 0] + p[:, 1] * u[:, 1])
        C = p[:, 0] ** 2 + p[:, 1] ** 2 - radius**2
        disc = B**2 - 4.0 * A * C
        with np.errstate(divide="ignore", invalid="ignore"):
            sq = np.sqrt(np.maximum(disc, 0.0))
            t1c = np.where(A > 1e-14, (-B - sq) / (2.0 * A), -np.inf)
            t2c = np.where(A > 1e-14, (-B + sq) / (2.0 * A), np.inf)
        side_hit = np.where(A > 1e-14, disc > 0, C <= 0)
        # end-cap slab |pz + t uz| <= h
        with np.errstate(divide="ignore", invalid="ignore"):
            ta = (-h - p[:, 2]) / u[:, 2]
            tb = (h - p[:, 2]) / u[:, 2]
        par = np.abs(u[:, 2]) < 1e-14
        t1s = np.where(par, -np.inf, np.minimum(ta, tb))
        t2s = np.where(par, np.inf, np.maximum(ta, tb))
        slab_ok = np.where(par, np.abs(p[:, 2]) <= h, True)
        l = np.minimum(t2c, t2s) - np.maximum(t1c, t1s)
        hit = side_hit & slab_ok & (l > 0)
        chords.append(l[hit])
        n_have += int(hit.sum())
        if n_have >= n:
            break
    return np.concatenate(chords)[:n]


@dataclass
class CellGeometry:
    """Pigment-bearing cell (or vacuole) geometry in micrometres.

    shape = 'cylinder' (Ancylonema, default) or 'sphere'.
    For a cylinder, `radius` is the cross-section radius and `length` the axial length.
    """
    shape: str = "cylinder"
    radius: float = 5.0          # um
    length: float = 20.0         # um (ignored for spheres)
    n_chords: int = 200_000
    seed: int = 1
    _chords: np.ndarray | None = field(default=None, repr=False)

    @property
    def volume_um3(self) -> float:
        if self.shape == "sphere":
            return 4.0 / 3.0 * np.pi * self.radius**3
        return np.pi * self.radius**2 * self.length

    @property
    def surface_um2(self) -> float:
        if self.shape == "sphere":
            return 4.0 * np.pi * self.radius**2
        return 2.0 * np.pi * self.radius**2 + 2.0 * np.pi * self.radius * self.length

    @property
    def projected_area_um2(self) -> float:
        """Orientation-averaged projected area, S/4 (Cauchy's theorem, convex body)."""
        return self.surface_um2 / 4.0

    @property
    def mean_chord_um(self) -> float:
        return 4.0 * self.volume_um3 / self.surface_um2

    @property
    def equivalent_sphere_diameter_um(self) -> float:
        return (6.0 * self.volume_um3 / np.pi) ** (1.0 / 3.0)

    def chords_um(self) -> np.ndarray:
        if self._chords is None:
            if self.shape == "sphere":
                self._chords = chords_sphere(self.radius, self.n_chords, self.seed)
            elif self.shape == "cylinder":
                self._chords = chords_cylinder(self.radius, self.length, self.n_chords, self.seed)
            else:
                raise ValueError(f"unknown shape {self.shape!r}")
        return self._chords

    def scaled(self, volume_fraction: float) -> "CellGeometry":
        """Same shape with volume scaled by `volume_fraction` (e.g. a vacuole)."""
        s = volume_fraction ** (1.0 / 3.0)
        return CellGeometry(self.shape, self.radius * s, self.length * s, self.n_chords, self.seed)

    def label(self, plain: bool = False) -> str:
        """Figure label (mathtext) or plain-text label; cylinders as diameter x length."""
        if plain:
            if self.shape == "sphere":
                return f"sphere d={2 * self.radius:g} um"
            return f"cylinder d={2 * self.radius:g} um x L={self.length:g} um"
        if self.shape == "sphere":
            return rf"sphere $d={2 * self.radius:g}\,\mu$m"
        return rf"cylinder $d\times L={2 * self.radius:g}\times{self.length:g}\,\mu$m"


def q_star(a_internal_per_m, geom: CellGeometry, method: str = "auto"):
    """Packaging factor Q*(lambda) for intracellular absorption coefficient a_i [m^-1].

    method = 'analytic' (spheres only), 'chord' (Monte-Carlo chords, any shape) or
    'auto' (analytic for spheres, chords otherwise).
    """
    a = np.asarray(a_internal_per_m, dtype=float)
    if method == "auto":
        method = "analytic" if geom.shape == "sphere" else "chord"
    if method == "analytic":
        if geom.shape != "sphere":
            raise ValueError("analytic Q* only exists for spheres; use method='chord'")
        return q_star_sphere(a * 2.0 * geom.radius * 1e-6)

    l = geom.chords_um() * 1e-6                           # m
    # Normalise by the SAMPLE mean chord so that Q* -> 1 exactly in the dilute
    # limit; it differs from the exact 4V/S by O(N^-1/2) ~ 0.1 % (checked in tests).
    lbar = l.mean()
    out = np.ones_like(a)
    flat = a.ravel()
    res = out.ravel()
    pos = flat > 0
    # <1 - exp(-a l)> / (a <l>), evaluated in chunks to bound memory
    for start in range(0, pos.sum(), 64):
        idx = np.flatnonzero(pos)[start:start + 64]
        x = flat[idx, None] * l[None, :]
        res[idx] = (-np.expm1(-x)).mean(axis=1) / (flat[idx] * lbar)
    return res.reshape(a.shape)


# --------------------------------------------------------------------------- #
# MAC -> in vivo MAC                                                           #
# --------------------------------------------------------------------------- #
def mac_vivo(mac_m2_per_kg, c_internal_kg_m3: float, geom: CellGeometry,
             vacuole_fraction: float = 1.0, background_abs_per_m=None):
    """In vivo (packaging-corrected) MAC.

    Parameters
    ----------
    mac_m2_per_kg : array, Napierian solution MAC of the pigment (Phase 1 output).
    c_internal_kg_m3 : pigment mass per unit CELL volume [kg m^-3]
                       (= pigment mass per cell / cell volume).
    geom : CellGeometry of the cell.
    vacuole_fraction : fraction of the cell volume holding the pigment; the pigment
                       concentration inside that compartment is c_internal / vacuole_fraction.
    background_abs_per_m : optional array, other absorption inside the compartment
                       (e.g. cell water). It shares the self-shading but is not
                       counted as pigment absorption.

    Returns
    -------
    mac_vivo [m^2 kg^-1], q_star (pigment packaging factor), a_internal [m^-1]
    """
    mac = np.asarray(mac_m2_per_kg, dtype=float)
    g = geom if vacuole_fraction >= 1.0 else geom.scaled(vacuole_fraction)
    a_pig = mac * (c_internal_kg_m3 / vacuole_fraction)
    a_bg = 0.0 if background_abs_per_m is None else np.asarray(background_abs_per_m, dtype=float)
    a_tot = a_pig + a_bg
    q = q_star(a_tot, g)
    # The total absorbed fraction is shared between pigment and background in
    # proportion to their absorption coefficients (same photon paths).
    return mac * q, q, a_pig


def mac_vivo_grid(mac_m2_per_kg, radii_um, c_internal_kg_m3, shape: str = "cylinder",
                  aspect_ratio: float = 2.0, vacuole_fraction: float = 1.0):
    """Packaging-corrected MAC on a (radius x concentration x wavelength) grid.

    radii_um        : iterable of characteristic cell sizes; for cylinders this is the
                      cross-section radius and length = aspect_ratio * 2 * radius.
    c_internal_kg_m3: iterable of intracellular pigment concentrations.
    Returns array of shape (n_radius, n_conc, n_wavelength) and the matching Q* array.
    """
    mac = np.asarray(mac_m2_per_kg, dtype=float)
    radii = np.atleast_1d(radii_um).astype(float)
    concs = np.atleast_1d(c_internal_kg_m3).astype(float)
    out = np.empty((radii.size, concs.size, mac.size))
    qs = np.empty_like(out)
    for i, r in enumerate(radii):
        geom = CellGeometry(shape, r, aspect_ratio * 2.0 * r)
        for j, c in enumerate(concs):
            out[i, j], qs[i, j], _ = mac_vivo(mac, c, geom, vacuole_fraction)
    return out, qs
