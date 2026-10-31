"""Conservative transverse finite volumes and an independent plug-flow solution.

rho cp u(y) T_x = k T_yy, T(x,0)=0, T_y(x,H)=0.
The axial coordinate is the evolution variable; this is NOT transient CFD.
Fixed geometry/velocity, constant properties, no axial diffusion or sources.
"""
from dataclasses import dataclass
import numpy as np
from scipy.linalg import eigh_tridiagonal


def positive(value, name, allow_zero=False):
    value = float(value)
    if not np.isfinite(value) or value < 0 or (value == 0 and not allow_zero):
        raise ValueError("Invalid " + name)
    return value


@dataclass(frozen=True)
class Channel:
    depth_m: float = 0.002
    speed_m_s: float = 0.025
    width_m: float = 0.1
    rho_water: float = 1000.0
    cp: float = 4217.0
    conductivity: float = 0.561
    rho_ice: float = 917.0
    latent_heat: float = 334000.0
    nu_m2_s: float = 1.0e-6
    velocity: str = "shear_free"

    def __post_init__(self):
        for name, value in vars(self).items():
            if name != "velocity":
                positive(value, name)
        if self.velocity not in ("plug", "shear_free"):
            raise ValueError("velocity must be plug or shear_free")

    @property
    def diffusivity(self):
        return self.conductivity / (self.rho_water * self.cp)

    @property
    def capacity_W_K(self):
        return self.rho_water * self.cp * self.speed_m_s * self.depth_m * self.width_m

    def dimensionless(self, temperature_cap_K):
        h, u = self.depth_m, self.speed_m_s
        return dict(Re_h=u*h/self.nu_m2_s, Pe_h=u*h/self.diffusivity,
                    Pr=self.nu_m2_s/self.diffusivity, Fr=u/np.sqrt(9.81*h),
                    Ste_cap=self.cp*temperature_cap_K/self.latent_heat,
                    thermal_length_m=u*h*h/self.diffusivity)


class CrossStream:
    def __init__(self, channel=Channel(), cells=96, stretch=1.0):
        if type(cells) is not int or not 4 <= cells <= 1024:
            raise ValueError("cells must be an integer in [4,1024]")
        positive(stretch, "stretch")
        if not 1 <= stretch <= 2:
            raise ValueError("stretch must lie in [1,2]")
        self.channel, self.cells = channel, cells
        self.edges = channel.depth_m*np.linspace(0, 1, cells+1)**stretch
        self.dy = np.diff(self.edges)
        self.y = (self.edges[1:]+self.edges[:-1])/2
        a, b = self.edges[:-1]/channel.depth_m, self.edges[1:]/channel.depth_m
        if channel.velocity == "plug":
            self.u = np.full(cells, channel.speed_m_s)
        else:
            self.u = 1.5*channel.speed_m_s*(a+b-(a*a+a*b+b*b)/3)
        # Capacity per unit width; exact cell-averaged hydraulic velocity.
        self.mass = channel.rho_water*channel.cp*self.u*self.dy
        self.root_mass = np.sqrt(self.mass)
        interior = channel.conductivity/np.diff(self.y)
        self.bed_conductance = channel.conductivity/self.y[0]
        diagonal = np.zeros(cells)
        diagonal[0] += self.bed_conductance
        diagonal[:-1] += interior
        diagonal[1:] += interior
        self.diagonal_K = diagonal
        self.offdiagonal_K = -interior
        self.eigenvalues, self.eigenvectors = eigh_tridiagonal(
            diagonal/self.mass, -interior/(self.root_mass[:-1]*self.root_mass[1:]))
        if self.eigenvalues[0] <= 0:
            raise ArithmeticError("Positive thermal generator required")

    def propagate(self, profile_K, distance_m):
        distance_m = positive(distance_m, "distance", True)
        profile = np.asarray(profile_K, dtype=float)
        if profile.shape != (self.cells,) or not np.isfinite(profile).all():
            raise ValueError("One finite initial temperature per transverse cell required")
        v = self.eigenvectors
        return (v@(np.exp(-self.eigenvalues*distance_m)*(v.T@(self.root_mass*profile))))/self.root_mass

    def bulk_kernel(self, distance_m):
        distance_m = positive(distance_m, "distance", True)
        if distance_m == 0:
            return self.mass/self.mass.sum()
        v = self.eigenvectors
        return self.root_mass*(v@(np.exp(-self.eigenvalues*distance_m)*(v.T@self.root_mass)))/self.mass.sum()

    def bulk(self, profile_K, distance_m):
        return float(self.bulk_kernel(distance_m)@np.asarray(profile_K))

    def heat_kernel_W_K(self, start_m, stop_m):
        if not 0 <= start_m < stop_m:
            raise ValueError("Require 0 <= region start < stop")
        return self.channel.capacity_W_K*(self.bulk_kernel(start_m)-self.bulk_kernel(stop_m))

    def wall_flux_W_m2(self, profile_K, distance_m):
        return float(self.bed_conductance*self.propagate(profile_K, distance_m)[0])

    def lowering_factor(self, start_m, stop_m, horizon_s):
        positive(horizon_s, "horizon")
        if not 0 <= start_m < stop_m:
            raise ValueError("Invalid region")
        c = self.channel
        return horizon_s/(c.rho_ice*c.latent_heat*c.width_m*(stop_m-start_m))

    def surface_response(self, start_m, stop_m, distance_m):
        """Temperature response per W/m2 of constant top flux on [start,stop]."""
        if not 0 <= start_m < stop_m or distance_m < 0:
            raise ValueError("Invalid surface-source interval")
        end = min(stop_m,distance_m)
        if end <= start_m:
            return np.zeros(self.cells)
        lam, v = self.eigenvalues, self.eigenvectors
        integral = np.exp(-lam*(distance_m-end))*(-np.expm1(-lam*(end-start_m)))/lam
        source = np.zeros(self.cells)
        source[-1] = 1/self.root_mass[-1]
        return (v@(integral*(v.T@source)))/self.root_mass

    def surface_heat_coefficient(self, source_start_m, source_stop_m, start_m, stop_m):
        """Regional bed power per unit top flux; includes source in energy ledger."""
        ra = self.surface_response(source_start_m,source_stop_m,start_m)
        rb = self.surface_response(source_start_m,source_stop_m,stop_m)
        overlap = max(0,min(stop_m,source_stop_m)-max(start_m,source_start_m))
        return float(self.channel.width_m*(self.mass@(ra-rb)+overlap))


def plug_cell_bulk_kernel(channel, edges, distance_m, modes=1024):
    """Analytic sine-series kernel integrated over piecewise-constant inlet cells.

    Independent of the finite-volume eigensystem. For x>0 truncation converges
    exponentially. The numerical comparison must also check mode refinement.
    """
    if channel.velocity != "plug":
        raise ValueError("Analytical kernel is for plug flow only")
    distance_m = positive(distance_m, "distance", True)
    edges = np.asarray(edges, dtype=float)
    if distance_m == 0:
        return np.diff(edges)/channel.depth_m
    wave = (np.arange(modes)+0.5)*np.pi/channel.depth_m
    decay = np.exp(-channel.diffusivity*wave*wave*distance_m/channel.speed_m_s)
    integral = np.cos(np.outer(wave, edges[:-1]))-np.cos(np.outer(wave, edges[1:]))
    return 2/(channel.depth_m**2)*np.sum((decay/(wave*wave))[:, None]*integral, axis=0)


def modal_witness(channel, station_x_m, start_m, stop_m, cap_K=0.2):
    """Smooth positive continuum counterexample, independent of FV/LP bounds.

    A finite sine combination annihilates all plug-flow bulk observations.
    A proven sine-ratio bound scales it to keep both inlet profiles in [0,cap].
    """
    if channel.velocity != "plug":
        raise ValueError("Continuum witness uses plug flow")
    stations = np.asarray(station_x_m, dtype=float)
    if len(stations) != len(set(stations)) or np.any(stations < 0):
        raise ValueError("Distinct nonnegative station positions required")
    if not 0 <= start_m < stop_m or cap_K <= 0:
        raise ValueError("Invalid target or cap")
    n = len(stations)+2
    k = (np.arange(n)+0.5)*np.pi/channel.depth_m
    lam = channel.diffusivity*k*k/channel.speed_m_s
    beta = 1/(channel.depth_m*k)
    e = np.exp(-np.outer(stations, lam))*beta
    objective = beta*(np.exp(-lam*start_m)-np.exp(-lam*stop_m))
    _, singular, vt = np.linalg.svd(e, full_matrices=True)
    rank = int(np.count_nonzero(singular > singular[0]*1e-13)) if len(singular) else 0
    if rank != len(stations):
        raise ArithmeticError("Observation matrix is numerically unresolved")
    null = vt[rank:].T
    direction = null@(null.T@objective)
    norm = np.linalg.norm(direction)
    if norm < 1e-13*np.linalg.norm(objective):
        raise ArithmeticError("Witness is zero or numerically unresolved")
    direction /= norm
    # |sin((2n+1)t)| <= (2n+1) sin(t), 0<=t<=pi/2.
    amplitude = cap_K/(2*np.sum(np.abs(direction)*(2*np.arange(n)+1)))
    base = np.zeros(n)
    base[0] = cap_K/2
    minus, plus = base-amplitude*direction, base+amplitude*direction
    return dict(mode_indices=list(range(n)), eigenvalues_per_m=lam.tolist(),
                minus_coefficients_K=minus.tolist(), plus_coefficients_K=plus.tolist(),
                calibration_difference_K=(e@(plus-minus)).tolist(),
                regional_power_difference_W=float(channel.capacity_W_K*objective@(plus-minus)),
                amplitude_K=float(amplitude), observation_singular_values=singular.tolist(),
                cap_K=cap_K, provenance="Analytical continuum sine modes; no FV or LP")
