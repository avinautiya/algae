"""Classical analytical benchmarks implemented independently in standard Python.

These are known solutions, not new physics. A VOF calculation must establish
that its geometry/regime and observation operator match each benchmark.
"""
from dataclasses import dataclass
from math import erf, exp, expm1, isfinite, pi, sqrt
from functools import cached_property


def valid(x, name, minimum=None, strict=False):
    if not isfinite(x) or (minimum is not None and (x < minimum or (strict and x == minimum))):
        raise ValueError(f"Invalid {name}: {x}")
    return x


def rectangular_jump(depth_up_m, discharge_m3_s, width_m, gravity=9.81, rho=1000.0):
    """Stationary, horizontal, hydrostatic rectangular jump; unit profile factors.

Negligible integrated bed/side friction between measurement sections. The
relation does not predict jump position, length, turbulence or wall heat flux.
"""
    for key, value in (("depth", depth_up_m), ("discharge", discharge_m3_s),
                       ("width", width_m), ("gravity", gravity), ("density", rho)):
        valid(value, key, 0.0, True)
    h, q = depth_up_m, discharge_m3_s/width_m
    u = q/h
    fr = u/sqrt(gravity*h)
    if fr <= 1.0:
        raise ValueError("Jump reference requires a supercritical upstream state")
    h2 = 0.5*h*(sqrt(1.0+8.0*fr*fr)-1.0)
    u2 = q/h2
    e1, e2 = h+u*u/(2*gravity), h2+u2*u2/(2*gravity)
    loss = (h2-h)**3/(4*h*h2)
    return dict(depth_up_m=h, depth_down_m=h2, discharge_m3_s=discharge_m3_s,
                width_m=width_m, velocity_up_m_s=u, velocity_down_m_s=u2,
                Fr_up=fr, Fr_down=u2/sqrt(gravity*h2), specific_energy_up_m=e1,
                specific_energy_down_m=e2, head_loss_m=loss,
                momentum_up_N=rho*width_m*(q*q/h+0.5*gravity*h*h),
                momentum_down_N=rho*width_m*(q*q/h2+0.5*gravity*h2*h2),
                mechanical_power_loss_W=rho*gravity*discharge_m3_s*loss)


@dataclass(frozen=True)
class DamBreak:
    """Initially stationary, flat-bed shallow-water dam break, left depth > right.

Supports a wet or initially dry downstream state. The hydrostatic, inviscid
depth-averaged reference is not an exact solution of short-wave VOF dynamics.
Use remote-from-gate regions and suitable shallow-water scaling for comparison.
"""
    left_depth_m: float = 1.0
    right_depth_m: float = 0.1
    gravity: float = 9.81

    def __post_init__(self):
        valid(self.left_depth_m, "left depth", 0.0, True)
        valid(self.right_depth_m, "right depth", 0.0)
        valid(self.gravity, "gravity", 0.0, True)
        if self.right_depth_m >= self.left_depth_m:
            raise ValueError("Require left depth > right depth")

    @cached_property
    def middle(self):
        if self.right_depth_m == 0.0:
            return 0.0, 2*sqrt(self.gravity*self.left_depth_m)
        lo, hi = self.right_depth_m, self.left_depth_m
        c0 = sqrt(self.gravity*self.left_depth_m)
        for _ in range(100):
            h = (lo+hi)/2
            ul = 2*(c0-sqrt(self.gravity*h))
            ur = (h-self.right_depth_m)*sqrt(0.5*self.gravity*(1/h+1/self.right_depth_m))
            if ul > ur:
                lo = h
            else:
                hi = h
        h = (lo+hi)/2
        return h, 2*(c0-sqrt(self.gravity*h))

    @cached_property
    def speeds(self):
        c0 = sqrt(self.gravity*self.left_depth_m)
        hs, us = self.middle
        if self.right_depth_m == 0.0:
            return (-c0, 2*c0)
        return (-c0, us-sqrt(self.gravity*hs), hs*us/(hs-self.right_depth_m))

    def point(self, x_m, time_s):
        valid(x_m, "x")
        valid(time_s, "time", 0.0)
        if time_s == 0:
            h = self.left_depth_m if x_m < 0 else self.right_depth_m
            return h, 0.0
        xi = x_m/time_s
        speeds = self.speeds
        c0 = sqrt(self.gravity*self.left_depth_m)
        if xi < speeds[0]:
            return self.left_depth_m, 0.0
        if xi < speeds[1]:
            h = (2*c0-xi)**2/(9*self.gravity)
            u = 2*(c0+xi)/3
            return h, h*u
        if self.right_depth_m == 0.0:
            return 0.0, 0.0
        if xi < speeds[2]:
            hs, us = self.middle
            return hs, hs*us
        return self.right_depth_m, 0.0

    def cell_average(self, left_m, right_m, time_s):
        """Average h and hu, splitting every wave; two-point Gauss is exact here.

The fan h is quadratic and hu cubic. This avoids comparing cell averages to
point values across a discontinuity, which would bias refinement measurements.
"""
        valid(left_m, "cell left")
        valid(right_m, "cell right", left_m, True)
        valid(time_s, "time", 0.0)
        waves = [s*time_s for s in self.speeds] if time_s else [0.0]
        cuts = [left_m]+sorted({x for x in waves if left_m < x < right_m})+[right_m]
        hsum = qsum = 0.0
        for a,b in zip(cuts,cuts[1:]):
            mid, half = (a+b)/2, (b-a)/2
            for x in (mid-half/sqrt(3),mid+half/sqrt(3)):
                h,q = self.point(x,time_s)
                hsum += half*h; qsum += half*q
        return hsum/(right_m-left_m), qsum/(right_m-left_m)


@dataclass(frozen=True)
class StefanReference:
    """One-phase, equal-density, conduction-only melting into temperate solid.

Wall x=0 is held at positive temperature excess; liquid occupies 0<x<s(t).
This is an equal-density verification limit, not an unequal-density ice-flow
solution. It verifies no VOF, advective or moving-mesh implementation by itself.
"""
    stefan_number: float
    rho_kg_m3: float = 1000.0
    cp_J_kg_K: float = 4217.0
    conductivity_W_m_K: float = 0.561
    latent_heat_J_kg: float = 334000.0

    def __post_init__(self):
        for name, value in vars(self).items():
            valid(value, name, 0.0, True)
        if not 1e-12 <= self.stefan_number <= 10.0:
            raise ValueError("Supported sensible/latent Stefan number range is [1e-12,10]")

    @property
    def diffusivity(self):
        return self.conductivity_W_m_K/(self.rho_kg_m3*self.cp_J_kg_K)

    @property
    def wall_excess_K(self):
        return self.stefan_number*self.latent_heat_J_kg/self.cp_J_kg_K

    @cached_property
    def lambda_value(self):
        lo, hi = 0.0, 2.0
        for _ in range(100):
            mid = (lo+hi)/2
            if sqrt(pi)*mid*exp(mid*mid)*erf(mid) > self.stefan_number:
                hi = mid
            else:
                lo = mid
        return (lo+hi)/2

    def front(self, time_s):
        valid(time_s, "time", 0.0)
        return 2*self.lambda_value*sqrt(self.diffusivity*time_s)

    def temperature(self, x_m, time_s):
        valid(x_m, "distance", 0.0)
        valid(time_s, "time", 0.0, True)
        if x_m >= self.front(time_s):
            return 0.0
        return self.wall_excess_K*(1-erf(x_m/(2*sqrt(self.diffusivity*time_s)))/erf(self.lambda_value))

    def interface_flux(self, time_s):
        valid(time_s, "time", 0.0, True)
        lam = self.lambda_value
        return self.conductivity_W_m_K*self.wall_excess_K*exp(-lam*lam)/(sqrt(pi*self.diffusivity*time_s)*erf(lam))

    def ledger_per_area(self, time_s):
        """Analytical cumulative wall input, liquid storage and latent energy, J/m²."""
        valid(time_s, "time", 0.0)
        lam, at = self.lambda_value, sqrt(self.diffusivity*time_s)
        ein = 2*self.conductivity_W_m_K*self.wall_excess_K*sqrt(time_s/(pi*self.diffusivity))/erf(lam)
        sensible = (self.rho_kg_m3*self.cp_J_kg_K*self.wall_excess_K*2*at
                    * (-expm1(-lam*lam))/(sqrt(pi)*erf(lam)))
        latent = self.rho_kg_m3*self.latent_heat_J_kg*self.front(time_s)
        return dict(wall_input_J_m2=ein, liquid_storage_J_m2=sensible, latent_J_m2=latent,
                    residual_J_m2=ein-sensible-latent)
