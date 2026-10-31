"""Derived hydraulic Schur complement and evolving inner bed-state operator.

The prototype retains the inner bed time derivative. Fluid hydraulics remain
quasisteady and uniformly subcritical. Upstream pond matching is truncated;
nonlinear branch selection is outside this implementation. See MATHEMATICS.md.
"""

from dataclasses import asdict, dataclass
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm


@dataclass(frozen=True)
class BaseConfig:
    alpha: float = 0.5
    beta: float = 0.5
    froude_coefficient: float = 0.3
    discharge: float = 1.0
    downstream_velocity: float = 1.0
    relative_ice_speed: float = 0.5
    left: float = -8.0
    right: float = 8.0
    points: int = 129


def hydraulic_coefficients(velocity, q, alpha, beta, froude):
    v = np.asarray(velocity, dtype=float)
    if np.any(v <= 0) or q <= 0 or not (0 <= alpha < 1) or beta <= 0:
        raise ValueError("Positive wet states and 0 <= alpha < 1, beta > 0 required")
    depth = (q / v) ** beta
    d = beta * depth / v - froude**2 * v  # minus head derivative H_V
    if np.any(d <= 0):
        raise ValueError("Strict subcriticality D=-H_V>0 is required")
    head_q = beta * depth / q
    resistance = q ** (alpha - 1) * v ** (3 - alpha)
    resistance_v = (3 - alpha) * resistance / v
    resistance_q = (alpha - 1) * resistance / q
    return d, head_q, resistance, resistance_v, resistance_q


def total_water_head(v, q, beta, froude):
    return (q / np.asarray(v)) ** beta + 0.5 * froude**2 * np.asarray(v)**2


@dataclass
class BaseState:
    config: BaseConfig
    x: np.ndarray
    velocity: np.ndarray
    bed: np.ndarray
    upstream_slope: float
    crest_velocity: float
    crest_curvature: float
    upstream_head: float
    evaluate: Callable

    def to_dict(self):
        c = self.config
        local_fr2 = c.froude_coefficient**2 * self.velocity**2 / (
            c.beta * (c.discharge / self.velocity)**c.beta
        )
        right_resistance=c.discharge**(c.alpha-1)*self.velocity[-1]**(3-c.alpha)
        normal_residual=self.upstream_slope-self.velocity[-1]**3/c.relative_ice_speed+right_resistance
        return {
            "configuration": asdict(c), "x": self.x.tolist(),
            "velocity": self.velocity.tolist(), "bed": self.bed.tolist(),
            "upstream_slope": self.upstream_slope,
            "crest_velocity": self.crest_velocity,
            "crest_curvature": self.crest_curvature,
            "upstream_total_head": self.upstream_head,
            "maximum_local_froude_squared": float(np.max(local_fr2)),
            "downstream_normal_flow_residual":float(normal_residual),
            "downstream_velocity_distance_to_asymptote":float(c.downstream_velocity-self.velocity[-1]),
            "upstream_base_kinetic_head":float(.5*c.froude_coefficient**2*self.velocity[0]**2),
            "scope": "Dimensionless finite positive-velocity subcritical interval",
        }


def build_base(config=BaseConfig()):
    """Integrate a nonuniform travelling seal with its smooth crest at X=0.

    The normal-flow downstream velocity and relative ice speed set p_minus.
    This is a truncated base, not proof of an infinite-domain matching limit.
    """
    c = config
    if c.points < 9 or not (c.left < 0 < c.right) or c.relative_ice_speed <= 0:
        raise ValueError("Need >=9 points, left<0<right and positive relative speed")
    vp, q, speed = c.downstream_velocity, c.discharge, c.relative_ice_speed
    hydraulic_coefficients(vp, q, c.alpha, c.beta, c.froude_coefficient)
    rp = q**(c.alpha - 1) * vp**(3 - c.alpha)
    pminus = vp**3 / speed - rp
    if pminus <= 0:
        raise ValueError("Parameters do not form a pond-to-flow seal: p_minus<=0")
    vc = (speed * pminus)**(1.0 / 3.0)
    if vc >= vp:
        raise ValueError("The crest velocity must lie below the downstream velocity")

    def rhs(x, state):
        v = float(state[0])
        d, _, r, _, _ = hydraulic_coefficients(
            v, q, c.alpha, c.beta, c.froude_coefficient
        )
        slope = pminus - v**3 / speed
        return [(slope + float(r)) / float(d), slope]

    options = dict(rtol=2e-11, atol=2e-13, dense_output=True, max_step=0.15)
    backward = solve_ivp(rhs, (0, c.left), (vc, 0), **options)
    forward = solve_ivp(rhs, (0, c.right), (vc, 0), **options)
    if not backward.success or not forward.success:
        raise RuntimeError("Nonuniform base integration failed")
    nl = max(2, min(c.points-3, round((c.points-1)*(-c.left)/(c.right-c.left))))
    nr = c.points-1-nl
    x = np.concatenate((np.linspace(c.left, 0, nl+1), np.linspace(0, c.right, nr+1)[1:]))

    def evaluate(locations):
        z = np.atleast_1d(np.asarray(locations, dtype=float))
        if np.any(z < c.left) or np.any(z > c.right):
            raise ValueError("Evaluation outside declared matching interval")
        result = np.empty((2, len(z)))
        negative = z < 0
        if np.any(negative): result[:, negative] = backward.sol(z[negative])
        if np.any(~negative): result[:, ~negative] = forward.sol(z[~negative])
        return result

    values = evaluate(x)
    dcrest = float(hydraulic_coefficients(vc, q, c.alpha, c.beta, c.froude_coefficient)[0])
    rcrest = q**(c.alpha-1) * vc**(3-c.alpha)
    curvature = -3*vc**2 * (rcrest/dcrest) / speed
    upstream_head = float(values[1,0] + total_water_head(
        values[0,0], q, c.beta, c.froude_coefficient
    ))
    return BaseState(c, x, values[0], values[1], pminus, vc, curvature,
                     upstream_head, evaluate)


def _moments(a, width):
    """Stable integrals of exp(-a s) and s exp(-a s) on a cell."""
    t = a * width
    if abs(t) < 1e-4:
        j0 = width * (1-t/2+t*t/6-t**3/24+t**4/120)
        j1 = width**2 * (0.5-t/3+t*t/8-t**3/30+t**4/144)
    else:
        j0 = -np.expm1(-t) / a
        j1 = (1-(1+t)*np.exp(-t)) / a**2
    return float(np.exp(-t)), float(j0), float(j1)


class HydraulicMap:
    """Linearized nonuniform hydraulic BVP with normal downstream matching.

    e=b-D*v+H_q*j; e_X=-R_V*v-R_q*j; v(L)=ell*j.
    The upstream total-head perturbation is e(left)=h.
    Coefficients are midpoint-constant per cell, bed perturbations linear.
    """
    def __init__(self, x, velocity, q, alpha, beta, froude):
        self.x = np.asarray(x, dtype=float)
        self.velocity = np.asarray(velocity, dtype=float)
        if self.x.ndim != 1 or self.velocity.shape != self.x.shape or len(self.x)<3:
            raise ValueError("Matching one-dimensional state arrays required")
        if np.any(np.diff(self.x)<=0): raise ValueError("Grid must increase")
        d, hq, r, rv, rq = hydraulic_coefficients(velocity,q,alpha,beta,froude)
        self.d, self.head_q = d, hq
        self.normal_velocity_derivative = (1-alpha)*self.velocity[-1]/((3-alpha)*q)
        a = rv/d
        f = rv*hq/d + rq
        z = hq[-1] - d[-1]*self.normal_velocity_derivative
        if z <= 0: raise ValueError("Nonpositive downstream head derivative")
        n = len(self.x)
        k = np.zeros((n,n)); p = np.zeros(n)
        k[-1,-1] = 1; p[-1] = z
        for i in range(n-2,-1,-1):
            width = self.x[i+1]-self.x[i]
            amid = 0.5*(a[i]+a[i+1]); fmid=0.5*(f[i]+f[i+1])
            weight,j0,j1 = _moments(amid,width)
            k[i] = weight*k[i+1]
            k[i,i] += amid*(j0-j1/width)
            k[i,i+1] += amid*j1/width
            p[i] = weight*p[i+1] + fmid*j0
        if not np.isfinite(p[0]) or p[0]<=0:
            raise ValueError("Hydraulic Schur complement is not positive")
        self.k, self.p = k,p
        self.direct = 1/p[0]
        self.bed_output = -k[0]/p[0]
        self.velocity_shape = (np.eye(n)-k)/d[:,None]
        self.velocity_flux = (hq-p)/d
        self.velocity_bed = self.velocity_shape + self.velocity_flux[:,None]*self.bed_output
        self.velocity_head = self.velocity_flux*self.direct

    def solve(self, head, bed):
        b = np.asarray(bed,dtype=float)
        if b.shape != self.x.shape: raise ValueError("Wrong bed-state length")
        j = float(self.direct*head+self.bed_output@b)
        e = self.k@b+self.p*j
        v = self.velocity_shape@b+self.velocity_flux*j
        return j,e,v


class SealModel:
    """Fixed-reference-frame linear bed/lake forecast around a nonuniform seal.

    b_T+a*b_X=-3 V_*^2*v; gamma*h_T=i-j. The bed reference is the
    prescribed travelling base, not the evolving crest. Upstream b(left)=0
    truncates a pond with zero leading-order incision. This matching limit
    remains to be established; this is NOT the full nonlinear seal closure.
    """
    def __init__(self, base, lake_capacity=2.0, subcritical_margin=0.02):
        if lake_capacity<=0: raise ValueError("Positive lake capacity required")
        if not (0<subcritical_margin<1): raise ValueError("Criticality margin must lie in (0,1)")
        self.base, self.lake_capacity = base,float(lake_capacity)
        self.subcritical_margin=float(subcritical_margin)
        c = base.config; n=len(base.x)
        if base.to_dict()["maximum_local_froude_squared"]>=1-self.subcritical_margin:
            raise ValueError("Base is outside the declared subcritical margin")
        hmap = HydraulicMap(base.x,base.velocity,c.discharge,c.alpha,c.beta,c.froude_coefficient)
        self.hydraulics=hmap
        dx=np.zeros((n,n))
        for i in range(1,n):
            gap=base.x[i]-base.x[i-1]
            dx[i,i]=1/gap; dx[i,i-1]=-1/gap
        melt_derivative=3*base.velocity**2
        full_l=-c.relative_ice_speed*dx-melt_derivative[:,None]*hmap.velocity_bed
        full_g=-melt_derivative*hmap.velocity_head
        self.bed_operator=full_l[1:,1:]
        self.head_to_bed=full_g[1:]
        self.output=np.concatenate(([hmap.direct],hmap.bed_output[1:]))
        self.operator=np.zeros((n,n))
        self.operator[0]=-self.output/self.lake_capacity
        self.operator[1:,0]=self.head_to_bed
        self.operator[1:,1:]=self.bed_operator
        self.input=np.zeros(n); self.input[0]=1/self.lake_capacity

    def initial_state(self, head=0.0, bed=None):
        b=np.zeros(len(self.base.x)) if bed is None else np.asarray(bed,dtype=float)
        if b.shape != self.base.x.shape or abs(b[0])>1e-14:
            raise ValueError("Bed state must match grid and have b(left)=0")
        return np.concatenate(([head],b[1:]))

    def bed_from_state(self,state):
        return np.concatenate(([0.0],np.asarray(state)[1:]))

    def head_input_kernel(self,times):
        """Prescribed-head kernel C exp(L t) G; excludes lake feedback."""
        c=self.hydraulics.bed_output[1:]
        return np.asarray([float(c@expm(self.bed_operator*float(t))@self.head_to_bed)
                           for t in times])

    def forecast(self,segments,initial=None,sample_step=0.1,relative_limit=0.1):
        """Piecewise-constant inflow forecast, using exact matrix exponentials.

        No clipping, hidden fitting, or automatic claim of physical validity.
        Integrals provide a surrogate reservoir-total-head ledger. The finite
        cutoff surface-stage correction is reported separately; roundoff in
        this surrogate is not proof of physical lake mass closure.
        """
        if sample_step<=0 or relative_limit<=0: raise ValueError("Positive steps/limits required")
        n=self.operator.shape[0]
        state=np.zeros(n) if initial is None else np.asarray(initial,dtype=float)
        if state.shape!=(n,): raise ValueError("Wrong initial-state shape")
        augmented=np.concatenate((state,[0.,0.,1.]))
        initial_head=float(state[0]); t=0.; snapshots=[(t,augmented.copy())]
        cache={}
        for duration,inflow in segments:
            duration,inflow=float(duration),float(inflow)
            if duration<=0 or self.base.config.discharge+inflow<0:
                raise ValueError("Positive durations and nonnegative total water input required")
            count=max(1,int(np.ceil(duration/sample_step))); dt=duration/count
            key=(dt,inflow)
            if key not in cache:
                a=np.zeros((n+3,n+3)); a[:n,:n]=self.operator
                a[:n,-1]=self.input*inflow
                a[n,-1]=inflow; a[n+1,:n]=self.output
                cache[key]=expm(a*dt)
            for _ in range(count):
                augmented=cache[key]@augmented
                t+=dt; snapshots.append((t,augmented.copy()))
        output=[]; c=self.base.config
        crest=int(np.argmin(abs(self.base.x)))
        initial_bed=self.bed_from_state(state)
        initial_velocity=self.hydraulics.solve(initial_head,initial_bed)[2]
        initial_stage=initial_head-c.froude_coefficient**2*self.base.velocity[0]*initial_velocity[0]
        for time,vector in snapshots:
            y=vector[:n]; b=self.bed_from_state(y)
            j,_,v=self.hydraulics.solve(y[0],b)
            derivative=np.gradient(b,self.base.x,edge_order=2)
            shift=-derivative[crest]/self.base.crest_curvature
            ledger=self.lake_capacity*(y[0]-initial_head)-(vector[n]-vector[n+1])
            kinetic=c.froude_coefficient**2*self.base.velocity[0]*v[0]
            stage=float(y[0]-kinetic)
            stage_ledger=self.lake_capacity*(stage-initial_stage)-(vector[n]-vector[n+1])
            ratios={"discharge":abs(j)/c.discharge,
                    "velocity":float(np.max(abs(v)/self.base.velocity)),
                    "bed_over_base_depth":float(np.max(abs(b)/(c.discharge/self.base.velocity)**c.beta))}
            positive=bool(c.discharge+j>0 and np.all(self.base.velocity+v>0))
            if positive:
                perturbed_fr2=c.froude_coefficient**2*(self.base.velocity+v)**2/(
                    c.beta*((c.discharge+j)/(self.base.velocity+v))**c.beta)
                maximum_fr2=float(np.max(perturbed_fr2))
            else: maximum_fr2=None
            critical_ok=positive and maximum_fr2<1-self.subcritical_margin
            valid=all(value<=relative_limit for value in ratios.values()) and critical_ok
            output.append({"time":time,"lake_head_perturbation":float(y[0]),
                           "discharge_perturbation":j,"discharge":c.discharge+j,
                           "crest_height_perturbation":float(b[crest]),
                           "linear_crest_location_shift":float(shift),
                           "total_head_ledger_residual":float(ledger),
                           "upstream_kinetic_head_perturbation":float(kinetic),
                           "upstream_surface_stage_perturbation":stage,
                           "surface_stage_storage_residual":float(stage_ledger),
                           "maximum_perturbed_local_froude_squared":maximum_fr2,
                           "inside_declared_subcritical_margin":bool(critical_ok),
                           "positive_discharge_and_velocity":positive,
                           "inflow_perturbation_integral":float(vector[n]),
                           "outflow_perturbation_integral":float(vector[n+1]),
                           "relative_perturbations":ratios,"within_declared_linear_screen":bool(valid)})
        return {"model":"finite_interval_linear_unsteady_seal",
                "units":"dimensionless","base":self.base.to_dict(),
                "lake_capacity":self.lake_capacity,"relative_screen":relative_limit,
                "subcritical_margin":self.subcritical_margin,
                "storage_scope":"Surrogate reservoir total-head state; finite-cutoff stage/pond matching not certified",
                "matching_status":"Canonical hydraulic normalization limit proved; full nonlinear dynamic matching and discrete cutoff convergence not established",
                "maximum_total_head_ledger_residual":max(abs(s["total_head_ledger_residual"]) for s in output),
                "maximum_surface_stage_storage_residual":max(abs(s["surface_stage_storage_residual"]) for s in output),
                "all_samples_within_declared_linear_screen":all(s["within_declared_linear_screen"] for s in output),
                "samples":output,"final_bed_perturbation":self.bed_from_state(augmented[:n]).tolist()}
