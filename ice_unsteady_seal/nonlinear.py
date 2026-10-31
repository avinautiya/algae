"""Bounded nonlinear finite-interval extension, away from dry/critical flow.

This implements the declared fixed-frame boundary problem. It does not select
between outer shock branches or resolve hydraulic waves. A BVP failure stops
the forecast rather than silently substituting a different physical model.
"""

import numpy as np
from scipy.integrate import solve_bvp
from scipy.interpolate import CubicSpline

from .model import hydraulic_coefficients, total_water_head


class NonlinearSealModel:
    def __init__(self, base, lake_capacity=2., hydraulic_tolerance=1e-8,subcritical_margin=.02):
        self.base=base; self.lake_capacity=float(lake_capacity)
        self.hydraulic_tolerance=float(hydraulic_tolerance)
        self.subcritical_margin=float(subcritical_margin)
        if self.lake_capacity<=0 or self.hydraulic_tolerance<=0:
            raise ValueError("Positive lake capacity and BVP tolerance required")
        if not (0<self.subcritical_margin<1): raise ValueError("Invalid criticality margin")

    def hydraulics(self,head_perturbation,bed_perturbation):
        base=self.base; c=base.config; x=base.x
        bed=np.asarray(bed_perturbation,dtype=float)
        if bed.shape!=x.shape or abs(bed[0])>1e-13:
            raise ValueError("Bed must match the grid and its upstream boundary")
        perturb=CubicSpline(x,bed,bc_type="natural")
        target=base.upstream_head+float(head_perturbation)
        resistance_right=c.discharge**(c.alpha-1)*base.velocity[-1]**(3-c.alpha)
        def rhs(locations,log_velocity,log_discharge):
            q=float(np.exp(log_discharge[0])); v=np.exp(log_velocity[0])
            original=base.evaluate(locations)[0]
            slope=base.upstream_slope-original**3/c.relative_ice_speed+perturb(locations,1)
            d,_,r,_,_=hydraulic_coefficients(v,q,c.alpha,c.beta,c.froude_coefficient)
            return np.asarray([(slope+r)/(d*v)])
        def boundary(left,right,log_discharge):
            q=float(np.exp(log_discharge[0])); vl=float(np.exp(left[0])); vr=float(np.exp(right[0]))
            inlet=base.bed[0]+total_water_head(vl,q,c.beta,c.froude_coefficient)-target
            normal=q**(c.alpha-1)*vr**(3-c.alpha)-resistance_right
            return np.asarray([inlet,normal])
        solution=solve_bvp(rhs,boundary,x,np.log(base.velocity)[None,:],
                           p=np.asarray([np.log(c.discharge)]),
                           tol=self.hydraulic_tolerance,max_nodes=max(4096,16*len(x)))
        if not solution.success:
            raise RuntimeError("Nonlinear hydraulic BVP failed: "+solution.message+
                               "; nodes="+str(len(solution.x))+"; max residual="+
                               str(float(np.max(solution.rms_residuals))))
        q=float(np.exp(solution.p[0])); v=np.exp(solution.sol(x)[0])
        hydraulic_coefficients(v,q,c.alpha,c.beta,c.froude_coefficient)
        fr2=c.froude_coefficient**2*v**2/(c.beta*(q/v)**c.beta)
        if np.max(fr2)>=1-self.subcritical_margin:
            raise ValueError("Nonlinear hydraulic state violates the criticality margin")
        return {"discharge":q,"velocity":v,"boundary_residual":float(np.max(abs(boundary(
            solution.sol(x[0]),solution.sol(x[-1]),solution.p)))),
            "maximum_bvp_residual":float(np.max(solution.rms_residuals)),
            "solver_nodes":len(solution.x),"maximum_local_froude_squared":float(np.max(fr2))}

    def rhs(self,state,inflow_perturbation):
        c=self.base.config
        bed=np.concatenate(([0.],np.asarray(state)[1:]))
        hydro=self.hydraulics(float(state[0]),bed)
        gradient=np.zeros_like(bed)
        gradient[1:]=np.diff(bed)/np.diff(self.base.x)
        rate=-c.relative_ice_speed*gradient-(hydro["velocity"]**3-self.base.velocity**3)
        j=hydro["discharge"]-c.discharge
        return np.concatenate(([(inflow_perturbation-j)/self.lake_capacity],rate[1:])),hydro

    def forecast(self,segments,initial=None,maximum_step=.02,relative_limit=.1):
        """RK4 on the declared spatial model; no automatic branch switching.

        Integrated inflow/outflow use the same RK stages as surrogate
        total-head storage; finite-cutoff stage corrections are separate.
        A CFL cap limits advection steps; it is not a proof of nonlinear
        stability. Refinement is required for a scientific prediction.
        """
        c=self.base.config; n=len(self.base.x)
        y=np.zeros(n) if initial is None else np.asarray(initial,dtype=float).copy()
        if y.shape!=(n,) or maximum_step<=0: raise ValueError("Invalid state or step")
        cfl=.35*np.min(np.diff(self.base.x))/c.relative_ice_speed
        initial_head=float(y[0]); t=0.; accumulated_in=0.; accumulated_out=0.
        samples=[]
        initial_stage=None
        def record(hydro):
            nonlocal initial_stage
            bed=np.concatenate(([0.],y[1:]))
            j=hydro["discharge"]-c.discharge
            perturbation=max(abs(j)/c.discharge,
                             float(np.max(abs(hydro["velocity"]-self.base.velocity)/self.base.velocity)),
                             float(np.max(abs(bed)/(c.discharge/self.base.velocity)**c.beta)))
            kinetic=.5*c.froude_coefficient**2*(hydro["velocity"][0]**2-self.base.velocity[0]**2)
            stage=float(y[0]-kinetic)
            if initial_stage is None: initial_stage=stage
            samples.append({"time":t,"lake_head_perturbation":float(y[0]),
                            "discharge":hydro["discharge"],"discharge_perturbation":j,
                            "total_head_ledger_residual":self.lake_capacity*(y[0]-initial_head)-accumulated_in+accumulated_out,
                            "upstream_surface_stage_perturbation":stage,
                            "upstream_kinetic_head_perturbation":float(kinetic),
                            "surface_stage_storage_residual":self.lake_capacity*(stage-initial_stage)-accumulated_in+accumulated_out,
                            "maximum_bvp_residual":hydro["maximum_bvp_residual"],
                            "maximum_relative_perturbation":perturbation,
                            "within_declared_small_perturbation_screen":perturbation<=relative_limit})
        _,first=self.rhs(y,0.); record(first)
        for duration,inflow in segments:
            duration,inflow=float(duration),float(inflow)
            if duration<=0 or c.discharge+inflow<0: raise ValueError("Invalid water input")
            count=max(1,int(np.ceil(duration/min(maximum_step,cfl)))); dt=duration/count
            for _ in range(count):
                k1,h1=self.rhs(y,inflow)
                k2,h2=self.rhs(y+.5*dt*k1,inflow)
                k3,h3=self.rhs(y+.5*dt*k2,inflow)
                k4,h4=self.rhs(y+dt*k3,inflow)
                y+=dt*(k1+2*k2+2*k3+k4)/6
                accumulated_in+=dt*inflow
                accumulated_out+=dt*((h1["discharge"]-c.discharge)+2*(h2["discharge"]-c.discharge)
                                     +2*(h3["discharge"]-c.discharge)+(h4["discharge"]-c.discharge))/6
                t+=dt
                _,hydro=self.rhs(y,inflow); record(hydro)
        return {"model":"finite_interval_nonlinear_unsteady_seal","units":"dimensionless",
                "scope":"Fixed reference; subcritical hydraulic BVP; no outer branch selection",
                "base":self.base.to_dict(),"lake_capacity":self.lake_capacity,
                "storage_scope":"Surrogate reservoir total head; finite-cutoff pond matching unverified",
                "maximum_total_head_ledger_residual":max(abs(s["total_head_ledger_residual"]) for s in samples),
                "samples":samples,"final_bed_perturbation":np.concatenate(([0.],y[1:])).tolist()}
