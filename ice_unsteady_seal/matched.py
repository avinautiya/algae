"""Canonical outer-gauge linear closure with bed translation removed.

Discretizes MATCHED_LINEAR_CLOSURE.md. Only alpha=beta=.5, q0=Vplus=1,
a=.5 is implemented. Finite cutoffs approximate the proved real-line law.
No nonlinear branch selection or environmental validation is asserted.
"""

import numpy as np
from scipy.linalg import expm

from .model import HydraulicMap


class MatchedSealModel:
    def __init__(self,base,lake_capacity=2.0):
        c=base.config
        if not np.allclose([c.alpha,c.beta,c.discharge,c.downstream_velocity,c.relative_ice_speed],
                           [.5,.5,1.,1.,.5],rtol=0,atol=1e-14):
            raise ValueError("Matched closure is implemented only for the canonical seal")
        if lake_capacity<=0: raise ValueError("Positive lake capacity required")
        self.base=base; self.lake_capacity=float(lake_capacity)
        h=HydraulicMap(base.x,base.velocity,1.,.5,.5,c.froude_coefficient)
        self.hydraulics=h
        n=len(base.x); dx=np.zeros((n,n))
        for i in range(1,n):
            width=base.x[i]-base.x[i-1]; dx[i,i]=1/width; dx[i,i-1]=-1/width
        m=3*base.velocity**2
        l0=-.5*dx-m[:,None]*h.velocity_shape
        gs=-m*(h.velocity_flux-base.velocity/5)
        # Both shape tails are zero in the real-line model; crop endpoints
        # are fixed to zero here. This needs domain/refinement convergence.
        self.transport_shape=l0[1:-1,1:-1]
        self.flux_to_shape=gs[1:-1]
        self.shape_output=h.bed_output[1:-1]
        self.direct=h.direct
        self.output=np.concatenate(([self.direct],self.shape_output))
        self.mu=1/self.lake_capacity-.3
        size=n-1  # relative head, N-2 interior shape values
        self.operator=np.zeros((size,size))
        self.operator[0]=-self.mu*self.output
        self.operator[1:,0]=self.flux_to_shape*self.direct
        self.operator[1:,1:]=self.transport_shape+self.flux_to_shape[:,None]*self.shape_output
        self.input=np.zeros(size); self.input[0]=1/self.lake_capacity
        w=base.velocity**3
        self.matching_diagnostics={
            "hydraulic_translation_identity_defect":float(h.bed_output@w+self.direct/2),
            "velocity_translation_identity_maximum_defect":float(np.max(abs(
                h.velocity_shape@w+.5*np.gradient(base.velocity,base.x,edge_order=2)))),
            "largest_positive_shape_forcing_coefficient":float(max(0.,np.max(self.flux_to_shape))),
            "scope":"Finite quadrature/crop defects; continuum identities are exact in the proved canonical limit"}

    def initial_state(self,relative_head=0.,shape=None):
        s=np.zeros(len(self.base.x)) if shape is None else np.asarray(shape,dtype=float)
        if s.shape!=self.base.x.shape or abs(s[0])+abs(s[-1])>1e-13:
            raise ValueError("Shape must match the grid and vanish at both crop endpoints")
        return np.concatenate(([relative_head],s[1:-1]))

    def renewal_kernel(self,times):
        """C exp(L0 t) g_s, the derived discharge-feedback kernel."""
        return np.asarray([float(self.shape_output@expm(self.transport_shape*float(t))@self.flux_to_shape)
                           for t in times])

    def forecast(self,segments,initial=None,sample_step=.1,relative_limit=.1):
        n=self.operator.shape[0]
        state=np.zeros(n) if initial is None else np.asarray(initial,dtype=float)
        if state.shape!=(n,) or sample_step<=0: raise ValueError("Invalid state/time step")
        # Auxiliary states: r=downstream bed offset, integral inflow,
        # integral outflow, constant 1. The outer shift/height is r/2.
        y=np.concatenate((state,[0.,0.,0.,1.]))
        initial_total_head=float(state[0]); t=0.; raw=[(t,y.copy())]; cache={}
        for duration,inflow in segments:
            duration,inflow=float(duration),float(inflow)
            if duration<=0 or 1+inflow<0: raise ValueError("Invalid forcing")
            steps=max(1,int(np.ceil(duration/sample_step))); dt=duration/steps
            key=(dt,inflow)
            if key not in cache:
                a=np.zeros((n+4,n+4)); a[:n,:n]=self.operator; a[:n,-1]=self.input*inflow
                a[n,:n]=-.6*self.output; a[n+1,-1]=inflow; a[n+2,:n]=self.output
                cache[key]=expm(a*dt)
            for _ in range(steps): y=cache[key]@y; t+=dt; raw.append((t,y.copy()))
        rows=[]; c=self.base.config; crest=int(np.argmin(abs(self.base.x)))
        for time,z in raw:
            s=np.concatenate(([0.],z[1:n],[0.])); offset=float(z[n]); j=float(self.output@z[:n])
            total_head=float(z[0]+offset/2)
            bed=s+offset*self.base.velocity**3
            replay,_,_=self.hydraulics.solve(total_head,bed)
            # Evaluate the matched state with its own discharge. Cropped
            # head matching can produce a distinct replay discharge.
            v=self.hydraulics.velocity_shape@bed+self.hydraulics.velocity_flux*j
            positive=bool(1+j>0 and np.all(self.base.velocity+v>0))
            fr2=float(np.max(c.froude_coefficient**2*(self.base.velocity+v)**2/
                              (.5*((1+j)/(self.base.velocity+v))**.5))) if positive else None
            amplitude=max(abs(j),float(np.max(abs(v)/self.base.velocity)),
                          float(np.max(abs(bed)/(1/self.base.velocity)**.5)))
            ledger=self.lake_capacity*(total_head-initial_total_head)-z[n+1]+z[n+2]
            rows.append({"time":time,"discharge":1+j,"discharge_perturbation":j,
                         "outer_seal_relative_head":float(z[0]),"reservoir_head_perturbation":total_head,
                         "outer_seal_height_shift":offset/2,"outer_seal_location_shift":offset/2,
                         "inner_crest_height_perturbation":float(s[crest]+offset/2),
                         "inner_crest_relative_head":float(z[0]-s[crest]),
                         "matched_limit_storage_ledger_residual":float(ledger),
                         "finite_crop_hydraulic_replay_difference":float(replay-j),
                         "maximum_relative_perturbation":amplitude,
                         "maximum_perturbed_local_froude_squared":fr2,
                         "inside_declared_screen":bool(positive and fr2<.98 and amplitude<=relative_limit)})
        return {"model":"canonical_outer_matched_linear_seal","units":"dimensionless",
                "base":self.base.to_dict(),"lake_capacity":self.lake_capacity,
                "matching_diagnostics":self.matching_diagnostics,
                "matching_scope":"Discretized proved canonical linear real-line model; nonlinear selection remains open",
                "storage_scope":"Head is pond stage in the ideal matched limit; finite crop/kinetic corrections remain",
                "maximum_matched_limit_storage_ledger_residual":max(abs(r['matched_limit_storage_ledger_residual']) for r in rows),
                "samples":rows,"final_shape":np.concatenate(([0.],y[1:n],[0.])).tolist()}
