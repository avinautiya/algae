"""Supported bare-ice v5 builder with physical divider width independent of dx.

The exact vendored hydraulic kernel is unchanged. This version replaces the
geometry/initial state before evolution, not an evolved simulation. No claim
of event convergence or environmental validity follows from construction.
"""
import math
import numpy as np
if __package__:
 from .vendor.swe2d_v5 import Sim as PhaseSim
else:
 from vendor.swe2d_v5 import Sim as PhaseSim

def cells(length,dx,name):
 if isinstance(length,bool) or not isinstance(length,(float,int)) or not math.isfinite(length) or length<=0:raise ValueError(name+' must be positive finite')
 n=round(length/dx)
 if n<1 or abs(length/dx-n)>1e-10:raise ValueError(name+' must align exactly with dx')
 return n

class FixedGeometrySim(PhaseSim):
 def __init__(self,*,dx=.0625,divider_width_m=.25,Wl=.25,Lpool=.5,Lr=.5,Ltail=.25,Q=.07,Tw=.5,heff=200.,cfl=.4,manning=.02,z_toe=-.08,z_byp=0.,pool_floor=-1.,tail_drop=1.5,tail_slope=.5,eta0=.15,min_lane_cells=4):
  if isinstance(dx,bool) or not isinstance(dx,(float,int)) or not math.isfinite(dx) or dx<=0:raise ValueError('dx must be positive finite')
  for name,value in [('z_toe',z_toe),('z_byp',z_byp),('pool_floor',pool_floor),('tail_drop',tail_drop),('tail_slope',tail_slope),('eta0',eta0)]:
   if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value):raise ValueError(name+' must be finite')
  if tail_drop<0 or tail_slope<0:raise ValueError('tail drop/slope must be nonnegative')
  self.divider_cells=cells(divider_width_m,dx,'divider width')
  nw=cells(Wl,dx,'lane width');npool=cells(Lpool,dx,'pool length');nr=cells(Lr,dx,'reach length');nt=cells(Ltail,dx,'tail length')
  if type(min_lane_cells) is not int or min_lane_cells<1:raise ValueError('min_lane_cells must be a positive integer')
  if nw<min_lane_cells:raise ValueError('lane is below declared minimum cell count; coarse geometry checks require explicit override')
  self.mesh_contract=dict(lane_cells=nw,declared_minimum_lane_cells=min_lane_cells,convergence_established=False)
  super().__init__(dx=dx,Lpool=Lpool,Lr=Lr,Ltail=Ltail,Wl=Wl,Q=Q,Tw=Tw,heff=heff,cfl=cfl,manning=manning,z_toe=z_toe,z_byp=z_byp,pool_floor=pool_floor,tail_drop=tail_drop,tail_slope=tail_slope,ramp_ratio=0.)
  self.nx=npool+nr+nt;self.nw=nw;self.ny=2*nw+self.divider_cells;self.i0=npool+1;self.i1=npool+nr
  shape=(self.nx+2,self.ny+2)
  self.zb=np.zeros(shape);self.zb[:self.i0]=pool_floor
  self.zd=np.zeros(shape);self.pdep=np.zeros(shape)
  self.wall=np.zeros(shape,dtype=bool);self.ice=np.zeros(shape,dtype=bool);self.toe=np.zeros(shape,dtype=bool);self.byp=np.zeros(shape,dtype=bool);self.solid=np.zeros(shape,dtype=bool)
  self.h=np.zeros(shape);self.hu=np.zeros(shape);self.hv=np.zeros(shape);self.hT=np.zeros(shape)
  toe=self.toe_cols();byp=self.byp_cols();divider=slice(nw+1,nw+self.divider_cells+1)
  for i in range(self.i0,self.nx+1):
   self.wall[i,divider]=True;self.zb[i,toe]=z_toe;self.zb[i,byp]=z_byp
   if i>self.i1:
    drop=tail_drop+tail_slope*(i-self.i1-.5)*dx;self.zb[i,toe]-=drop;self.zb[i,byp]-=drop
   else:
    self.toe[i,toe]=True;self.byp[i,byp]=True;self.ice[i,toe]=True;self.ice[i,byp]=True
  self.zb[:,0]=self.zb[:,1];self.zb[:,-1]=self.zb[:,-2];self.zb[0]=self.zb[1];self.zb[-1]=self.zb[-2]
  self.regions={};self.region_vol={}
  # PhaseSim's pre-evolution reset replaces the initial reference and records
  # initial_condition_definition. No completed evolution is erased.
  self.reset_state(eta0)
  self.geometry_contract=dict(builder_version='fixed-physical-divider-v1',lane_width_m=Wl,divider_width_m=divider_width_m,total_pool_width_m=2*Wl+divider_width_m,pool_length_m=Lpool,reach_length_m=Lr,tail_length_m=Ltail,initial_toe_z_m=z_toe,initial_bypass_z_m=z_byp,pool_floor_m=pool_floor,tail_drop_m=tail_drop,tail_slope=tail_slope,inlet_discharge_m3_s=Q,inlet_flux_per_width_m2_s=Q/(2*Wl+divider_width_m),upstream_boundary='uniform discharge per unit width across fixed full pool width',ramp_ratio=0.,coordinate_frame='idealised_cartesian_origin0_m')

 def byp_cols(self):return slice(self.nw+self.divider_cells+1,self.ny+1)
 def set_deposit(self,*args,**kwargs):raise ValueError('this fixed-geometry builder supports bare ice only; deposits require a new geometry/source contract')
 def remove_debris(self):raise ValueError('no debris is supported by this builder')
