"""Classical passive-network absorption bound for uncertain initial heat.
Zero perturbation at both Dirichlet boundaries, positive fixed conductances/C.
No nonlinear/unknown-phase or moving-geometry guarantee.
"""
import numpy as np

def right_absorption_weights(edge_conductances):
 g=np.asarray(edge_conductances,float)
 if g.ndim!=1 or len(g)<2 or not np.isfinite(g).all()or np.any(g<=0):raise ValueError('positive finite edges required')
 r=1/g
 return np.cumsum(r)[:-1]/sum(r)

def initial_heat_bound(edge_conductances,node_capacities,error_radii,latent_heat=334000.):
 p=right_absorption_weights(edge_conductances);c=np.asarray(node_capacities,float);e=np.asarray(error_radii,float)
 if c.shape!=p.shape or e.shape!=p.shape or not np.isfinite(c).all()or not np.isfinite(e).all()or np.any(c<=0)or np.any(e<0)or not np.isfinite(latent_heat)or latent_heat<=0:raise ValueError('finite capacities and temperature-error radii required')
 q=float(np.dot(c*p,e))
 return dict(initial_state_right_heat_error_bound_J_m2=q,melt_equivalent_bound_mmwe=q/latent_heat,right_absorption_weights=p.tolist(),
  scope='All finite horizons, conditional on fixed passive linear graph and identical boundary perturbations0; no field-confidence or nonlinear-phase claim.')
