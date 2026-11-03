"""Continuous calibrated-heat ratio enclosure over an explicit linear polytope.
Charnes-Cooper linear-fractional programming; established mathematics.
No claim of actual-field coverage or outward-rounded arithmetic.
"""
import numpy as np
from scipy.optimize import linprog

def enclosure(A,b,bounds,train,future,training_target,tolerance=0.):
 c=np.asarray(train,float);f=np.asarray(future,float);B=np.asarray(bounds,float);n=len(c)
 a=np.asarray(A,float).reshape((-1,n));b=np.asarray(b,float)
 if c.ndim!=1 or f.shape!=c.shape or B.shape!=(n,2)or b.shape!=(len(a),):raise ValueError('shape mismatch')
 if not all(np.isfinite(x).all()for x in [a,b,c,f,B])or np.any(B[:,0]>B[:,1]):raise ValueError('finite bounded polytope required')
 if not np.isfinite(training_target)or training_target<=0 or not np.isfinite(tolerance)or tolerance<0:raise ValueError('positive TRAIN target and nonnegative tolerance required')
 opts=dict(A_ub=a if len(a)else None,b_ub=b if len(a)else None,bounds=B.tolist(),method='highs')
 feasible=linprog(np.zeros(n),**opts)
 if feasible.status==2:return dict(status='infeasible_assumptions',coverage_guarantee=False)
 if not feasible.success:raise RuntimeError(feasible.message)
 denominator=linprog(c,**opts);negative=linprog(f,**opts)
 if not denominator.success or denominator.fun<=1e-12:raise ValueError('TRAIN heat not uniformly positive')
 if not negative.success or negative.fun < -1e-9:raise ValueError('future signed heat may be negative: specify phase branch before melt conversion')
 rows=[np.r_[row,-rhs]for row,rhs in zip(a,b)]
 for i,(lo,hi)in enumerate(B):
  e=np.zeros(n);e[i]=1.;rows.extend([np.r_[e,-hi],np.r_[-e,lo]])
 lp=dict(A_ub=np.asarray(rows),b_ub=np.zeros(len(rows)),A_eq=np.array([np.r_[c,0.]]),b_eq=[1.],bounds=[(None,None)]*n+[(0,None)],method='highs')
 witnesses=[]
 for sign in [1.,-1.]:
  r=linprog(np.r_[sign*f,0.],**lp)
  if not r.success or r.x[-1]<=0:raise RuntimeError('fractional enclosure failed')
  g=r.x[:-1]/r.x[-1];ratio=float(f@g/(c@g))
  violation=max(float(np.max(a@g-b))if len(a)else 0.,float(np.max(B[:,0]-g)),float(np.max(g-B[:,1])))
  if violation>1e-6:raise RuntimeError('witness fails original-polytope check')
  witnesses.append(dict(parameters=g.tolist(),ratio=ratio,calibrated_heat_mmwe=float(training_target*ratio),max_original_constraint_violation=violation))
 low,high=[w['calibrated_heat_mmwe']for w in witnesses];radius=(high-low)/2
 return dict(status='continuous_polytope_enclosure',lower_mmwe=low,upper_mmwe=high,minimax_point_mmwe=(low+high)/2,radius_mmwe=radius,abstain=bool(radius>tolerance),witnesses=witnesses,
  mathematical_scope='Exact optimization formulation for declared linear polytope and positive TRAIN denominator; floating-point LP, not outward-certified bounds. The real process need not belong to this model set.',historical_novelty=False,field_coverage_guarantee=False)

def coefficient_box_enclosure(A,b,bounds,train_lower,train_upper,future_lower,future_upper,training_target,tolerance=0.):
 """Conservative independent coefficient-box outer bound for g>=0.
 Shared sensor correlations are discarded, never asserted independent physically.
 Does not repair an invalid feasible model set or unknown nonlinear phase.
 """
 B=np.asarray(bounds,float);cl=np.asarray(train_lower,float);cu=np.asarray(train_upper,float);fl=np.asarray(future_lower,float);fu=np.asarray(future_upper,float)
 if B.ndim!=2 or np.any(B[:,0]<0):raise ValueError('nonnegative parameters required for coefficient ordering')
 if not(cl.shape==cu.shape==fl.shape==fu.shape==(len(B),))or any(not np.isfinite(x).all()for x in [cl,cu,fl,fu])or np.any(cl>cu)or np.any(fl>fu):raise ValueError('finite ordered coefficient boxes required')
 lower=enclosure(A,b,bounds,cu,fl,training_target,tolerance)
 upper=enclosure(A,b,bounds,cl,fu,training_target,tolerance)
 if lower['status']=='infeasible_assumptions':return lower
 lo=lower['lower_mmwe'];hi=upper['upper_mmwe'];radius=max(0.,(hi-lo)/2)
 return dict(status='coefficient_box_outer_enclosure',lower_mmwe=lo,upper_mmwe=hi,minimax_point_mmwe=(lo+hi)/2,radius_mmwe=radius,abstain=bool(radius>tolerance),
  lower_witness=lower['witnesses'][0],upper_witness=upper['witnesses'][1],joint_error_correlations_retained=False,field_coverage_guarantee=False,
  scope='Outer bound conditional on every true regressor/goal coefficient lying inside supplied boxes and true parameters belonging to P; exact joint feasibility and nonlinear phase unproved.')
