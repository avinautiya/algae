"""Positive thermal-edge inference by integrated consistent-P1 energy balance.

Established inverse heat/RC-network mathematics. Source temperatures are noisy
regressors, so least-squares curvature is not a parameter confidence interval.
"""
import hashlib,json,math
import numpy as np
from scipy.optimize import lsq_linear

def fit(depths_m,times_s,temperatures_C,*,capacity_J_m3K,window_s=86400.,
        conductivity_bounds_W_mK=(.05,6.),sensor_bound_C=.4):
    z=np.asarray(depths_m,dtype=float);t=np.asarray(times_s,dtype=float);T=np.asarray(temperatures_C,dtype=float)
    if z.ndim!=1 or len(z)<3 or np.any(np.diff(z)<=0)or t.ndim!=1 or len(t)<3 or np.any(np.diff(t)<=0):raise ValueError('ordered depth/time coordinates required')
    if T.shape!=(len(t),len(z)) or not np.isfinite(T).all()or not np.isfinite(z).all()or not np.isfinite(t).all():raise ValueError('finite shape-matched profiles required')
    if not math.isfinite(capacity_J_m3K) or capacity_J_m3K<=0 or not math.isfinite(window_s)or window_s<=0:raise ValueError('positive capacity/window required')
    kl,ku=conductivity_bounds_W_mK
    if not all(math.isfinite(x) for x in (kl,ku,sensor_bound_C)) or not 0<kl<ku or sensor_bound_C<0:raise ValueError('positive physical bounds required')
    A=[];y=[];records=[];first=0
    while first<len(t)-1:
        end=min(int(np.searchsorted(t,t[first]+window_s)),len(t)-1)
        duration=t[end]-t[first]
        integrals=np.sum((T[first:end]+T[first+1:end+1])*.5*np.diff(t[first:end+1])[:,None],axis=0)
        change=T[end]-T[first]
        for j in range(1,len(z)-1):
            hL=z[j]-z[j-1];hR=z[j+1]-z[j]
            rhs=capacity_J_m3K*(hL*change[j-1]/6+(hL+hR)*change[j]/3+hR*change[j+1]/6)
            row=np.zeros(len(z)-1);row[j-1]=integrals[j-1]-integrals[j];row[j]=-(integrals[j]-integrals[j+1])
            # Normalize each integrated balance to W/m², retaining partial-window support.
            A.append(row/duration);y.append(rhs/duration)
            records.append(dict(first=first,end=end,sensor=j,duration_s=float(duration),hL=float(hL),hR=float(hR)))
        first=end
    A=np.asarray(A);y=np.asarray(y);h=np.diff(z)
    solution=lsq_linear(A,y,bounds=(kl/h,ku/h),tol=1e-12,max_iter=1000)
    if not solution.success:raise RuntimeError('bounded energy inversion did not converge')
    G=solution.x;residual=A@G-y;bounds=[]
    for r in records:
        j=r['sensor'];duration=r['duration_s']
        bounds.append(capacity_J_m3K*sensor_bound_C*(r['hL']+r['hR'])/duration+2*sensor_bound_C*(G[j-1]+G[j]))
    s=np.linalg.svd(A,compute_uv=False)
    fingerprint=hashlib.sha256(json.dumps(dict(z=z.tolist(),t=t.tolist(),T=T.tolist(),C=capacity_J_m3K,window=window_s),separators=(',',':')).encode()).hexdigest()
    return dict(conductances_W_m2K=G.tolist(),conductivities_W_mK=(G*h).tolist(),
        conductivity_bounds_W_mK=list(conductivity_bounds_W_mK),capacity_J_m3K=capacity_J_m3K,
        interval_seconds=window_s,training_sha256=fingerprint,balance_RMSE_W_m2=float(np.sqrt(np.mean(residual**2))),
        max_absolute_balance_residual_W_m2=float(np.max(np.abs(residual))),
        conservative_regressor_sensor_bounds_W_m2=bounds,
        balances_within_conditional_sensor_bound=int(np.sum(np.abs(residual)<=np.asarray(bounds))),
        total_balances=len(y),singular_values=s.tolist(),rank=int(np.linalg.matrix_rank(A)),
        unique_linear_solution_conditional_on_exact_inputs=bool(np.linalg.matrix_rank(A)==len(G)),
        physical_bound_active=[bool(abs(g*x-kl)<1e-6 or abs(g*x-ku)<1e-6)for g,x in zip(G,h)],
        C_assumed_not_identified=True,errors_in_regressors=True,confidence_interval=False,
        historical_mathematical_novelty_claim=False,used_melt_target=False)

def lowest_edge_enclosure(depths_m,times_s,temperatures_C,*,capacity_J_m3K,
                          window_s=86400.,conductivity_bounds_W_mK=(.05,6.),sensor_bound_C=.4):
    """Necessary LP relaxation, not exact feasible errors-in-variables inference.

    Keeps sensor-error contributions to BOTH sides of consistent-P1 balances.
    All statements are conditional on piecewise-linear time reconstruction,
    fixed geometry/C, constant edge conductances and no unmodelled source.
    The returned heat is across the deepest MEASURED span, not ice melt.
    """
    from scipy.optimize import linprog
    # Reuse full validation, but the NNLS estimate never restricts the LP.
    diagnostic=fit(depths_m,times_s,temperatures_C,capacity_J_m3K=capacity_J_m3K,
        window_s=window_s,conductivity_bounds_W_mK=conductivity_bounds_W_mK,sensor_bound_C=sensor_bound_C)
    z=np.asarray(depths_m,float);t=np.asarray(times_s,float);T=np.asarray(temperatures_C,float)
    A=[];b=[];first=0;edges=len(z)-1
    while first<len(t)-1:
        end=min(int(np.searchsorted(t,t[first]+window_s)),len(t)-1);duration=t[end]-t[first]
        integrals=np.sum((T[first:end]+T[first+1:end+1])*.5*np.diff(t[first:end+1])[:,None],axis=0)
        change=T[end]-T[first]
        for j in range(1,len(z)-1):
            hL=z[j]-z[j-1];hR=z[j+1]-z[j]
            y=capacity_J_m3K*(hL*change[j-1]/6+(hL+hR)*change[j]/3+hR*change[j+1]/6)/duration
            row=np.zeros(edges);row[j-1]=(integrals[j-1]-integrals[j])/duration;row[j]=-(integrals[j]-integrals[j+1])/duration
            reg=np.zeros(edges);reg[j-1]=reg[j]=2*sensor_bound_C
            radius=capacity_J_m3K*sensor_bound_C*(hL+hR)/duration
            A.extend([row-reg,-row-reg]);b.extend([y+radius,-y+radius])
        first=end
    kl,ku=conductivity_bounds_W_mK;bounds=[(kl/h,ku/h)for h in np.diff(z)]
    objective=np.zeros(edges);objective[-1]=1.
    low=linprog(objective,A_ub=A,b_ub=b,bounds=bounds,method='highs')
    high=linprog(-objective,A_ub=A,b_ub=b,bounds=bounds,method='highs')
    if low.status==2 or high.status==2:return dict(status='necessary_model_constraints_infeasible',icemelt_prediction=False)
    if not low.success or not high.success:raise RuntimeError('LP did not produce verified finite extrema')
    return dict(status='necessary_relaxation_bounded',lowest_conductance_W_m2K=[float(low.fun),float(-high.fun)],
        lowest_conductivity_W_mK=[float(low.fun*(z[-1]-z[-2])),float(-high.fun*(z[-1]-z[-2]))],
        NNLS_diagnostic=diagnostic,icemelt_prediction=False,
        exact_joint_sensor_error_feasibility_established=False,confidence_interval=False,
        assumptions='Fixed geometry/C, piecewise-linear samples, positive constant edges, zero unmodelled heat; necessary outer relaxation only')
