"""Training-only diffusivity inference. Established methods, not new physics.

No melt target enters either estimator. Data fingerprints describe supplied
training data; they cannot certify a caller has never inspected other data.
"""
import hashlib,json,math

def number(x,name,positive=False):
    if type(x) not in (float,int) or not math.isfinite(x) or (positive and x<=0):
        raise ValueError(name+' must be finite numeric'+(' and positive' if positive else ''))
    return float(x)

def prepare(times,observations,sensor_ids,boundary_sensor_ids,train_end,evaluation_start):
    t=[number(x,'time') for x in times]
    end=number(train_end,'train end'); start=number(evaluation_start,'evaluation start')
    if len(t)<2 or any(b<=a for a,b in zip(t,t[1:])):raise ValueError('strictly increasing training times required')
    if t[-1]>end or start<=end:raise ValueError('training must end before evaluation; no supplied evaluation rows')
    ids=list(sensor_ids)
    if not ids or len(set(ids))!=len(ids) or any(not isinstance(x,str) or not x for x in ids):raise ValueError('unique sensor labels required')
    if set(ids)&set(boundary_sensor_ids):raise ValueError('imposed boundary sensors cannot be scored in fitting')
    y=[[number(x,'observed temperature') for x in row] for row in observations]
    if len(y)!=len(t) or any(len(row)!=len(ids) for row in y):raise ValueError('temperature shape mismatch')
    source=dict(times_s=t,observations_C=y,sensor_ids=ids,train_end_s=end,evaluation_start_s=start)
    fingerprint=hashlib.sha256(json.dumps(source,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
    return t,y,ids,fingerprint

def fit_grid(times_s,observations_C,*,sensor_ids,boundary_sensor_ids,train_end_s,
             evaluation_start_s,alpha_grid_m2_s,predict_training,sensor_absolute_bound_C=.4):
    """Select alpha by squared distance to the declared observation intervals.

    Callback accepts alpha and returns predictions exactly at supplied training
    times/sensors. It must simulate using training forcing/profile data only.
    The finite grid/selection rule is supplied before evaluation scoring.
    Bound is conditional sensor input, not a Gaussian sigma or field CI.
    """
    t,y,ids,fingerprint=prepare(times_s,observations_C,sensor_ids,boundary_sensor_ids,train_end_s,evaluation_start_s)
    delta=number(sensor_absolute_bound_C,'sensor bound')
    if delta<0:raise ValueError('sensor bound must be nonnegative')
    grid=[number(a,'diffusivity',True) for a in alpha_grid_m2_s]
    if not grid or any(b<=a for a,b in zip(grid,grid[1:])):raise ValueError('strictly increasing explicit diffusivity grid required')
    results=[]
    for alpha in grid:
        prediction=[[number(x,'predicted temperature') for x in row] for row in predict_training(alpha)]
        if len(prediction)!=len(t) or any(len(row)!=len(ids) for row in prediction):raise ValueError('training callback shape mismatch')
        residuals=[p-o for pp,oo in zip(prediction,y) for p,o in zip(pp,oo)]
        loss=math.fsum(max(abs(r)-delta,0.)**2 for r in residuals)/len(residuals)
        results.append(dict(alpha_m2_s=alpha,interval_squared_distance_C2=loss,
                            nominal_rmse_C=math.sqrt(math.fsum(r*r for r in residuals)/len(residuals)),
                            maximum_absolute_residual_C=max(abs(r) for r in residuals)))
    # Exact ties are resolved by smallest diffusivity, not by evaluation error.
    winner=min(results,key=lambda x:(x['interval_squared_distance_C2'],x['alpha_m2_s']))
    best=winner['interval_squared_distance_C2']
    return dict(selected_alpha_m2_s=winner['alpha_m2_s'],selection='minimum training interval loss; exact ties smallest alpha',
                selected_at_grid_boundary=winner is results[0] or winner is results[-1],
                loss_identifies_unique_grid_point=sum(x['interval_squared_distance_C2']==best for x in results)==1,
                fully_interval_compatible_grid_m2_s=[x['alpha_m2_s'] for x in results if x['interval_squared_distance_C2']==0],
                grid_results=results,training_sha256=fingerprint,training_rows=len(t),fitted_sensor_ids=ids,
                training_end_s=train_end_s,evaluation_start_s=evaluation_start_s,
                sensor_absolute_bound_C=delta,thermal_parameter_fit_uses_melt=False,
                uncertainty_scope='conditional sensor intervals only; boundary, geometry and model errors additional',
                statistical_confidence_interval=False)

def derivative_crh(times_s,profiles_C,depths_m,*,train_end_s,evaluation_start_s,min_abs_curvature_K_m2=0.):
    """Zero-intercept LS T_t=alpha*T_zz, using interior time/depth nodes.

    Three-point quadratic formulas support unequal spacings. Temperature
    derivatives amplify noise; this comparator is not a robust causal estimate.
    Boundary nodes participate as neighbors, never as fitted residual rows.
    """
    t=[number(x,'time') for x in times_s];z=[number(x,'depth') for x in depths_m]
    if len(t)<3 or len(z)<3 or any(b<=a for a,b in zip(t,t[1:])) or any(b<=a for a,b in zip(z,z[1:])):
        raise ValueError('three increasing time/depth nodes required')
    if t[-1]>number(train_end_s,'train end') or number(evaluation_start_s,'evaluation start')<=train_end_s:
        raise ValueError('training/evaluation ordering invalid')
    y=[[number(x,'temperature') for x in row] for row in profiles_C]
    if len(y)!=len(t) or any(len(row)!=len(z) for row in y):raise ValueError('profile shape mismatch')
    threshold=number(min_abs_curvature_K_m2,'curvature threshold')
    if threshold<0:raise ValueError('negative curvature threshold')
    xvals=[];yvals=[]
    for i in range(1,len(t)-1):
        h0=t[i]-t[i-1];h1=t[i+1]-t[i]
        for j in range(1,len(z)-1):
            dz0=z[j]-z[j-1];dz1=z[j+1]-z[j]
            Tt=-h1*y[i-1][j]/(h0*(h0+h1))+(h1-h0)*y[i][j]/(h0*h1)+h0*y[i+1][j]/(h1*(h0+h1))
            Tzz=2*(y[i][j-1]/(dz0*(dz0+dz1))-y[i][j]/(dz0*dz1)+y[i][j+1]/(dz1*(dz0+dz1)))
            if abs(Tzz)>threshold:xvals.append(Tzz);yvals.append(Tt)
    den=math.fsum(x*x for x in xvals)
    if den==0:raise ValueError('diffusivity unidentifiable from zero retained curvature')
    alpha=math.fsum(x*y for x,y in zip(xvals,yvals))/den
    if not math.isfinite(alpha) or alpha<=0:raise ValueError('nonphysical inferred diffusivity; do not clip')
    return dict(alpha_m2_s=alpha,method='zero-intercept derivative regression on training-only interior nodes',
                retained_rows=len(xvals),curvature_sum_squares=den,
                residual_rmse_K_s=math.sqrt(math.fsum((b-alpha*a)**2 for a,b in zip(xvals,yvals))/len(xvals)),
                uses_melt_target=False,thermal_conductivity_identified=False,
                explanation='Temperature identifies diffusivity k/C; bulk heat capacity remains an independent input.')
