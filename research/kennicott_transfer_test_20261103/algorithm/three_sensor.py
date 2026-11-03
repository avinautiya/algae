"""Three-sensor integrated thermal inference, no melt targets in fitting.

Established P1 finite-element energy/NNLS mathematics; not a novelty claim.
No basal-gap extrapolation is hidden. Absolute capacity remains external.
"""
import bisect,math

class UnsupportedModel(ValueError):pass

def validate(z,t,T,maximum_gap_seconds=3600.):
 if len(z)!=3 or not all(math.isfinite(v)for v in z)or not z[0]<z[1]<z[2]:raise ValueError('exactly3 ordered finite sensor depths required')
 if len(t)!=len(T)or len(t)<3 or any(not math.isfinite(v)for v in t)or any(b<=a for a,b in zip(t[:-1],t[1:])):raise ValueError('strictly ordered, finite native clock required')
 if any(len(row)!=3 or any(not math.isfinite(v)for v in row)for row in T):raise ValueError('complete finite3sensor rows required; no imputation')
 if maximum_gap_seconds is not None and(maximum_gap_seconds<=0 or any(b-a>maximum_gap_seconds for a,b in zip(t[:-1],t[1:]))):raise UnsupportedModel('native gap exceeds declared limit')

def integrated_rows(z,t,T,balance_window_seconds,maximum_gap_seconds=3600.):
 validate(z,t,T,maximum_gap_seconds)
 if not math.isfinite(balance_window_seconds)or balance_window_seconds<=0:raise ValueError('positive balance window')
 hl,hr=z[1]-z[0],z[2]-z[1];mass=[hl/6,(hl+hr)/3,hr/6];A=[];y=[];meta=[];a=0
 while a<len(t)-1:
  b=min(bisect.bisect_left(t,t[a]+balance_window_seconds),len(t)-1);duration=t[b]-t[a]
  integ=[math.fsum((t[j+1]-t[j])*(T[j][i]+T[j+1][i])/2 for j in range(a,b))for i in range(3)]
  A.append([(integ[0]-integ[1])/hl/duration,(integ[2]-integ[1])/hr/duration]);y.append(sum(m*(T[b][i]-T[a][i])for i,m in enumerate(mass))/duration)
  meta.append({'first_index':a,'last_index':b,'duration_seconds':duration});a=b
 return A,y,meta

def nnls2(A,y):
 """Exact two-variable NNLS candidates: interior, two faces, origin."""
 if len(A)!=len(y)or not A:raise ValueError('nonempty design')
 p=math.fsum(x[0]*x[0]for x in A);q=math.fsum(x[0]*x[1]for x in A);r=math.fsum(x[1]*x[1]for x in A);u=math.fsum(x[0]*v for x,v in zip(A,y));v=math.fsum(x[1]*v for x,v in zip(A,y));det=p*r-q*q
 candidates=[(0.,0.)]
 if p>0:candidates.append((max(u/p,0.),0.))
 if r>0:candidates.append((0.,max(v/r,0.)))
 if det>1e-15*p*r:
  c=((u*r-v*q)/det,(v*p-u*q)/det)
  if min(c)>=0:candidates.append(c)
 def loss(c):return math.fsum((a*c[0]+b*c[1]-v)**2 for(a,b),v in zip(A,y))
 best=min(candidates,key=lambda c:(loss(c),sum(c),c));trace=p+r;disc=math.hypot(p-r,2*q);hi=max((trace+disc)/2,0.);lo=max((trace-disc)/2,0.)
 rho=q/math.sqrt(p*r)if p*r>0 else None
 return best,{'singular_values':[math.sqrt(hi),math.sqrt(lo)],'singular_ratio':math.sqrt(lo/hi)if hi>0 else 0.,'column_correlation':rho,'squared_residual':loss(best),'candidate_count':len(candidates)}

def fit(z,t,T,*,balance_window_seconds=21600.,minimum_singular_ratio=1e-3,maximum_gap_seconds=3600.,sensor_absolute_bound_C=None):
 A,y,meta=integrated_rows(z,t,T,balance_window_seconds,maximum_gap_seconds);alpha,diag=nnls2(A,y)
 residual=[sum(a*b for a,b in zip(row,alpha))-v for row,v in zip(A,y)];supported=diag['singular_ratio']>=minimum_singular_ratio and min(alpha)>0
 report={'alpha_m2_s':list(alpha),'diagnostics':diag,'balance_RMSE_K_m_per_s':math.sqrt(sum(r*r for r in residual)/len(residual)),'balance_rows':meta,'positive_parameter_fit':min(alpha)>0,'practically_identifiable_under_declared_threshold':diag['singular_ratio']>=minimum_singular_ratio,'supported_point_estimate':supported,'absolute_capacity_identified':False,'melt_target_used':False,'parameter_uncertainty_interval':None,'historical_mathematical_novelty_claim':False}
 if sensor_absolute_bound_C is not None:
  if not math.isfinite(sensor_absolute_bound_C)or sensor_absolute_bound_C<0:raise ValueError('finite nonnegative sensor bound required')
  hl,hr=z[1]-z[0],z[2]-z[1];d=sensor_absolute_bound_C
  radii=[d*(hl+hr)/m['duration_seconds']+2*d*(alpha[0]/hl+alpha[1]/hr)for m in meta]
  report['necessary_sensor_outer_bound_consistent_rows']=sum(abs(x)<=r for x,r in zip(residual,radii));report['necessary_sensor_outer_bound_total_rows']=len(meta)
  report['joint_coherent_sensor_error_feasibility_established']=False
 return report

def slab_moment_K_m(Tmiddle,Tdeep,spacing_m,gap_m,basal_T_C):
 """Integral of P1 temperature from deepest-edge midpoint to ice."""
 if not all(math.isfinite(v)for v in[Tmiddle,Tdeep,spacing_m,gap_m,basal_T_C])or spacing_m<=0 or gap_m<0:raise ValueError('finite declared gap and temperatures')
 return spacing_m*Tmiddle/8+3*spacing_m*Tdeep/8+gap_m*(Tdeep+basal_T_C)/2

def reconstruct(z,t,T,*,ice_depth_m,basal_temperature_C,balance_window_seconds=21600.,minimum_singular_ratio=1e-3,maximum_gap_seconds=3600.,sensor_absolute_bound_C=None):
 """Full-block temperature-conditioned reconstruction, not future-weather forecast.

 alpha inferred from all temperatures in this block; no camera used. Signed
 flux across measured deepest-edge midpoint minus basal-slab sensible storage.
 Constant C makes Q/C dimensionK*m and allows an honest output-scale quotient.
 """
 if not math.isfinite(ice_depth_m)or ice_depth_m<z[2]:raise UnsupportedModel('declared ice depth below deepest sensor required')
 if not math.isfinite(basal_temperature_C):raise ValueError('explicit finite basal temperature required')
 f=fit(z,t,T,balance_window_seconds=balance_window_seconds,minimum_singular_ratio=minimum_singular_ratio,maximum_gap_seconds=maximum_gap_seconds,sensor_absolute_bound_C=sensor_absolute_bound_C)
 if not f['supported_point_estimate']:raise UnsupportedModel('rank/excitation/positive-parameter gate failed')
 gap=ice_depth_m-z[2];hr=z[2]-z[1];G=f['alpha_m2_s'][1]/hr;trace=[]
 for j in range(len(t)-1):
  duration=t[j+1]-t[j];incoming=G*duration*((T[j][1]-T[j][2])+(T[j+1][1]-T[j+1][2]))/2
  dH=slab_moment_K_m(T[j+1][1],T[j+1][2],hr,gap,basal_temperature_C)-slab_moment_K_m(T[j][1],T[j][2],hr,gap,basal_temperature_C)
  trace.append({'start_seconds':t[j],'end_seconds':t[j+1],'normalized_top_heat_K_m':incoming,'normalized_debris_storage_change_K_m':dH,'normalized_interface_heat_K_m':incoming-dH})
 return {'fit':f,'normalized_heat_trace':trace,'normalized_total_interface_heat_K_m':math.fsum(x['normalized_interface_heat_K_m']for x in trace),'negative_interface_heat_steps':sum(x['normalized_interface_heat_K_m']<0 for x in trace),'basal_gap_m':gap,'basal_temperature_C':basal_temperature_C,'initial_ice_cold_content_identified':False,'ice_cold_content_evolution_solved':False,'camera_or_melt_target_used':False,'temperature_fit_uses_full_block_covariates':True,'future_weather_forecast':False,'no_silent_energy_clipping':True}

def aggregate(trace,start_seconds,end_seconds):
 """Exact trace-edge support only: no invented daily/14-day target endpoints."""
 if not start_seconds<end_seconds:raise ValueError('ordered support')
 selected=[r for r in trace if r['start_seconds']>=start_seconds and r['end_seconds']<=end_seconds]
 if not selected or selected[0]['start_seconds']!=start_seconds or selected[-1]['end_seconds']!=end_seconds or any(a['end_seconds']!=b['start_seconds']for a,b in zip(selected[:-1],selected[1:])):raise UnsupportedModel('target endpoints must bind exact complete native trace edges')
 return math.fsum(r['normalized_interface_heat_K_m']for r in selected)

def scale_quotient(train_observed_mmwe,train_normalized_heat,evaluation_normalized_heat):
 if not all(math.isfinite(v)for v in[train_observed_mmwe,train_normalized_heat,evaluation_normalized_heat])or train_observed_mmwe<=0 or train_normalized_heat<=0:raise UnsupportedModel('positive independent TRAIN target and normalized heat required')
 if evaluation_normalized_heat<0:raise UnsupportedModel('negative signed heat is not unsupported melt silently clipped')
 return train_observed_mmwe*evaluation_normalized_heat/train_normalized_heat

def absolute_melt_mmwe(model,capacity_J_m3K,*,temperate_ice_confirmed,all_other_interface_terms_zero_confirmed,latent_J_kg=333700.):
 if temperate_ice_confirmed is not True or all_other_interface_terms_zero_confirmed is not True or model['basal_temperature_C']!=0:raise UnsupportedModel('independent temperate interface/complete energy balance required')
 if model['negative_interface_heat_steps']:raise UnsupportedModel('negative interface flux requires solved cold content or refreezing')
 if not math.isfinite(capacity_J_m3K)or capacity_J_m3K<=0 or not math.isfinite(latent_J_kg)or latent_J_kg<=0:raise ValueError('external positive capacity and latent heat required')
 return capacity_J_m3K*model['normalized_total_interface_heat_K_m']/latent_J_kg

def fit_past(z,t,T,*,asof_seconds,lookback_seconds,minimum_history_seconds,balance_window_seconds=21600.,minimum_singular_ratio=1e-3,maximum_gap_seconds=3600.,sensor_absolute_bound_C=None):
 """Parameter adaptation sees only supplied clock labels<=asof, never future."""
 if len(t)!=len(T)or len(t)<3 or any(b<=a for a,b in zip(t[:-1],t[1:])):raise ValueError('ordered source clock required')
 if not all(math.isfinite(v)and v>0 for v in[lookback_seconds,minimum_history_seconds])or minimum_history_seconds>lookback_seconds:raise ValueError('explicit positive history<=lookback required')
 stop=bisect.bisect_right(t,asof_seconds);start=bisect.bisect_left(t,asof_seconds-lookback_seconds)
 if stop-start<3 or t[stop-1]-t[start]<minimum_history_seconds:raise UnsupportedModel('insufficient past thermal history; no invented warmup flux')
 r=fit(z,t[start:stop],T[start:stop],balance_window_seconds=balance_window_seconds,minimum_singular_ratio=minimum_singular_ratio,maximum_gap_seconds=maximum_gap_seconds,sensor_absolute_bound_C=sensor_absolute_bound_C)
 r.update({'asof_seconds':asof_seconds,'fit_first_source_seconds':t[start],'fit_last_source_seconds':t[stop-1],'future_covariates_used_in_parameter_fit':False})
 return r
