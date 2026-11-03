"""Post-test explanatory FOSD falsifier; no refit, no selected beta, no new targets."""
import argparse,bisect,datetime as dt,json,math
from pathlib import Path
from postscoreaudit import rows,read,native_air_area

def distribution(values,weights):
 if len(values)!=len(weights)or not values or any(not math.isfinite(v)for v in values+weights)or any(w<0 for w in weights)or sum(weights)<=0:
  raise ValueError('finite nonnegative nonzero weights required')
 d={}
 for x,w in zip(values,weights):
  if w>0:d[x]=d.get(x,0.)+w
 return sorted((x,w/sum(weights))for x,w in d.items())

def fosd(test,train,tolerance=1e-12):
 """Test dominates TRAIN iff F_test(x)<=F_train(x) at every support atom."""
 thresholds=sorted({x for x,w in test}|{x for x,w in train});a=dict(test);b=dict(train);ct=cr=0.;checks=[]
 for x in thresholds:
  ct+=a.get(x,0);cr+=b.get(x,0);checks.append((x,ct-cr))
 worst=max(checks,key=lambda t:t[1]);strict=min(checks,key=lambda t:t[1])
 return {'test_first_order_stochastically_dominates_train':worst[1]<=tolerance,'maximum_CDF_violation':worst[1],'maximum_violation_threshold':worst[0],'minimum_CDF_difference':strict[1],'checked_support_thresholds':len(thresholds),'numerical_tolerance':tolerance}

def increasing_convex_order(test,train,tolerance=1e-12):
 """Finite-support stop-loss characterization of increasing convex order.

 Tail differences are piecewise linear between union knots, constant below
 the lowest knot, and zero above the highest; checking knots suffices.
 """
 thresholds=sorted({x for x,w in test}|{x for x,w in train});checks=[]
 for t in thresholds:
  gap=math.fsum(w*max(x-t,0)for x,w in test)-math.fsum(w*max(x-t,0)for x,w in train)
  checks.append((t,gap))
 worst=min(checks,key=lambda x:x[1]);nontrivial=checks[:-1]
 return {'test_increasing_convex_dominates_train':worst[1]>=-tolerance,'minimum_stop_loss_difference':worst[1],'minimum_threshold':worst[0],'minimum_below_maximum_support':min((x[1]for x in nontrivial),default=0.),'checked_support_thresholds':len(thresholds),'numerical_tolerance':tolerance}

def icx_moisture_error_radius(test,train):
 """Worst pointwise covariate errors: shift test down, TRAIN up.

 For increasing convex objectives those simultaneous extremal shifts are
 worst-case. Does not assume an independent probability distribution.
 """
 def gate(delta):return increasing_convex_order([(x-delta,w)for x,w in test],[(x+delta,w)for x,w in train],tolerance=1e-14)['test_increasing_convex_dominates_train']
 if not gate(0.):return {'maximum_symmetric_pointwise_theta_error_radius':0.,'status':'nominal_gate_failed'}
 lo=0.;hi=.25
 if gate(hi):raise ValueError('search did not bracket failure')
 for _ in range(60):
  middle=(lo+hi)/2
  if gate(middle):lo=middle
  else:hi=middle
 return {'maximum_symmetric_pointwise_theta_error_radius':lo,'numerical_bracket':[lo,hi],'stop_loss_numeric_tolerance':1e-14,'gate_at_typical_sensor_accuracy_point031':gate(.031),'condition':'Observed air weights fixed; all test moisture may shift down by delta and all TRAIN moisture up by delta; no clipping; covers arbitrary pointwise errors of this size. Common constant offset cancels.'}

def expfactor(test,train,beta):
 # A common reference cancels, stable for moderate declared beta and theta.
 ref=min(x for x,w in train+test)
 return math.fsum(w*math.exp(beta*(x-ref))for x,w in test)/math.fsum(w*math.exp(beta*(x-ref))for x,w in train)

def run(root,directory):
 root=Path(root);directory=Path(directory);p=read(directory/'PREDICTIONS.json');r=read(directory/'RESULTS.json');n=read(directory/'NORMALIZED_PREDICTIONS.json')
 features=rows(root/'data/moisture/Tower2_CANONICAL_FEATURES.csv');labels=[dt.datetime.fromisoformat(x['SOURCE_LOCAL_LABEL'])for x in features]
 air=[(dt.datetime.fromisoformat(x['SOURCE_LOCAL_LABEL']),float(x['TA_3']))for x in rows(root/'data/moisture/Tower2_NATIVE_AIR_TA3.csv')]
 theta=[float(x['SM3_raw_m3_m3'])for x in features[:-1]];weights=[native_air_area(a,b,air)for a,b in zip(labels[:-1],labels[1:])]
 endpoints=[labels.index(dt.datetime.fromisoformat(x))for x in n['epochs']]
 dist=[distribution(theta[a:b],weights[a:b])for a,b in zip(endpoints[:-1],endpoints[1:])];train=dist[0];periods=[]
 for j,test in enumerate(dist[1:]):
  gate=fosd(test,train);icx=increasing_convex_order(test,train);factor=expfactor(test,train,5.4);air_pred=p['models']['air_pdd'][j]['predicted_mmwe'];primary_pred=p['models']['anomaly_exp_air'][j]['predicted_mmwe']
  if not math.isclose(primary_pred,air_pred*factor,rel_tol=1e-12):raise AssertionError('fixed beta primary decomposition failed')
  base_error=r['model_metrics']['air_pdd']['period_errors'][j]['signed_error_mmwe']
  grid=[{'beta':i/10,'exponential_multiplier_relative_to_air_PDD':expfactor(test,train,i/10)}for i in range(401)]
  periods.append({'start':n['epochs'][j+1],'end':n['epochs'][j+2],**gate,'increasing_convex_order':icx,'ICX_conditional_moisture_uncertainty_gate':icx_moisture_error_radius(test,train),'ICX_and_positive_PDD_error_rules_out_improvement_for_all_beta_nonnegative':icx['test_increasing_convex_dominates_train']and base_error>=0,'PDD_signed_observed_error_mmwe':base_error,'air_exposure_weighted_mean_theta':sum(x*w for x,w in test),'train_air_exposure_weighted_mean_theta':sum(x*w for x,w in train),'locked_beta':5.4,'locked_primary_multiplier':factor,'locked_decomposition_residual_mmwe':primary_pred-air_pred*factor,'FOSD_and_positive_PDD_error_rules_out_improvement_for_all_beta_nonnegative':gate['test_first_order_stochastically_dominates_train']and base_error>=0,'descriptive_beta_grid_no_selection':grid,'grid_min_multiplier':min(x['exponential_multiplier_relative_to_air_PDD']for x in grid),'grid_min_positive_beta_multiplier':min(x['exponential_multiplier_relative_to_air_PDD']for x in grid[1:]),'no_beta_fit_performed':True})
 return {'scope':'Post-test mechanistic/statistical explanation of issued models, not new validation or new historical theorem','periods':periods,'all_periods_FOSD':all(x['test_first_order_stochastically_dominates_train']for x in periods),'minimum_all_periods_ICX_moisture_error_radius':min(x['ICX_conditional_moisture_uncertainty_gate']['maximum_symmetric_pointwise_theta_error_radius']for x in periods),'conditioned_on_observed_covariates_not_certified_true_moisture':True,'all_periods_increasing_convex_order':all(x['increasing_convex_order']['test_increasing_convex_dominates_train']for x in periods),'all_periods_positive_PDD_error':all(x['PDD_signed_observed_error_mmwe']>=0 for x in periods),'new_target_access':False,'parameter_selection':False,'causal_moisture_effect_established':False}

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',required=True);p.add_argument('--directory',required=True);p.add_argument('--output',required=True);a=p.parse_args();r=run(a.root,a.directory)
 with Path(a.output).open('x')as f:json.dump(r,f,indent=2,allow_nan=False);f.write('\n')
 print(json.dumps({'all_periods_FOSD':r['all_periods_FOSD'],'periods':[{k:v for k,v in x.items()if k!='descriptive_beta_grid_no_selection'}for x in r['periods']]},indent=2))
