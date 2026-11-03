"""Retrospective bounded-error loss audit; no fitting or deployment selection.

true cumulative readings = recorded readings + e_j. Shared middle reading must
enter training and evaluation increments with opposite signs. Worst-case losses
include this dependence exactly. Source uncertainty interpretation is unknown.
"""
import argparse,hashlib,itertools,json,math
from pathlib import Path
MODELS=('constant_rate','positive_degree_day','net_shortwave_dose')
ASSOCIATIONS=('per_reading','per_interval')
def parameters(x,y,f,delta,association):
 if not all(type(v) in (int,float) and math.isfinite(v) and v>=0 for v in (x,y,f,delta)):raise ValueError('finite nonnegative inputs required')
 if association not in ASSOCIATIONS:raise ValueError('explicit uncertainty semantics required')
 r=2*delta if association=='per_reading' else delta
 if x<=r or y<r:raise ValueError('physical nonnegative increments not established over full uncertainty box')
 return (2*f+2)*delta if association=='per_reading' else (f+1)*delta

def worst_error(x,y,f,delta,association):
 r=parameters(x,y,f,delta,association);nominal=f*x-y
 return dict(nominal_signed_error_mmwe=nominal,signed_error_bounds_mmwe=[nominal-r,nominal+r],worst_absolute_error_mmwe=abs(nominal)+r)

def dominance(x,y,fa,fb,delta,association):
 parameters(x,y,fa,delta,association);parameters(x,y,fb,delta,association)
 # Ra²-Rb²=(fa-fb)*true_X*((fa+fb)*true_X-2*true_Y).
 # true_X is positive over the full box. Remaining factor is affine.
 s=fa+fb;center=s*x-2*y
 radius=(2*s+4)*delta if association=='per_reading' else (s+2)*delta
 lo,hi=center-radius,center+radius
 vals=[(fa-fb)*lo,(fa-fb)*hi];lower,upper=min(vals),max(vals)
 if fa==fb:status='identical_predictions'
 elif upper<0:status='a_strictly_better_everywhere'
 elif lower>0:status='b_strictly_better_everywhere'
 elif upper<=0:status='a_no_worse_with_possible_ties'
 elif lower>=0:status='b_no_worse_with_possible_ties'
 else:status='ranking_changes_with_admissible_measurement_error'
 return dict(model_a_minus_b_squared_loss_sign_factor_bounds=[lower,upper],status=status)

def fixed_forecast(x,y,f,delta,association):
 parameters(x,y,f,delta,association)
 target_radius=(2 if association=='per_reading' else 1)*delta
 residual=f*x-y
 return dict(prediction_mmwe=f*x,worst_absolute_error_mmwe=abs(residual)+target_radius,
  signed_error_bounds_mmwe=[residual-target_radius,residual+target_radius])

def run(source,out):
 if out.exists():raise ValueError('refuse overwriting a previous audit')
 raw=source.read_bytes();data=json.loads(raw);groups={}
 for r in data['rows']:
  if r['stake']=='H1D' and r['model'] in MODELS:
   key=(r['mean_label_end_offset_h'],tuple(r['visit_hours']))
   if r['model'] in groups.setdefault(key,{}):raise ValueError('duplicate fixed support')
   groups[key][r['model']]=r
 expected={(offset,hours) for offset in (0.,.5,1.) for hours in itertools.product((0.,12.,24.),repeat=3)}
 if set(groups)!=expected or any(set(g)!=set(MODELS) for g in groups.values()):raise ValueError('all81 fixed timing scenarios required')
 rows=[];summary={}
 for association in ASSOCIATIONS:
  losses={m:[] for m in MODELS};pairs={f'{a}__{b}':[] for a,b in itertools.combinations(MODELS,2)}
  scenario_minima={m:0 for m in MODELS};fixed_losses={m:[] for m in MODELS}
  for key,g in sorted(groups.items()):
   supports={(r['calibration_observed_mmwe'],r['heldout_observed_mmwe']) for r in g.values()}
   if len(supports)!=1:raise ValueError('unmatched field measurements')
   x,y=supports.pop();delta=90.
   methods={m:worst_error(x,y,g[m]['test_train_dose_ratio'],delta,association) for m in MODELS}
   for m in MODELS:
    losses[m].append(methods[m]['worst_absolute_error_mmwe'])
    fixed_losses[m].append(fixed_forecast(x,y,g[m]['test_train_dose_ratio'],delta,association)['worst_absolute_error_mmwe'])
   best=min(v['worst_absolute_error_mmwe'] for v in methods.values())
   for m in MODELS:scenario_minima[m]+=methods[m]['worst_absolute_error_mmwe']==best
   comparisons={}
   for a,b in itertools.combinations(MODELS,2):
    label=f'{a}__{b}';result=dominance(x,y,g[a]['test_train_dose_ratio'],g[b]['test_train_dose_ratio'],delta,association)
    comparisons[label]=result;pairs[label].append(result['status'])
   rows.append(dict(association=association,clock=[key[0],list(key[1])],method_losses=methods,pairwise=comparisons))
  worst={m:max(v) for m,v in losses.items()};best=min(worst.values())
  summary[association]=dict(fixed_already_issued_forecast_worst_errors_mmwe={m:max(v) for m,v in fixed_losses.items()},worst_over_measurements_and_all_clocks_mmwe=worst,
   retrospective_minimax_methods=[m for m,v in worst.items() if v==best],
   scenarios_with_smallest_conditional_worst_loss=scenario_minima,
   pairwise_status_counts={pair:{status:items.count(status) for status in sorted(set(items))} for pair,items in pairs.items()},
   strict_pairwise_dominance_all_clocks={pair:all(s=='a_strictly_better_everywhere' for s in ss) or all(s=='b_strictly_better_everywhere' for s in ss) for pair,ss in pairs.items()})
 result=dict(source_results_sha256=hashlib.sha256(raw).hexdigest(),stake='H1D',native_source_error_m=.1,native_source_density_kg_m3=900,conditional_delta_mmwe=90,
  uncertainty_semantics_author_confirmed=False,intervals_are_statistical_confidence_intervals=False,
  observations_or_models_retuned=False,retrospective_audit_not_deployment_selection=True,
  primary_estimand='input-error sensitivity: f*true_training_increment minus true_test_increment; not an already-issued fixed forecast',
  secondary_estimand='fixed forecast f*recorded_training_increment minus true_test_increment',rows=rows,summary=summary)
 out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');return summary
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('source',type=Path);p.add_argument('out',type=Path);a=p.parse_args();print(json.dumps(run(a.source,a.out),indent=2))
