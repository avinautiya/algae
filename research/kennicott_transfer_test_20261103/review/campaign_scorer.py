"""Independent strict same-operator campaign scorer; receives released typed rows.
No archive loader, duration guessing, density conversion, fit or target selection.
"""
import hashlib,json,math

def canonical_sha(obj):return hashlib.sha256(json.dumps(obj,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def score(predictions,targets,receipt):
 if receipt.get('access_role')!='evaluation' or receipt.get('source_binding_resolved')is not True:raise ValueError('Source binding and evaluation release required')
 if receipt.get('frozen_predictions_sha256')!=canonical_sha(predictions):raise ValueError('Prediction freeze mismatch')
 if receipt.get('released_typed_targets_sha256')!=canonical_sha(targets):raise ValueError('Released target hash mismatch')
 operator=predictions.get('target_operator');other=targets.get('target_operator')
 if operator!=other:raise ValueError('Target operator mismatch')
 if not isinstance(operator,dict)or operator.get('kind')not in ('campaign_total','campaign_mean_daily_rate'):raise ValueError('Explicit native operator required')
 if operator.get('quantity') not in ('ice_ablation','water_equivalent_ablation','net_ice_ablation_loss_conditional') or not operator.get('unit')or not operator.get('source_operator_reference'):raise ValueError('Quantity, units, source reference required')
 if operator['kind']=='campaign_mean_daily_rate' and operator.get('rate_definition_verified')is not True:raise ValueError('Mean daily-rate definition unresolved')
 tr=targets.get('rows',[]);pr=predictions.get('models',{})
 if not tr or not pr:raise ValueError('No declared rows/models')
 keys=[x['case_id']for x in tr]
 if len(keys)!=len(set(keys)):raise ValueError('Duplicate target cases')
 for row in tr:
  if not all(row.get(k)for k in ('campaign_id','source_stake_id','source_thermal_binding_id','period_start','period_end')):raise ValueError('Unbound native footprint/clock')
  if not math.isfinite(row['observed']):raise ValueError('Missing target')
 results={}
 for name,rows in pr.items():
  if len(rows)!=len(tr)or len(set(r['case_id']for r in rows))!=len(rows)or set(r['case_id']for r in rows)!=set(keys):raise ValueError('Prediction-case mismatch')
  by={r['case_id']:r for r in rows};errs=[]
  for r in tr:
   p=by[r['case_id']]
   if p.get('source_thermal_binding_id')!=r['source_thermal_binding_id']:raise ValueError('Thermal footprint binding mismatch')
   if not math.isfinite(p['predicted']):raise ValueError('Missing prediction; retain abstention separately')
   e=p['predicted']-r['observed'];errs.append({'case_id':r['case_id'],'campaign_id':r['campaign_id'],'predicted':p['predicted'],'observed':r['observed'],'signed_error':e,'absolute_error':abs(e)})
  campaigns=sorted(set(r['campaign_id']for r in tr));cm={c:sum(x['absolute_error']for x in errs if x['campaign_id']==c)/sum(x['campaign_id']==c for x in errs)for c in campaigns}
  results[name]={'MAE_native_unit':sum(x['absolute_error']for x in errs)/len(errs),'mean_signed_error_native_unit':sum(x['signed_error']for x in errs)/len(errs),'campaign_balanced_MAE_native_unit':sum(cm.values())/len(cm),'per_campaign_MAE_native_unit':cm,'all_rows':errs}
 return {'target_operator':operator,'models':results,'n_native_rows':len(tr),'n_campaign_groups':len(set(r['campaign_id']for r in tr)),'independent_glacier_count_not_inferred':True,'statistical_CI':None,'scope':'Same-operator source-bound evaluation only; no independent-validation/novelty/superiority inference is automatic'}
