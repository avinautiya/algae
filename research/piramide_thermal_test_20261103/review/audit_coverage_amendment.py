"""After-release integrity audit of V3's availability-only endpoint amendment.
No camera outcomes are read; prediction artifacts contain known TRAIN calibration.
"""
import argparse,datetime as dt,hashlib,json
from pathlib import Path

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def ident(row):return (row['method'],row['geometry']['id'],row['hold'],row['capacity_id'])
def nearest_finite(metadata,original):
 center=dt.datetime.fromisoformat(original);eligible=[]
 for row in metadata:
  t=dt.datetime.fromisoformat(row['source_local_label'])
  if row['distance_isfinite'] is True and abs((t-center).total_seconds())<=43200:eligible.append(t)
 if not eligible:raise ValueError('no finite eligible endpoint; remain unscorable')
 return min(eligible,key=lambda t:(abs((t-center).total_seconds()),t)).isoformat()

def audit(current_freeze,prior_freeze,current_normalized,prior_normalized,current_predictions,prior_predictions,v2_failure,missingness_audit,*,evaluation_released=False):
 if evaluation_released is not True:raise ValueError('parent must announce evaluation release before artifact audit')
 f0,f1,n0,n1,p0,p1,failed,missing=[load(q) for q in [prior_freeze,current_freeze,prior_normalized,current_normalized,prior_predictions,current_predictions,v2_failure,missingness_audit]]
 checks=[]
 def check(ok,description):checks.append(dict(pass_=bool(ok),check=description))
 amendment=f1['coverage_amendment'];old=f0['config'];new=f1['config']
 check(amendment['original_freeze_sha256']==sha(prior_freeze),'original physical freeze preserved and linked')
 check(amendment['missingness_audit_sha256']==sha(missingness_audit),'availability audit linked by hash')
 check(failed['status']=='unscorable' and failed['evaluation_error_computed'] is False and failed['score_written'] is False,'V2 missing endpoint has no scored error')
 check(not (Path(prior_predictions).parent/'RESULTS.json').exists(),'V2 has no silently substituted result artifact')
 check(missing['missing_value_native_nan_confirmed'] is True and missing['new_target_values_exposed'] is False and missing['prediction_error_computed'] is False,'frozen availability audit used missingness only')
 selected=nearest_finite(missing['nearby_timestamps_finiteness_only'],missing['original_end'])
 check(selected==missing['selected_end']==amendment['selected_end']==new['end'],'independent nearest-finite/tie-earlier endpoint selection')
 check(old['end']==missing['original_end']==amendment['original_end'],'original endpoint retained in deviation record')
 check({k:v for k,v in old.items() if k not in ['end','camera_evaluation']}=={k:v for k,v in new.items() if k not in ['end','camera_evaluation']},'only physical-horizon and camera-coverage config changed')
 camera_expected=json.loads(json.dumps(old['camera_evaluation']))
 for interval in camera_expected['intervals']:
  if interval['role']=='evaluation':
   check(interval['end']==old['end'],'old authorized EVAL endpoint matches original horizon')
   interval['end']=new['end']
 check(camera_expected==new['camera_evaluation'],'only camera EVAL end permission changed; TRAIN roles/operator unchanged')

 for category in ['source_hashes','input_hashes','analysis_source_hashes','analysis_input_hashes']:
  a,b=f0.get(category,{}),f1.get(category,{})
  check(all(b.get(k)==v for k,v in a.items()),'all previously frozen '+category+' unchanged')
 check(p0['primary_candidate']==p1['primary_candidate'] and p0['baseline_ids']==p1['baseline_ids'] and p0['secondary_model_ids']==p1['secondary_model_ids'],'same primary and full model family')
 check(len(p0['models'])==len(p1['models'])==39 and set(p0['models'])==set(p1['models']),'all39 issued methods retained')
 check(p0['calibration']==p1['calibration'] and p0['calibrations']==p1['calibrations'],'all training targets/doses/gains/diagnostics unchanged')
 check(p0['target_operator']==p1['target_operator'],'native observation operator/density/footprint unchanged')
 rows0={ident(r):r for r in n0['results']};rows1={ident(r):r for r in n1['results']}
 check(len(rows0)==len(rows1)==36 and set(rows0)==set(rows1),'all36 physical scenario IDs retained')
 immutable=['method','geometry','hold','capacity_id','physicalcapacity_J_m3K','capacity_J_m3K','alpha_m2_s','debris_conductivity_W_mK','alpha_refitted','training_normalized_basal_heat_K_m','raw_training_melt_mmwe']
 for key,newrow in rows1.items():
  oldrow=rows0[key]
  check(all(newrow[k]==oldrow[k] for k in immutable),'identical parameters/TRAIN heat '+str(key))
  a,b=oldrow['normalized_heat_trace'],newrow['normalized_heat_trace']
  check(len(b)<=len(a) and b==a[:len(b)],'causal physical trace exact unchanged prefix '+str(key))
 for key in ['positive_buried_T_DT1','positive_air_TA3_linear_gapfill18h']:
  check(p0['feature_doses'][key]['train']==p1['feature_doses'][key]['train'],'PDD training forcing/dose unchanged '+key)
 return dict(all_checks_pass=all(c['pass_'] for c in checks),checks=checks,original_end=old['end'],amended_end=new['end'],endpoint_shift_s=(dt.datetime.fromisoformat(new['end'])-dt.datetime.fromisoformat(old['end'])).total_seconds(),v2_formally_unscorable=True,v3_is_transparent_protocol_amendment_not_original_freeze_success=True,all_parameters_and_train_gains_unchanged=all(c['pass_'] for c in checks if any(w in c['check'] for w in ['parameters/TRAIN','training targets/doses/gains'])),raw_camera_values_read=False,performance_based_endpoint_selection_established=False,statistical_confidence_interval=False)
if __name__=='__main__':
 parser=argparse.ArgumentParser(description=__doc__)
 for key in ['current-freeze','prior-freeze','current-normalized','prior-normalized','current-predictions','prior-predictions','v2-failure','missingness-audit','out']:parser.add_argument('--'+key,type=Path,required=True)
 parser.add_argument('--evaluation-released',action='store_true');a=parser.parse_args()
 if a.out.exists():raise ValueError('refuse overwriting audit')
 result=audit(a.current_freeze,a.prior_freeze,a.current_normalized,a.prior_normalized,a.current_predictions,a.prior_predictions,a.v2_failure,a.missingness_audit,evaluation_released=a.evaluation_released)
 a.out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps({k:result[k] for k in ['all_checks_pass','endpoint_shift_s','all_parameters_and_train_gains_unchanged']},indent=2))
