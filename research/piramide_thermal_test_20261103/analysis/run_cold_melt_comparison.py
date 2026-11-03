"""Two-stage comparison. Publish predictions before requesting evaluation truth."""
import argparse,csv,datetime as dt,hashlib,json,math,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'data/camera'),str(ROOT/'review')]
from calibration import fit_training_gain
from camera_loader import load_after_freeze,score_interval_mmwe
from scorer import write_prediction_seal,score_after_freeze

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):return json.loads(Path(p).read_text())
def dump(p,v):
    with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')
def positive_linear_integral(y0,y1,seconds):
    if y0>=0 and y1>=0:return seconds*(y0+y1)/2
    if y0<=0 and y1<=0:return 0.
    return seconds*max(y0,y1)**2/(2*abs(y1-y0))

def feature_dose(path,column,start,end,*,hold,max_gap_s):
    with Path(path).open() as stream:rows=list(csv.DictReader(stream))
    a,b=dt.datetime.fromisoformat(start),dt.datetime.fromisoformat(end)
    t=[dt.datetime.fromisoformat(r['SOURCE_LOCAL_LABEL']) for r in rows];v=[float(r[column]) for r in rows]
    if any(not math.isfinite(x) for x in v):raise ValueError('nonfinite feature data')
    if not a<b or a<t[0] or b>t[-1]:raise ValueError('feature support missing')
    total=0.;gaps=0.;maximum=0.
    for i in range(len(t)-1):
        x,y=max(a,t[i]),min(b,t[i+1])
        if y<=x:continue
        span=(t[i+1]-t[i]).total_seconds();seconds=(y-x).total_seconds()
        if span<=0 or span>max_gap_s:raise ValueError('unsupported feature gap; do not fill silently')
        maximum=max(maximum,span)
        if span>600:gaps+=seconds
        if hold=='linear':
            left=v[i]+(v[i+1]-v[i])*(x-t[i]).total_seconds()/span
            right=v[i]+(v[i+1]-v[i])*(y-t[i]).total_seconds()/span
            total+=positive_linear_integral(left,right,seconds)
        elif hold=='left':total+=max(v[i],0.)*seconds
        else:raise ValueError('unknown feature reconstruction')
    return dict(dose=total,positive_temperature_seconds=True,max_source_gap_s=maximum,
                duration_in_source_intervals_over_10min_s=gaps,
                reconstruction=hold,source_sha256=sha(path),column=column)

def fit(normalized,receipt,output,spec_path):
    if output.exists():raise ValueError('new output directory required')
    p=load(normalized);spec=load(spec_path);receipt_obj=load(receipt)
    if receipt_obj['access_role']!='train':raise ValueError('training receipt required')
    if sha(normalized)!=receipt_obj['predictions_sha256'] or Path(normalized).resolve()!=Path(receipt_obj['predictions_path']).resolve():
        raise ValueError('normalized predictions are not the artifact bound to the release receipt')
    frozen=load(receipt_obj['frozen_spec_path'])
    if sha(spec_path)!=frozen['analysis_spec_sha256']:raise ValueError('analysis specification changed after freeze')
    for name,h in {**frozen['analysis_source_hashes'],**frozen['analysis_input_hashes']}.items():
        if sha(ROOT/name)!=h:raise ValueError('analysis source/input changed after freeze: '+name)
    a,b,z=p['start_local'],p['split_local'],p['end_local']
    if [a,b]!=spec['primary_train'] or [b,z]!=spec['primary_evaluation']:raise ValueError('calendar split differs from frozen analysis')
    methods={'training_temperature_forward_fit':'forward','training_derivative_OLS_NOT_source_NTLS':'derivative_OLS'}
    expected={prefix+'__'+g['id']+'__'+h+'__'+cc['id'] for prefix in methods.values() for g in frozen['config']['geometries'] for h in frozen['config']['holds'] for cc in frozen['config']['capacities']}
    seen=set()
    for row in p['results']:
        if row['method'] not in methods:raise ValueError('unrecognized thermal method')
        ident=methods[row['method']]+'__'+row['geometry']['id']+'__'+row['hold']+'__'+row['capacity_id']
        if ident in seen:raise ValueError('duplicate thermal method/scenario')
        seen.add(ident)
    if seen!=expected:raise ValueError('thermal family incomplete or changed')
    # This is the first outcome access: evaluation endpoints are masked out by
    # the loader. No evaluation truth is passed into any selection or fitting.
    frame=load_after_freeze(ROOT/'data/camera',1,receipt)
    y=score_interval_mmwe(frame,a,b)
    train_seconds=(dt.datetime.fromisoformat(b)-dt.datetime.fromisoformat(a)).total_seconds()
    test_seconds=(dt.datetime.fromisoformat(z)-dt.datetime.fromisoformat(b)).total_seconds()
    models={};calibrations={};method_ids=[]
    for row in p['results']:
        if row['status']!='predicted' or not row['no_refreezing_melt_conversion_eligible']:
            raise ValueError('frozen thermal method refused/negative basal flux; no silent family exclusion')
        ident=methods[row['method']]+'__'+row['geometry']['id']+'__'+row['hold']+'__'+row['capacity_id']
        h0=row['training_normalized_basal_heat_K_m'];h1=row['evaluation_normalized_basal_heat_K_m']
        c=fit_training_gain([dict(start_s=0.,end_s=train_seconds,observed_mmwe=y['observed_ice_ablation_mmwe'],dose=h0)],
            train_end_s=train_seconds,evaluation_start_s=train_seconds,dose_quantity='normalized basal heat K m')
        c['scaled_capacity_diagnostic_NOT_physical_capacity_J_m3K']=c['gain_mmwe_per_dose']*spec['latent_heat_J_kg']
        c['postprocessing_gain_relative_to_physical_component']=c['gain_mmwe_per_dose']*spec['latent_heat_J_kg']/row['physicalcapacity_J_m3K']
        c['gain_is_not_identified_physical_heat_capacity']=True
        # The physical-range flag is never used to choose a more favorable score.
        models[ident]=[dict(start=b,end=z,predicted_mmwe=c['gain_mmwe_per_dose']*h1)]
        calibrations[ident]=c;method_ids.append(ident)
    model_id='constant_rate'
    models[model_id]=[dict(start=b,end=z,predicted_mmwe=y['observed_ice_ablation_mmwe']*test_seconds/train_seconds)]
    source_features={}
    for ident,path,column,hold,gap in [('positive_buried_T_DT1',ROOT/'data/decoded/debtemp_Tower1.csv','DT1','left',1800.),
            ('positive_air_TA3_linear_gapfill18h',ROOT/'data/decoded/trh_Tower1.csv','TA_3','linear',64800.)]:
        h0=feature_dose(path,column,a,b,hold=hold,max_gap_s=gap);h1=feature_dose(path,column,b,z,hold=hold,max_gap_s=gap)
        c=fit_training_gain([dict(start_s=0,end_s=train_seconds,observed_mmwe=y['observed_ice_ablation_mmwe'],dose=h0['dose'])],
            train_end_s=train_seconds,evaluation_start_s=train_seconds,dose_quantity='positive temperature K s')
        models[ident]=[dict(start=b,end=z,predicted_mmwe=c['gain_mmwe_per_dose']*h1['dose'])]
        calibrations[ident]=c;source_features[ident]=dict(train=h0,evaluation=h1)
    primary=spec['primary_candidate'];baselines=spec['baseline_ids']
    secondary=sorted(set(models)-{primary,*baselines})
    operator=dict(quantity='native_ice_ablation_mmwe',clock_basis='source_local_naive_no_utc_claim',
        observation_method=receipt_obj['approved_operator'],ice_density_kg_m3=917.,
        camera_source_manifest_sha256=receipt_obj['camera_manifest_sha256'],footprint='Tower1 marked ablation stake')
    prediction=dict(source_id='zenodo.11581073_Tower1_2022',models=models,primary_candidate=primary,
        baseline_ids=baselines,secondary_model_ids=secondary,target_operator=operator,
        primary_metric='absolute_total_ablation_error_mmwe',tie_tolerance_mmwe=spec['tie_tolerance_mmwe'],
        heldout_targets_seen_before_prediction_freeze=False,conditional_hindcast=True,
        calibration=y,calibrations=calibrations,feature_doses=source_features,
        normalized_predictions_sha256=sha(normalized),analysis_spec_sha256=sha(spec_path),
        source_and_model_limitations=spec['limitations'])
    output.mkdir(parents=True);dump(output/'PREDICTIONS.json',prediction)
    write_prediction_seal(output/'PREDICTIONS.json',output/'PREDICTION_SEAL.json')
    return dict(prediction_sha256=sha(output/'PREDICTIONS.json'),models=len(models),evaluation_targets_opened=False)

def score(directory,receipt,commit):
    path=directory/'PREDICTIONS.json';p=load(path);r=load(receipt)
    if r['access_role']!='evaluation' or sha(path)!=r['predictions_sha256']:raise ValueError('final evaluation lock required')
    def truth():
        f=load_after_freeze(ROOT/'data/camera',1,receipt)
        interval=p['models'][p['primary_candidate']][0]
        v=score_interval_mmwe(f,interval['start'],interval['end'])
        return dict(source_id=p['source_id'],target_operator=p['target_operator'],
            intervals=[dict(start=interval['start'],end=interval['end'],observed_mmwe=v['observed_ice_ablation_mmwe'])],
            uncertainty_semantics='unknown',source_diagnostics=v)
    result=score_after_freeze(path,directory/'PREDICTION_SEAL.json',expected_prediction_sha256=r['predictions_sha256'],freeze_commit=commit,truth_loader=truth)
    dump(directory/'RESULTS.json',result);return result

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('mode',choices=['fit','score'])
    parser.add_argument('--normalized',type=Path);parser.add_argument('--receipt',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--spec',type=Path);parser.add_argument('--commit')
    args=parser.parse_args()
    result=fit(args.normalized,args.receipt,args.out,args.spec) if args.mode=='fit' else score(args.out,args.receipt,args.commit)
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
