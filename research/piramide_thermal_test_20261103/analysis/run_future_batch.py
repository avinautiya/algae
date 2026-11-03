"""Fixed-gain future-period hindcast; never recalibrates after V3."""
import argparse,json,math,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'review'),str(ROOT/'data/camera')]
from run_cold_melt_comparison import sha,load,dump,feature_dose
from camera_loader import load_after_freeze,score_interval_mmwe
from batch_scorer import write_prediction_seal,score_after_freeze,pooled_with_prior_result

def build(normalized_path,freeze_path,output):
    if output.exists():raise ValueError('new prediction directory required')
    d=load(normalized_path);f=load(freeze_path);source_path=ROOT/'model/frozen_piramide_v3/melt_comparison/PREDICTIONS.json'
    if d['freeze_sha256']!=sha(freeze_path):raise ValueError('normalized batch does not bind this source freeze')
    p0=load(source_path)
    if sha(source_path)!=d['source_fixed_gain_predictions_sha256']:raise ValueError('fixed TRAIN-gain source changed')
    for name,h in {**f['analysis_source_hashes'],**f['analysis_input_hashes']}.items():
        if sha(ROOT/name)!=h:raise ValueError('analysis input/code changed after batch freeze: '+name)
    blocks=d['future_intervals'];models={};seen=set()
    if len(blocks)!=3:raise ValueError('exactly three source-eligible future blocks required')
    for row in d['results']:
        name=row['model_id']
        if name in seen or name not in p0['calibrations']:raise ValueError('changed/duplicate thermal family')
        seen.add(name);gain=p0['calibrations'][name]['gain_mmwe_per_dose']
        if [(x['start'],x['end'])for x in row['future_blocks']]!=[(x['start'],x['end'])for x in blocks]:raise ValueError('thermal blocks disagree')
        models[name]=[dict(start=x['start'],end=x['end'],predicted_mmwe=gain*x['normalized_basal_heat_K_m'])for x in row['future_blocks']]
    if len(seen)!=36:raise ValueError('entire frozen 36-case thermal family required')
    prior_config=load(ROOT/'model/frozen_piramide_v3/FREEZE.json')['config']
    from datetime import datetime
    duration=lambda a,b:(datetime.fromisoformat(b)-datetime.fromisoformat(a)).total_seconds()
    rate=p0['calibration']['observed_ice_ablation_mmwe']/duration(prior_config['start'],prior_config['split'])
    models['constant_rate']=[dict(start=x['start'],end=x['end'],predicted_mmwe=rate*duration(x['start'],x['end']))for x in blocks]
    sources={}
    for name,path,column,hold,gap in [('positive_buried_T_DT1',ROOT/'data/decoded/debtemp_Tower1.csv','DT1','left',1800.),
            ('positive_air_TA3_linear_gapfill18h',ROOT/'data/decoded/trh_Tower1.csv','TA_3','linear',64800.)]:
        gain=p0['calibrations'][name]['gain_mmwe_per_dose'];rows=[];feature=[]
        for b in blocks:
            x=feature_dose(path,column,b['start'],b['end'],hold=hold,max_gap_s=gap)
            rows.append(dict(start=b['start'],end=b['end'],predicted_mmwe=gain*x['dose']));feature.append(x)
        models[name]=rows;sources[name]=feature
    if set(models)!=set(p0['models']):raise ValueError('changed complete model family')
    p=dict(source_id=p0['source_id'],models=models,primary_candidate=p0['primary_candidate'],baseline_ids=p0['baseline_ids'],
        secondary_model_ids=p0['secondary_model_ids'],target_operator=p0['target_operator'],
        primary_metric='mean_absolute_period_ablation_error_mmwe',tie_tolerance_mmwe=p0['tie_tolerance_mmwe'],
        heldout_targets_seen_before_prediction_freeze=False,source_fixed_gain_predictions_sha256=sha(source_path),
        normalized_predictions_sha256=sha(normalized_path),fixed_gain_source_period=p0['calibration'],
        refitted_parameters=False,conditional_hindcast=True,feature_doses=sources,
        prior_V3_result_sha256=sha(ROOT/'model/frozen_piramide_v3/melt_comparison/RESULTS.json'),
        first_V3_air_baseline_loss_preserved=True,unsupported_fourth_period=d['excluded_intervals'])
    output.mkdir(parents=True);dump(output/'PREDICTIONS.json',p);write_prediction_seal(output/'PREDICTIONS.json',output/'PREDICTION_SEAL.json')
    return dict(prediction_sha256=sha(output/'PREDICTIONS.json'),models=len(models),future_periods=3,future_camera_outcomes_opened=False,parameters_refitted=False)

def score(directory,receipt_path,commit):
    path=directory/'PREDICTIONS.json';p=load(path);r=load(receipt_path)
    if r['access_role']!='evaluation' or sha(path)!=r['predictions_sha256']:raise ValueError('final future prediction lock required')
    def truth():
        frame=load_after_freeze(ROOT/'data/camera',1,receipt_path);intervals=[]
        for b in p['models'][p['primary_candidate']]:
            x=score_interval_mmwe(frame,b['start'],b['end'])
            intervals.append(dict(start=b['start'],end=b['end'],observed_mmwe=x['observed_ice_ablation_mmwe']))
        return dict(source_id=p['source_id'],target_operator=p['target_operator'],intervals=intervals,uncertainty_semantics='unknown')
    result=score_after_freeze(path,directory/'PREDICTION_SEAL.json',expected_prediction_sha256=r['predictions_sha256'],freeze_commit=commit,truth_loader=truth)
    dump(directory/'RESULTS.json',result)
    pooled=pooled_with_prior_result(result,ROOT/'model/frozen_piramide_v3/melt_comparison/RESULTS.json',expected_prior_result_sha256=p['prior_V3_result_sha256'])
    dump(directory/'SECONDARY_POOLED_WITH_V3.json',pooled)
    return result

def main():
    a=argparse.ArgumentParser(description=__doc__);a.add_argument('mode',choices=['build','score']);a.add_argument('--normalized',type=Path)
    a.add_argument('--freeze',type=Path);a.add_argument('--receipt',type=Path);a.add_argument('--out',type=Path,required=True);a.add_argument('--commit')
    x=a.parse_args();r=build(x.normalized,x.freeze,x.out)if x.mode=='build'else score(x.out,x.receipt,x.commit)
    print(json.dumps(r,indent=2))
if __name__=='__main__':main()
