"""Locked Tower2 stake transfer. No evaluation outcomes before prediction seal."""
import argparse,json,hashlib,sys
from pathlib import Path
from anomaly_index import ROOT,sha,load_site,dose,model_grids
sys.path[:0]=[str(ROOT),str(ROOT/'moisture/review')]
from paths import require_piramide
OLD=require_piramide();sys.path.insert(0,str(OLD/'data/camera'))
from camera_loader import load_after_freeze,score_interval_mmwe

MODELS=['anomaly_exp_air','air_pdd','buried_pdd','gradient_index','constant_rate','shifted_ratio_air','centered_softplus_air','calendar_exp_air','moisture_only_exp']
def load(p):return json.loads(Path(p).read_text())
def write(p,d):
    with Path(p).open('x')as f:json.dump(d,f,indent=2,allow_nan=False);f.write('\n')

def prepare(directory):
    if directory.exists():raise ValueError('new freeze directory required')
    source=ROOT/'data/moisture/TOWER2_CAMERA_FREEZE_SPEC.json';c=load(source);fit=ROOT/'moisture/DEVELOPMENT_FIT.json'
    if sha(fit)!=c['T1_development_fit_sha256']:raise ValueError('development fit/source mismatch')
    params=c['frozen_model_parameters_from_T1development']
    learned={m:d['selected']['parameter']for m,d in load(fit)['fits'].items()}
    if set(params)!=set(MODELS)or any(params[m]!=p for m,p in learned.items()):raise ValueError('method or fitted hyperparameter changed')
    f=load_site('Tower2');expected=c['camera_evaluation']['intervals']
    if f['epochs']!=[expected[0]['start']]+[r['end']for r in expected]:raise ValueError('feature/camera endpoints differ')
    required=load(ROOT/'data/moisture/REQUIRED_FREEZE_FILES.json')['files']
    inputs={str((ROOT/'data/moisture'/n).relative_to(ROOT)):h for n,h in required.items()}
    inputs['moisture/DEVELOPMENT_FIT.json']=sha(fit)
    inputs['data/moisture/TOWER2_CAMERA_FREEZE_SPEC.json']=sha(source)
    src=['moisture/anomaly_index.py','moisture/experiment.py','moisture/review/moisture_scorer.py','paths.py','moisture/PROTOCOL.md']
    frozen=dict(primary_candidate='anomaly_exp_air',models=MODELS,parameters=params,development_fit_sha256=sha(fit),
        development_parameter_source='Tower1_released_periods_only',primary_metric='mean_absolute_period_ablation_error_mmwe',
        camera_evaluation=c['camera_evaluation'],source_hashes={p:sha(ROOT/p)for p in src},input_hashes=inputs,
        target_loader_sha256=sha(OLD/'data/camera/camera_loader.py'),
        external_source_hashes={p:sha(OLD/p)for p in ['data/decoded/debtemp_Tower2.csv','data/decoded/DECODE_VERIFICATION.json','data/camera/SEALED_SOURCE_MANIFEST.json','data/camera/Ablation_Tower2.mat']},
        prior_release_commit='c06d62cbe0eca34379d63d3c22d0c4f8a01c2f21',
        target_scope='calibrated stake transfer inside the same campaign; not independent glacier/weather replication',
        Tower2_camera_values_opened=False,novel_physics_claim=False,award_or_publication_guarantee=False)
    directory.mkdir(parents=True);write(directory/'FREEZE.json',frozen)

def verify(f):
    for p,h in {**f['source_hashes'],**f['input_hashes']}.items():
        if sha(ROOT/p)!=h:raise ValueError('frozen source/input changed: '+p)
    if sha(OLD/'data/camera/camera_loader.py')!=f['target_loader_sha256']:raise ValueError('target loader changed')
    for p,h in f['external_source_hashes'].items():
        if sha(OLD/p)!=h:raise ValueError('external source changed: '+p)

def normalized(directory,commit):
    if len(commit)!=40:raise ValueError('published initial freeze required')
    f=load(directory/'FREEZE.json');verify(f);features=load_site('Tower2')
    doses={m:dose(features,m,f['parameters'][m]).tolist()for m in MODELS}
    d=dict(preexecution_commit=commit,freeze_sha256=sha(directory/'FREEZE.json'),epochs=features['epochs'],doses=doses,
        parameters=f['parameters'],development_fit_sha256=f['development_fit_sha256'],source_max_air_gap_seconds=features['source_max_air_gap_seconds'],
        camera_target_values_opened=False,refitted_hyperparameters=False)
    write(directory/'NORMALIZED_PREDICTIONS.json',d)
    return dict(sha256=sha(directory/'NORMALIZED_PREDICTIONS.json'),models=len(doses),Tower2_camera_values_opened=False)

def calibrate(directory,receipt):
    from moisture_scorer import write_prediction_seal
    f=load(directory/'FREEZE.json');verify(f);d=load(directory/'NORMALIZED_PREDICTIONS.json');r=load(receipt)
    if r['access_role']!='train'or sha(directory/'NORMALIZED_PREDICTIONS.json')!=r['predictions_sha256']:raise ValueError('TRAIN release not bound to normalized predictions')
    if sha(directory/'FREEZE.json')!=r['frozen_spec_sha256']or d['freeze_sha256']!=r['frozen_spec_sha256']:raise ValueError('specification mismatch')
    frame=load_after_freeze(OLD/'data/camera',2,receipt);epochs=d['epochs'];observed=score_interval_mmwe(frame,epochs[0],epochs[1])
    y=observed['observed_ice_ablation_mmwe']
    if y<=0:raise ValueError('nonpositive TRAIN ablation; no invented positive gain')
    models={};gains={}
    for name,h in d['doses'].items():
        if h[0]<=0:raise ValueError('positive training dose required')
        gain=y/h[0];gains[name]=gain
        models[name]=[dict(start=a,end=b,predicted_mmwe=gain*q)for a,b,q in zip(epochs[1:-1],epochs[2:],h[1:])]
    op=dict(quantity='native_ice_ablation_mmwe',clock_basis='source_local_naive_no_utc_claim',observation_method=r['approved_operator'],
        ice_density_kg_m3=917.,camera_source_manifest_sha256=r['camera_manifest_sha256'],footprint='Tower2 marked ablation stake')
    p=dict(source_id='zenodo.11581073_Tower2_2022',primary_candidate=f['primary_candidate'],baseline_ids=[m for m in MODELS if m!=f['primary_candidate']],secondary_model_ids=[],
        models=models,target_operator=op,primary_metric=f['primary_metric'],tie_tolerance_mmwe=1e-6,
        heldout_targets_seen_before_prediction_freeze=False,development_fit_sha256=f['development_fit_sha256'],
        development_parameters={m:f['parameters'][m] for m in model_grids()},development_parameter_source=f['development_parameter_source'],
        local_train_gains=gains,local_train_target=observed,normalized_predictions_sha256=sha(directory/'NORMALIZED_PREDICTIONS.json'),
        no_Tower2_hyperparameter_training=True,prediction_is_empirical_index=True,physical_energy_conservation_claim=False,
        source_limits=['Contemporaneous observed moisture/temperature: conditional hindcast, not a future-weather forecast.',
            'Only constant additive moisture offsets cancel; gain/drift/nonlinear/media errors do not.',
            'Camera precision unknown; no independent weather/glacier replication from serial same-campaign periods.',
            'Different SM3 depths and within-accuracy moisture readings remain source limits.'])
    write(directory/'PREDICTIONS.json',p);write_prediction_seal(directory/'PREDICTIONS.json',directory/'PREDICTION_SEAL.json')
    return dict(prediction_sha256=sha(directory/'PREDICTIONS.json'),models=len(models),Tower2_evaluation_values_opened=False)

def score(directory,receipt,commit):
    from moisture_scorer import score_after_freeze
    p=load(directory/'PREDICTIONS.json');r=load(receipt)
    if r['access_role']!='evaluation'or sha(directory/'PREDICTIONS.json')!=r['predictions_sha256']:raise ValueError('final EVAL lock required')
    def truth():
        frame=load_after_freeze(OLD/'data/camera',2,receipt);rows=[]
        for interval in p['models'][p['primary_candidate']]:
            y=score_interval_mmwe(frame,interval['start'],interval['end'])
            rows.append(dict(start=interval['start'],end=interval['end'],observed_mmwe=y['observed_ice_ablation_mmwe']))
        return dict(source_id=p['source_id'],target_operator=p['target_operator'],intervals=rows,uncertainty_semantics='unknown')
    d=score_after_freeze(directory/'PREDICTIONS.json',directory/'PREDICTION_SEAL.json',expected_prediction_sha256=r['predictions_sha256'],freeze_commit=commit,truth_loader=truth)
    write(directory/'RESULTS.json',d);return d

if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__);a.add_argument('mode',choices=['freeze','normalized','calibrate','score']);a.add_argument('directory',type=Path);a.add_argument('--commit');a.add_argument('--receipt',type=Path);x=a.parse_args()
    if x.mode=='freeze':prepare(x.directory);print('freeze prepared; no Tower2 camera values opened')
    else:print(json.dumps(normalized(x.directory,x.commit)if x.mode=='normalized'else calibrate(x.directory,x.receipt)if x.mode=='calibrate'else score(x.directory,x.receipt,x.commit),indent=2))
