"""Replay ONLY the frozen training fit to serialize ridge weights; no selection changes."""
import json
import numpy as np
from pathlib import Path
from run_benchmark import load_profile,masks,ridge_features,ridge_fit,ridge_predict,sha,put
ROOT=Path(__file__).resolve().parent
r=json.loads((ROOT/'results/FIELD_TEMPERATURE_RESULTS.json').read_text());out={};maxdiff=0.
for p in r['profile_results']:
 t,dates,z,v=load_profile(p['split']);train,test=masks(t);idx=np.flatnonzero(train|test);X=ridge_features(v,idx);Y=v[idx,1:-1]
 fit=ridge_fit(X[train[idx]],Y[train[idx]],p['ridge_penalty'])
 prediction=ridge_predict(fit,X)
 import csv
 with (ROOT/'results'/p['prediction_csv']).open() as f:rows=list(csv.DictReader(f))
 saved=np.array([float(row['ridge_C']) for row in rows]).reshape(len(idx),len(z)-2)
 delta=float(np.max(np.abs(saved-prediction)));maxdiff=max(maxdiff,delta)
 if delta>1e-12:raise ValueError('training-only replay differs from original ridge predictions')
 out[p['profile_id']]={k:a.tolist() for k,a in fit.items()}
put(ROOT/'results/RIDGE_PARAMETERS.json',dict(scope='deterministic serialization of original selected training-only fit; no model selection or prediction change',executed_runner_sha256=sha(ROOT/'run_benchmark.py'),max_replay_prediction_difference_C=maxdiff,models=out,feature_order=[f'{side}_lag{lag}h' for lag in [0,1,3,6,12,24] for side in ['upper','lower']]))
print('Ridge parameters serialized; max difference',maxdiff)
