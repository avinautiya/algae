"""Post-freeze native camera loader. No CLI prints outcome values.

Only raw exposed-stake distance is an authorized scoring operator. Source
README distance[m] + final-paper exposed-stake operator establish m ICE;
provided ablation column units are conflicting/unverified and never used.
"""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd

OPERATOR='exposed_stake_distance_endpoint_difference_m_ice_times_917_to_mmwe'

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def verify_receipt(receipt_path,camera_dir):
    r=json.loads(Path(receipt_path).read_text())
    if r.get('schema')!='piramide-camera-unseal-v1' or r.get('freeze_published') is not True:raise ValueError('Published prediction freeze is required before reading outcomes')
    if not isinstance(r.get('publication_reference'),str) or not r['publication_reference'].strip():raise ValueError('Publication reference required')
    for pathkey,hashkey in [('frozen_spec_path','frozen_spec_sha256'),('predictions_path','predictions_sha256')]:
        if sha(r[pathkey])!=r[hashkey]:raise ValueError('Frozen artifact hash mismatch')
    manifest=Path(camera_dir)/'SEALED_SOURCE_MANIFEST.json'
    if sha(manifest)!=r.get('camera_manifest_sha256'):raise ValueError('Camera source manifest mismatch')
    if r.get('approved_operator')!=OPERATOR:raise ValueError('Unverified target operator or units')
    if r.get('clock_interpretation')!='same_source_naive_local_no_UTC_assumption':raise ValueError('Explicit native-clock interpretation required')
    spec=json.loads(Path(r['frozen_spec_path']).read_text())
    required={'approved_operator':r['approved_operator'],'clock_interpretation':r['clock_interpretation'],'camera_manifest_sha256':r['camera_manifest_sha256'],'intervals':r.get('intervals')}
    if spec.get('camera_evaluation')!=required:raise ValueError('Camera endpoints/operator are not bound into frozen spec')
    role=r.get('access_role')
    if role not in ['train','evaluation']:raise ValueError('Explicit train/evaluation access role required')
    stage={'train':'normalized_training','evaluation':'final_evaluation'}[role]
    if r.get('prediction_stage')!=stage:raise ValueError('Prediction publication stage does not authorize this role')
    endpoints=r.get('intervals')
    if not isinstance(endpoints,list) or not endpoints:raise ValueError('Predeclared evaluation endpoints required')
    for x in endpoints:
        if set(x)!=set(['tower','start','end','role']) or x['role'] not in ['train','evaluation'] or x['tower'] not in [1,2,3] or pd.Timestamp(x['end'])<=pd.Timestamp(x['start']):raise ValueError('Invalid frozen interval')
    return r

def authorized_labels(receipt,tower):
    intervals=[x for x in receipt['intervals'] if x['tower']==tower and x['role']==receipt['access_role']]
    if not intervals:raise ValueError('Tower not authorized for this access role')
    return sorted(set(pd.Timestamp(x[k]) for x in intervals for k in ['start','end']))

def load_after_freeze(camera_dir,tower,receipt_path):
    d=Path(camera_dir);r=verify_receipt(receipt_path,d)
    labels=authorized_labels(r,tower)
    file=f'Ablation_Tower{tower}.mat';m=json.loads((d/'SEALED_SOURCE_MANIFEST.json').read_text())
    member=next(x for x in m['members'] if x['file']==file)
    if sha(d/file)!=member['sha256']:raise ValueError('Camera member hash mismatch')
    from matio import load_from_mat
    raw=load_from_mat(d/file,raw_data=True)['tracks']
    a=raw.properties['any']
    names=[str(x.item()) for x in a['varNames'][0,0].ravel()]
    if names!=['distance','ablation']:raise ValueError('Unexpected native camera schema')
    rt=a['rowTimes'][0,0]; tz=rt.properties.get('tz',rt.properties.get('tmz'))
    if np.asarray(tz).size:raise ValueError('Clock changed from audited empty timezone')
    epoch=np.asarray(rt.properties['data']).ravel()
    if not np.isfinite(epoch).all() or np.iscomplexobj(epoch) or not (epoch==np.rint(epoch)).all():raise ValueError('Unsupported native clock')
    times=epoch.astype('int64').astype('datetime64[ms]')
    # First dereference of native outcome data occurs HERE, after receipt checks.
    columns=a['data'][0,0].ravel()
    native_index=pd.DatetimeIndex(times,name='SOURCE_LOCAL_LABEL')
    positions=native_index.get_indexer(labels)
    if (positions<0).any():raise ValueError('Predeclared endpoint absent; no interpolation')
    distance_column=np.asarray(columns[names.index('distance')]).ravel()
    if len(distance_column)!=len(times):raise ValueError('Outcome/time shape mismatch')
    distance=distance_column[positions]
    frame=pd.DataFrame({'distance_m_ice':distance},index=native_index[positions])
    frame.attrs['access_role']=r['access_role']
    frame.attrs['release_scope']='Only frozen endpoints for authorized role; no full target series returned'
    if not frame.index.is_unique or not frame.index.is_monotonic_increasing:raise ValueError('Invalid native time axis')
    return frame

def score_interval_mmwe(frame,start,end):
    """Exact native endpoint difference. Missing endpoints never interpolated.

    Density917kg/m³ times metres ICE gives kg/m² == mm WATER equivalent.
    Negative differences are returned and flagged, never clipped or reset.
    """
    s,e=pd.Timestamp(start),pd.Timestamp(end)
    if s.tz is not None or e.tz is not None or e<=s:raise ValueError('Require increasing naive source-clock endpoints')
    if s not in frame.index or e not in frame.index:raise ValueError('Frozen endpoint absent from camera record')
    x0,x1=float(frame.loc[s,'distance_m_ice']),float(frame.loc[e,'distance_m_ice'])
    if not np.isfinite([x0,x1]).all():raise ValueError('Native endpoint missing/nonfinite')
    subset=frame.loc[s:e,'distance_m_ice']
    steps=np.diff(subset.to_numpy(dtype=float))
    return {'observed_ice_ablation_mmwe':(x1-x0)*917.,'operator':OPERATOR,
            'negative_endpoint_difference':bool(x1<x0),
            'negative_released_increments':int(np.sum(steps<0)),
            'missing_released_points':int(subset.isna().sum()),
            'diagnostic_scope':'Released endpoints only, not complete native chronology',
            'uncertainty':'Camera precision/reset semantics not independently established; no confidence interval'}
