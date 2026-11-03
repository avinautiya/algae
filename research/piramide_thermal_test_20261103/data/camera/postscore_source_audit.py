"""Additive post-release source audit; never run before a valid role receipt.

Recomputes raw endpoint conversion independently and compares the supplied
ablation column's endpoint scaling. It NEVER changes the frozen operator.
Only the authorized role's endpoints are accessed, not the full outcome series.
"""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd
from camera_loader import verify_receipt,authorized_labels,sha,OPERATOR

def audit(camera_dir,tower,receipt):
    folder=Path(camera_dir);r=verify_receipt(receipt,folder)
    labels=authorized_labels(r,tower)
    m=json.loads((folder/'SEALED_SOURCE_MANIFEST.json').read_text())
    filename=f'Ablation_Tower{tower}.mat';entry=next(x for x in m['members'] if x['file']==filename)
    if sha(folder/filename)!=entry['sha256']:raise ValueError('Source bytes differ from sealed member')
    from matio import load_from_mat
    obj=load_from_mat(folder/filename,raw_data=True)['tracks'];a=obj.properties['any']
    names=[str(x.item()) for x in a['varNames'][0,0].ravel()]
    if names!=['distance','ablation']:raise ValueError('Unexpected schema')
    rt=a['rowTimes'][0,0];epoch=np.asarray(rt.properties['data']).ravel()
    if not np.isfinite(epoch).all() or np.iscomplexobj(epoch) or not (epoch==np.rint(epoch)).all():raise ValueError('Unsupported epoch')
    index=pd.DatetimeIndex(epoch.astype('int64').astype('datetime64[ms]'))
    positions=index.get_indexer(labels)
    if (positions<0).any():raise ValueError('Missing frozen endpoint')
    # First native outcome read; all roles/spec/source checks have completed.
    cells=a['data'][0,0].ravel()
    values={name:np.asarray(cells[names.index(name)]).ravel()[positions] for name in names}
    if not all(np.isfinite(v).all() for v in values.values()):raise ValueError('Nonfinite released endpoint')
    endpoints=[{'source_local_label':str(t),'distance_source_m':float(values['distance'][i]),
                'provided_ablation_source_units_unverified':float(values['ablation'][i])}
               for i,t in enumerate(labels)]
    rows=[]
    for interval in r['intervals']:
        if interval['tower']!=tower or interval['role']!=r['access_role']:continue
        i=labels.index(pd.Timestamp(interval['start']));j=labels.index(pd.Timestamp(interval['end']))
        delta_distance=float(values['distance'][j]-values['distance'][i])
        delta_supplied=float(values['ablation'][j]-values['ablation'][i])
        ratio=delta_supplied/delta_distance if delta_distance!=0 else None
        rows.append({'start':interval['start'],'end':interval['end'],'raw_distance_change_m_ice':delta_distance,
                    'independent_ice_conversion_mmwe':delta_distance*917.,
                    'provided_ablation_change_source_units_unverified':delta_supplied,
                    'provided_delta_over_distance_delta':ratio,
                    'negative_distance_change':delta_distance<0,
                    'frozen_target_operator':OPERATOR,'operator_changed':False,
                    'candidate_expected_scalings_not_selected':{'m_water_per_m_ice':.917,'mm_water_per_m_ice':917.},
                    'reset_assessment':'Only released endpoints inspected; interior reset/manual correction cannot be ruled out'})
    return {'scope':f'Authorized {r["access_role"]} endpoints only','receipt_sha256':sha(receipt),
            'source_sha256':entry['sha256'],'source_fields_have_units':False,
            'endpoints':endpoints,'intervals':rows,
            'uncertainty':'No quantitative camera precision established. Scaling agreement cannot establish camera accuracy.',
            'validation_claim':'Source/operator audit only; no new prediction comparison or model tuning'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('camera_dir');p.add_argument('receipt');p.add_argument('--tower',type=int,default=1);p.add_argument('--output',required=True);args=p.parse_args()
    report=audit(args.camera_dir,args.tower,args.receipt)
    Path(args.output).write_text(json.dumps(report,indent=2)+'\n')
    print('Role-authorized source audit written; no target values printed.')
