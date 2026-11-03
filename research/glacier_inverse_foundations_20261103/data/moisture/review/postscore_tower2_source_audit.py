"""Run ONLY after parent announces the requested Tower2 role release.

The old sealed-camera loader validates published hashes/role before native
outcome access. This additive wrapper does not change frozen model code.
"""
import argparse,hashlib,json,sys,math
from pathlib import Path
ROOT=Path(__file__).parents[3];sys.path.insert(0,str(ROOT))
from paths import require_piramide
ANCHORS=['2022-12-02T03:00:00','2022-12-16T03:00:00','2022-12-30T03:00:00','2023-01-13T03:00:00','2023-01-27T03:00:00','2023-02-10T03:00:00']

def run(role,receipt_path,output):
    if role not in ['train','evaluation']:raise ValueError('Explicitrole required')
    receipt=json.loads(Path(receipt_path).read_text())
    if receipt.get('access_role')!=role:raise ValueError('Requestedrole differsfrompublishedreceipt')
    expected=list(zip(ANCHORS[:-1],ANCHORS[1:]));expected=expected[:1] if role=='train' else expected[1:]
    actual=[(x['start'],x['end']) for x in receipt['intervals'] if x['role']==role and x['tower']==2]
    if actual!=expected:raise ValueError('FrozenTower2epochs changed; no substitution permitted')
    old=require_piramide();camera=old/'data/camera';sys.path.insert(0,str(camera))
    from postscore_source_audit import audit
    raw=audit(camera,2,receipt_path)
    if [(x['start'],x['end']) for x in raw['intervals']]!=expected:raise ValueError('Sourceaudit scope mismatch')
    unit_checks=[]
    for x in raw['intervals']:
        distance=x['raw_distance_change_m_ice'];mwe=x['provided_ablation_change_source_units_unverified']
        raw_mmwe=distance*917.;provided_mmwe=mwe*1000.
        if not math.isfinite(raw_mmwe) or not math.isfinite(provided_mmwe):raise ValueError('Nonfinitefrozenendpoint; originalrefusalrule retained')
        unit_checks.append({'start':x['start'],'end':x['end'],'distance_delta_times917_mmwe':raw_mmwe,'provided_ablation_delta_times1000_mmwe':provided_mmwe,
          'absolute_unit_difference_mmwe':abs(raw_mmwe-provided_mmwe),'unit_scaling_consistent':math.isclose(raw_mmwe,provided_mmwe,rel_tol=1e-11,abs_tol=1e-9),
          'operator_changed':False,'negative_endpoint_difference':x['negative_distance_change']})
    result={'scope':'OnlypublishedroleTower2endpointsourceaudit; no model fit or alternate scoringoperator',
      'role':role,'receipt_sha256':hashlib.sha256(Path(receipt_path).read_bytes()).hexdigest(),'source_audit':raw,'unit_checks':unit_checks,
      'excluded_Tower1Feb24_and_all_Tower3_outcomes_accessed':False,'camera_endpoint_precision':'Unknown; no confidenceinterval',
      'corrections_reset_scope':'Onlyauthorizedendpoints; fullinteriorcorrectionandresetcontinuity unverified',
      'moisture_QA':'ProcesssoiltemperatureQCflagsareproxychecks, notabsolutefieldaccuracycertification; independentdriftuncertaintyremains'}
    with Path(output).open('x')as f:json.dump(result,f,indent=2);f.write('\n')
    print('Role-authorizedTower2sourceaudit saved; no target values printed.')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('role',choices=['train','evaluation']);p.add_argument('--receipt',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    name='TRAIN_RELEASE_RECEIPT.json'if a.role=='train'else'EVALUATION_RELEASE_RECEIPT.json'
    receipt=a.receipt or ROOT/'moisture/frozen_tower2_final'/name
    run(a.role,receipt,a.output)
