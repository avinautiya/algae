"""Bind a source-role release to previously published immutable artifacts.

Publication is verified by the caller against GitHub; this script does not
authenticate an arbitrary supplied URL. Camera outcomes are never imported.
"""
import argparse,hashlib,json
from pathlib import Path
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('--spec',type=Path,required=True);a.add_argument('--predictions',type=Path,required=True)
    a.add_argument('--publication',required=True);a.add_argument('--role',choices=['train','evaluation'],required=True)
    a.add_argument('--out',type=Path,required=True);x=a.parse_args()
    spec=json.loads(x.spec.read_text());camera=spec['camera_evaluation']
    if x.out.exists():a.error('refuse overwriting release receipt')
    r=dict(schema='piramide-camera-unseal-v1',freeze_published=True,publication_reference=x.publication,
        frozen_spec_path=str(x.spec.resolve()),frozen_spec_sha256=sha(x.spec),predictions_path=str(x.predictions.resolve()),
        predictions_sha256=sha(x.predictions),access_role=x.role,
        prediction_stage='normalized_training' if x.role=='train' else 'final_evaluation',**camera)
    with x.out.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps(dict(role=x.role,predictions_sha256=r['predictions_sha256'],camera_outcomes_read=False)))
if __name__=='__main__':main()
