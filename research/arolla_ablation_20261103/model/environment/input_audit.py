"""Target-free audit: no model fit, prediction, or held-out target scoring."""
import argparse, hashlib, json
from datetime import datetime
from pathlib import Path

def audit(path):
    raw=Path(path).read_bytes(); data=json.loads(raw)
    # All possible visit hours for earliest primary visit through latest final visit.
    lo=datetime.fromisoformat('2010-07-28T00:00:00+02:00')
    hi=datetime.fromisoformat('2010-09-12T00:00:00+02:00')
    rows=[r for r in data['forcing'] if lo <= datetime.fromisoformat(r['source_label_local']) < hi]
    keys=sorted(set().union(*(r['features'] for r in rows)))
    counts={k:{'present':sum(r['features'].get(k) is not None for r in rows),
               'missing':sum(r['features'].get(k) is None for r in rows)} for k in keys}
    return {'source_sha256':hashlib.sha256(raw).hexdigest(), 'scope':'target-free source availability only; no prediction or score',
            'source_label_window':'2010-07-28 00:00 <= label < 2010-09-12 00:00 UTC+02',
            'hourly_rows':len(rows),'counts':counts,
            'used_by_baselines':{'PDD':['T'],'constant_rate':['elapsed time'],'net_shortwave':['SWIN','SWOUT']},
            'phase_gate':'Native post-July4 snow depth is unavailable. None -> 0 forbidden.',
            'property_provenance':data['debris']}

if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('source'); p.add_argument('--output'); a=p.parse_args()
    s=json.dumps(audit(a.source),indent=2)+'\n'
    if a.output: Path(a.output).write_text(s)
    else: print(s,end='')
