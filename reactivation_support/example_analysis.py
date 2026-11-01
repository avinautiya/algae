"""Small pinned-snapshot run solely to exercise the completed-output tools."""
from pathlib import Path
import json,platform
from .snapshot import cases
from .campaign_analysis import summarize
from .plot_campaign import plot

if __name__=='__main__':
    p=Path(__file__).parent/'records';p.mkdir(exist_ok=True)
    r=cases.run('strip',dx=.05,dtm=6*3600,days=.5,removal_day=.3)
    if platform.system()=='Darwin':
        r['snapshot_maxrss_raw_labelled_MB']=r['maxrss_MB']
        r['maxrss_MB']/=1024
        r['maxrss_note']='Normalized Darwin ru_maxrss bytes; pinned source assumes Linux KiB.'
    r['evidence_class']='pinned c70c553 development fixture; failed verification gate; NOT final evaluation'
    f=p/'example_development_run.json';f.write_text(json.dumps(r,indent=2,allow_nan=False))
    a=summarize(r,threshold=.001,persistence_hours=24)
    (p/'example_analysis.json').write_text(json.dumps(a,indent=2))
    plot([f],Path(__file__).parent/'figures/example_development',label='Pinned-snapshot development fixture — NOT final evaluation')
    print('requested removal day',r['removal_day'],'actual export day',r['export'][0]['t']/86400)
    print('analysis ordinary-flow start',a['ordinary_from_day'])
