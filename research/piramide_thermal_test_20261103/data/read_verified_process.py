"""Read decoded process CSV only after manifest/hash/cadence validation."""
import hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd

def read_process(directory,family='debtemp',tower='Tower1'):
    d=Path(directory);m=json.loads((d/'DECODE_VERIFICATION.json').read_text())
    row=next(x for x in m['checks'] if x['family']==family and x['tower']==tower)
    p=d/row['file']
    if hashlib.sha256(p.read_bytes()).hexdigest()!=row['sha256']:raise ValueError('Process CSV hash mismatch')
    frame=pd.read_csv(p,index_col='SOURCE_LOCAL_LABEL',parse_dates=True)
    if len(frame)!=row['rows'] or list(frame.columns)!=row['columns']:raise ValueError('Shape/column mismatch')
    if frame.index.tz is not None:raise ValueError('No authoritative native timezone is present; UTC relabeling forbidden')
    if not frame.index.is_unique or not frame.index.is_monotonic_increasing:raise ValueError('Invalid time axis')
    if family in ['debtemp','sm'] and not np.all(np.diff(frame.index.values).astype('timedelta64[s]').astype(int)==1800):raise ValueError('Unexpected native cadence')
    return frame
