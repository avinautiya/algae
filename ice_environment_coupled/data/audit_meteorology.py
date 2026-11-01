"""Read observed AWS fields without inventing face radiation or stream forcing."""
import csv,datetime,hashlib,json,math
from collections import Counter
from pathlib import Path

FIELDS='year doy time timeDEC WS WD TA RH SWIN TS LWIN SWOUT LWOUT RAD_CS PRECIP'.split()
def load(path):
    records=[]
    with path.open(newline='') as stream:
        reader=csv.reader(stream,delimiter='\t')
        header=[x.strip() for x in next(reader)]
        if header!=FIELDS:raise ValueError('unrecognized AWS schema')
        last=None
        for row in reader:
            if not row:continue
            if len(row)!=len(FIELDS):raise ValueError('ragged AWS row')
            d={k:float(v.strip()) for k,v in zip(FIELDS,row)}
            year,doy,hhmm=(int(d[k]) for k in ('year','doy','time'))
            if any(d[k]!=int(d[k]) for k in ('year','doy','time')):raise ValueError('noninteger timestamp')
            hh,mm=divmod(hhmm,100)
            if not 0<=hh<24 or not 0<=mm<60:raise ValueError('invalid HHMM')
            date=datetime.datetime(year,1,1)+datetime.timedelta(days=doy-1,hours=hh,minutes=mm)
            if date.year!=year or (last is not None and date<=last):raise ValueError('invalid/nonmonotone date')
            if abs(d['timeDEC']-(doy+hh/24+mm/1440))>1e-5:raise ValueError('decimal day mismatch')
            missing=[k for k,v in d.items() if not math.isfinite(v)]
            flags=[]
            if math.isfinite(d['RH']) and not 0<=d['RH']<=100:flags.append('humidity_outside_percent_bounds')
            if math.isfinite(d['WS']) and d['WS']<0:flags.append('negative_wind_speed')
            if math.isfinite(d['PRECIP']) and d['PRECIP']<0:flags.append('negative_precipitation')
            records.append(dict(timestamp_local_label=date.isoformat(),values={k:v if math.isfinite(v) else None for k,v in d.items()},
                                missing_fields=missing,negative_shortwave=d['SWIN']<0,quality_flags=flags))
            last=date
    if not records:raise ValueError('empty AWS data')
    return records

def audit(root):
    p=root/'meteodata_AWS2013.txt';rows=load(p)
    summary={}
    for key in ('TA','RH','WS','SWIN','LWIN','PRECIP'):
        vals=[r['values'][key] for r in rows if r['values'][key] is not None]
        summary[key]=dict(valid=len(vals),missing=len(rows)-len(vals),min=min(vals) if vals else None,max=max(vals) if vals else None)
    result=dict(source='CliffEBM published Lirung AWS example',file_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
                rows=len(rows),first=rows[0]['timestamp_local_label'],last=rows[-1]['timestamp_local_label'],
                negative_shortwave_rows=sum(r['negative_shortwave'] for r in rows),fields=summary,
                humidity_out_of_bounds_rows=sum('humidity_outside_percent_bounds' in r['quality_flags'] for r in rows),
                incomplete_date_labels={day:n for day,n in Counter(r['timestamp_local_label'][:10] for r in rows).items() if n!=24},
                policy=['Timezone is unspecified here: labels are not converted to UTC.',
                        'Negative observed shortwave is retained and flagged, not silently clipped.',
                        'Missing values are preserved; no imputation or automatic field forecast.',
                        'Station radiation is not face radiation; slope/aspect/shading and radiation partition required.',
                        'Stream runoff, water temperature, pressure and observation geometry are not supplied by these fields.'])
    (root/'AWS_AUDIT.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result,indent=2))
    return result
if __name__=='__main__':audit(Path(__file__).parent)
