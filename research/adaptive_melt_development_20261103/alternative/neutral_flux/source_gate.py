"""Target-free native support audit; no reconstruction and no ablation read."""
import csv,datetime as dt,json,math,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
def stamp(x):return dt.datetime.fromisoformat(x)
def valid(row,tier):
    try: finite=all(math.isfinite(float(row[k])) for k in ('H_up_W_m2','LE_up_W_m2','pressure_Pa','gradient_Ri_at1m_approx'))
    except (ValueError,TypeError): return False
    if not finite or row['status']!='conditional_neutral_proxy_only':return False
    if tier in ('near_neutral','measurement_resolved') and abs(float(row['gradient_Ri_at1m_approx']))>.03:return False
    if tier=='measurement_resolved' and not all(row[k]=='True' for k in ('wind_slope_sign_resolved_given_stated_accuracy','T_slope_sign_resolved_given_stated_accuracy')):return False
    return True

def support(rows,start,end,tier,max_gap=1800):
    if end<=start:raise ValueError('Nonpositive interval')
    total=(end-start).total_seconds();supported=0.;segments=0
    for a,b in zip(rows[:-1],rows[1:]):
        ta,tb=stamp(a['SOURCE_LOCAL_LABEL']),stamp(b['SOURCE_LOCAL_LABEL'])
        if tb<=ta:raise ValueError('Non-increasing source labels')
        overlap=max(0.,(min(tb,end)-max(ta,start)).total_seconds())
        if overlap and (tb-ta).total_seconds()<=max_gap and valid(a,tier) and valid(b,tier):supported+=overlap;segments+=1
    return {'period_seconds':total,'supported_seconds':supported,'unsupported_seconds':total-supported,'supported_time_fraction':supported/total,'supported_segments':segments,'complete':abs(supported-total)<1e-6}

def run():
    out={'scope':'Process features only; no native camera outcomes or target errors read','definition':'Both bounding native profiles admissible; gap <=1800 s; no interpolation across missing/invalid samples. |Ri|<=.03 is a nominal sensitivity screen, not proof of Monin-Obukhov validity. Measurement-resolved adds both wind/temperature slope signs bounded away from zero under supplied accuracy assumptions.','source_sha256':{},'sites':{}}
    for site in ('Tower1','Tower2'):
        p=ROOT/'turbulence'/f'{site}_NEUTRAL_PROFILE_FEATURES.csv';rows=list(csv.DictReader(p.open()));out['source_sha256'][str(p)]=hashlib.sha256(p.read_bytes()).hexdigest()
        labels=['2022-12-02T03:00:00','2022-12-16T03:00:00','2022-12-30T00:00:00' if site=='Tower1' else '2022-12-30T03:00:00','2023-01-13T03:00:00','2023-01-27T03:00:00','2023-02-10T03:00:00']
        out['sites'][site]=[{'start':a,'end':b,'role':'TRAIN' if k==0 else 'already_inspected_DEVELOPMENT','tiers':{tier:support(rows,stamp(a),stamp(b),tier)for tier in ('nominal_profile','near_neutral','measurement_resolved')}}for k,(a,b)in enumerate(zip(labels[:-1],labels[1:]))]
    out['decision']='NO-GO for complete native full-period flux-based melt prediction; no targets scored. A separately labelled reconstructed proxy would require a new missingness and stability specification.'
    dest=Path(__file__).parent/'NATIVE_SUPPORT_GATE.json';dest.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({s:[{k:round(v['supported_time_fraction'],4)for k,v in x['tiers'].items()}for x in rows]for s,rows in out['sites'].items()},indent=2))
if __name__=='__main__':run()
