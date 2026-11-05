"""Additional published equations on the same new stake windows; no refitting.

DETI variants are explicitly reconstructions: threshold and fractional lag
conventions were not recoverable as an exact author-code implementation.
"""
from pathlib import Path
from zipfile import ZipFile
from datetime import datetime
import json,math
import calendar
import numpy as np
from decode import raw_weather
from run_test import integrate,HOUR
R=Path(__file__).resolve().parent

def main():
    r=json.loads((R/'RESULTS.json').read_text());rows=[];refused=[]
    with ZipFile(R.parent/'Debris_properties_estimations.zip') as z:time,fields=raw_weather(z.read('Debris_properties_estimations/Data/Calibration/Meteo_AWS_filled_partition.mat'))
    for row in r['rows']:
        if row['method']!='robust_candidate':continue
        start=calendar.timegm(datetime.fromisoformat(row['date']).timetuple())*1000;end=start+24*HOUR;h=row['h_m'];label_end=row['label_end']
        lag=max(0,int(math.floor(21.54*h-1.193+.5)))
        try:lagged,weights=integrate(time,fields,start-lag*HOUR,end-lag*HOUR,label_end)
        except ValueError as e:
            refused.append(dict(site=row['site'],date=row['date'],label_end=label_end,reason=str(e)));continue
        current,w=integrate(time,fields,start,end,label_end);assert np.allclose(w,weights)
        issued={'KP1_2026_S22_fixed_effects':10**(.6215-1.4572*h)*max(row['mean_daily_T'],0)/917}
        tf=.016*h**(-.621);srf=.0079*math.exp(-11.21*h)
        rate=tf*lagged[:,0]+srf*.87*lagged[:,1]
        for threshold in [0,1]:
            for gate in ['original_current_gate','modified_lagged_gate']:
                active=(current[:,0] if gate=='original_current_gate' else lagged[:,0])>threshold
                issued[f'DETI_reconstruction_{gate}_TT{threshold}']=float(weights@(rate*active)/917)
        for method,p in issued.items():rows.append(dict(row,method=method,predicted_m_ice=p,error_m=p-row['observed_m_ice']))
        # Retain comparator/candidate rows on exactly this lag-supported subset.
        for other in r['rows']:
            if (other['site'],other['date'],other['label_end'])==(row['site'],row['date'],row['label_end']):rows.append(other)
    scores=[]
    for label_end in [False,True]:
        for method in sorted({x['method'] for x in rows}):
            sub=[x for x in rows if x['method']==method and x['label_end']==label_end];sites={}
            for site in [1,2,3]:
                group=[q for q in sub if q['site']==site]
                if group:sites[str(site)]=dict(n=len(group),MAE_m=float(np.mean([abs(q['error_m']) for q in group])))
            if sites:scores.append(dict(method=method,label_end=label_end,sites=sites,equal_site_MAE_m=float(np.mean([s['MAE_m'] for s in sites.values()])),pooled_day_MAE_m=float(np.mean([abs(q['error_m']) for q in sub]))))
    result=dict(scope='new2022-2023direct stake campaign; matched lag-supported daily windows; published fixed coefficients and explicitly labelled DETI reconstruction variants',sources={'KP1':'Winter-Billington2020Table4, DCG-MIP2026SupplementS22','DETI':'Carenzo2016Eq5-8 and section4.3; DCG-MIP2026S28-S31'},parameters={'lag_slope':21.54,'lag_intercept':-1.193,'TF1':.016,'TF2':-.621,'SRF1':.0079,'SRF2':-11.21,'original_albedo':.13,'TT_scenarios':[0,1],'lag_rounding':'nonnegative nearest integer, ties up'},caveat='DETI reproduction is not author-code equivalence; do not call it an exact reproduced operational model. No scenario selected by score.',rows=rows,scores=scores,refused=refused)
    with (R/'PUBLISHED_RESULTS.json').open('x') as f:json.dump(result,f,indent=2,allow_nan=False)
    print(json.dumps([{k:v for k,v in s.items() if k!='sites'} for s in scores],indent=2))

if __name__=='__main__':main()
