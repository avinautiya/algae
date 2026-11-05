"""Independent arithmetic audit, inverse unit checks and weighting diagnostics."""
from pathlib import Path
from zipfile import ZipFile
import json,hashlib,calendar
import numpy as np
from decode import table
R=Path(__file__).resolve().parent

def main():
    result=json.loads((R/'RESULTS.json').read_text());freeze=(R.parent/'PIRAMIDE_FREEZE.json').read_bytes()
    assert hashlib.sha256(freeze).hexdigest()==result['model_freeze_sha256'];models=json.loads(freeze)['models']
    assert all('Piramide' not in m['donors'] for m in models.values())
    native={}
    with ZipFile(R.parent/'Debris_properties_estimations.zip') as z:
        for site in [1,2,3]:
            t,c=table(z.read(f'Debris_properties_estimations/Data/Ablation_Tower{site}.mat'));native[site]=(t,c)
    checked=0
    for row in result['rows']:
        from datetime import datetime
        t,c=native[row['site']];start=calendar.timegm(datetime.fromisoformat(row['date']).timetuple())*1000;end=calendar.timegm(datetime.fromisoformat(row['end']).timetuple())*1000
        i=np.flatnonzero(t==start);j=np.flatnonzero(t==end);assert len(i)==len(j)==1
        # Verify native converted ablation independently against distance,
        # which also prevents a factor1000 README-unit error.
        obs=float((c['ablation'][j[0]]-c['ablation'][i[0]])/.917)
        assert abs(obs-row['observed_m_ice'])<1e-12
        if row['method'] in models:
            m=models[row['method']];p=sum(row['features'][f]/s*b for f,s,b in zip(m['features'],m['feature_scale'],m['beta']))*m['target_scale']
        else:p=10**(.62-1.46*row['h_m'])*max(row['mean_daily_T'],0)/917
        assert abs(p-row['predicted_m_ice'])<1e-12;checked+=1
    score_count=0
    for score in result['scores']:
        sub=[r for r in result['rows'] if r['method']==score['method'] and r['label_end']==score['label_end'] and (score['stratum']=='all_complete_days' or r['source_snow_screen'])]
        values=[]
        for site,s in score['sites'].items():
            v=[abs(r['predicted_m_ice']-r['observed_m_ice']) for r in sub if str(r['site'])==site]
            assert len(v)==s['n'] and abs(np.mean(v)-s['MAE_m'])<1e-12;values.append(np.mean(v))
        assert abs(np.mean(values)-score['equal_site_MAE_m'])<1e-12;score_count+=1
    diagnostics=[]
    for method in models.keys()|{'published_KP1_2020_fixed_effects'}:
        sub=[r for r in result['rows'] if r['method']==method and not r['label_end']]
        diagnostics.append(dict(method=method,n=len(sub),pooled_day_MAE_m=float(np.mean([abs(r['error_m']) for r in sub]))))
    matched=json.loads((R/'PUBLISHED_RESULTS.json').read_text());assert len(matched['rows'])>0
    for label_end in [False,True]:
        groups={}
        for method in {r['method'] for r in matched['rows']}:
            groups[method]={(r['site'],r['date']) for r in matched['rows'] if r['method']==method and r['label_end']==label_end}
        assert all(keys==groups['robust_candidate'] for keys in groups.values())
    out=dict(frozen_model_hash_verified=True,no_recipient_training=True,native_units_identity_checked=True,predictions_recomputed=checked,score_summaries_recomputed=score_count,matched_published_support_verified=True,pooled_day_diagnostic=diagnostics,interpretation='one external campaign with measured stakes; frozen coefficients, audited endpoint definitions; not an independently preregistered multi-glacier confirmation')
    (R/'AUDIT.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))

if __name__=='__main__':main()
