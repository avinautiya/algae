"""Independent arithmetic and training-membership audit of saved predictions."""
import argparse, hashlib, json, math
from pathlib import Path

def main(path):
    d=json.loads(path.read_text());maxdiff=0.;completed=0;nominal=None
    for name,s in d['scenarios'].items():
        if s['status']!='completed':continue
        completed+=1;rows=s['accepted_intervals'];train=[r for r in rows if r['role']=='train'];test=[r for r in rows if r['role']=='test']
        assert set(s['calibration']['training_ids'])=={r['id'] for r in train}
        assert not ({r['id'] for r in train}&{r['id'] for r in test})
        for site,models in s['test_metrics'].items():
            selected=[r for r in test if site=='pooled' or r['station']==site]
            for model,m in models.items():
                es=[r['predictions_m'][model]-r['observed_lowering_m'] for r in selected]
                calculated=dict(MAE_m=sum(abs(e) for e in es)/len(es),RMSE_m=(sum(e*e for e in es)/len(es))**.5,
                                bias_m=sum(es)/len(es),observed_total_m=sum(r['observed_lowering_m'] for r in selected),predicted_total_m=sum(r['predictions_m'][model] for r in selected))
                assert len(es)==m['n']
                for k,v in calculated.items():maxdiff=max(maxdiff,abs(v-m[k]));assert math.isclose(v,m[k],rel_tol=1e-12,abs_tol=1e-12)
        if name=='clock+0h_T+0.00C_D1.30m':nominal=s
    invariance=0.
    for clock in [-2,0,2]:
        for temp in [-.05,0.,.05]:
            scenarios=[d['scenarios'][f'clock{clock:+d}h_T{temp:+.2f}C_D{diameter:.2f}m'] for diameter in [1.08,1.30,1.52]]
            if any(s['status']!='completed' for s in scenarios):continue
            for a,b in zip(scenarios,scenarios[1:]):
                aa={r['id']:r for r in a['accepted_intervals']};bb={r['id']:r for r in b['accepted_intervals']};assert aa.keys()==bb.keys()
                for key in aa:
                    for model in ['constant_rate_fitted','constant_h_fitted','db04_scale_fitted']:
                        delta=abs(aa[key]['predictions_m'][model]-bb[key]['predictions_m'][model]);invariance=max(invariance,delta);assert delta<1e-12
    result=dict(status='PASS',completed_scenarios=completed,max_metric_difference=maxdiff,
                fitted_predictions_max_DH_invariance_difference_m=invariance,
                input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                note='Engineering scenarios are dependent; no statistical significance or field-superiority claim')
    path.with_name('ARITHMETIC_AUDIT.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('results',type=Path);a=p.parse_args();main(a.results)
