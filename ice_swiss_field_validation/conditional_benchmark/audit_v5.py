"""Saved-output audit: arithmetic, local geometry and common-cohort feature effect."""
import argparse,hashlib,json,math
from pathlib import Path

MODELS=['constant_rate_fitted','constant_h_fitted','db04_unfitted','db03_unfitted','gnielinski_smooth_unfitted','db04_scale_fitted']

def mae(rows,model):
    return math.fsum(abs(r['predictions_m'][model]-r['observed_lowering_m']) for r in rows)/len(rows)

def main(v5root,v4root):
    summaries={};audits={}
    for method in ['linear','left','right']:
        path=v5root/method/'CONDITIONAL_RESULTS.json';data=json.loads(path.read_text())
        old=json.loads((v4root/method/'CONDITIONAL_RESULTS.json').read_text())
        metric_difference=geometry_difference=0.;count=0;comparisons={}
        for key,scenario in data['scenarios'].items():
            if scenario['status']!='completed':continue
            count+=1;rows=scenario['accepted_intervals'];train=[r for r in rows if r['role']=='train'];test=[r for r in rows if r['role']=='test']
            assert set(scenario['calibration']['training_ids'])=={r['id'] for r in train}
            assert not({r['id'] for r in train}&{r['id'] for r in test})
            for row in rows:
                for cell in row['hydraulic_cells']:
                    # Independent rearrangement: Re=4Q/[nu(w+2h)].
                    expected=4*cell['discharge_m3_s']/(1.8e-6*(cell['channel_width_m']+2*cell['local_depth_m']))
                    geometry_difference=max(geometry_difference,abs(cell['Re']-expected)/expected)
                    assert math.isclose(cell['Re'],expected,rel_tol=1e-12)
            for site,scores in scenario['test_metrics'].items():
                selected=[r for r in test if site=='pooled' or r['station']==site]
                for model,reported in scores.items():
                    errors=[r['predictions_m'][model]-r['observed_lowering_m'] for r in selected]
                    recalculated=dict(MAE_m=sum(abs(e) for e in errors)/len(errors),bias_m=sum(errors)/len(errors),RMSE_m=math.sqrt(sum(e*e for e in errors)/len(errors)))
                    for name,value in recalculated.items():
                        metric_difference=max(metric_difference,abs(value-reported[name]));assert math.isclose(value,reported[name],rel_tol=1e-12,abs_tol=1e-12)
            if math.isclose(scenario['channel_width_m'],2.):
                prior=old['scenarios'][key.replace('W2.00m','D1.30m')]['accepted_intervals']
                trainids={r['id'] for r in train};testids={r['id'] for r in test}
                priortrain=[r for r in prior if r['id'] in trainids]
                priortest=[dict(r,predictions_m=dict(r['predictions_m'])) for r in prior if r['id'] in testids]
                assert len(priortrain)==len(train) and len(priortest)==len(test)
                fitted={}
                for model,basis in [('constant_rate_fitted','duration_s'),('constant_h_fitted','constant_h_basis_m'),('db04_scale_fitted','db04_basis_m')]:
                    coefficient=max(0.,math.fsum(r[basis]*r['observed_lowering_m'] for r in priortrain)/math.fsum(r[basis]**2 for r in priortrain));fitted[model]=coefficient
                    for row in priortest:row['predictions_m'][model]=coefficient*row[basis]
                comparisons[key]=dict(common_train_ids=sorted(trainids),common_test_ids=sorted(testids),
                    v4_refit_same_training=dict(fitted=fitted,MAE_m={site:{model:mae([r for r in priortest if site=='pooled' or r['station']==site],model) for model in MODELS} for site in ['P3','P5','pooled']}),
                    v5_local_MAE_m={site:{model:scenario['test_metrics'][site][model]['MAE_m'] for model in MODELS} for site in ['P3','P5','pooled']})
        advantages={}
        for site in ['P3','P5','pooled']:
            advantages[site]={}
            for baseline in MODELS[:-1]:
                values=[s['test_metrics'][site][baseline]['MAE_m']-s['test_metrics'][site]['db04_scale_fitted']['MAE_m'] for s in data['scenarios'].values() if s['status']=='completed']
                advantages[site][baseline]=dict(improved=sum(v>0 for v in values),scenarios=len(values),min_advantage_m=min(values),max_advantage_m=max(values))
        audits[method]=dict(status='PASS',completed_scenarios=count,max_metric_difference=metric_difference,max_relative_Re_geometry_identity_difference=geometry_difference,source_sha256=hashlib.sha256(path.read_bytes()).hexdigest())
        summaries[method]=dict(engineering_advantages=advantages,common_training_and_test_geometry_comparison=comparisons)
    (v5root/'V5_ARITHMETIC_AND_GEOMETRY_AUDIT.json').write_text(json.dumps(audits,indent=2)+'\n')
    (v5root/'V5_COMMON_COHORT_FEATURE_AUDIT.json').write_text(json.dumps(summaries,indent=2)+'\n')
    print(json.dumps(audits,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('v5root',type=Path);p.add_argument('v4root',type=Path);a=p.parse_args();main(a.v5root,a.v4root)
