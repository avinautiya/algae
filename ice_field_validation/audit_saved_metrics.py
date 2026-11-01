"""Independent saved-event arithmetic audit; never fits or reruns a model."""
import hashlib
import json
import math
from pathlib import Path
import statistics

SOURCE = Path(__file__).resolve().parent/'event_benchmark'
DEST = Path(__file__).resolve().parent
BASELINE = 'training_pooled_constant_normal_rate'


def metrics(events, model):
    errors = [r['predicted_duration_h'][model]-r['observed_duration_h'] for r in events]
    absolute = list(map(abs, errors))
    relative = [a/r['observed_duration_h'] for a,r in zip(absolute, events)]
    n = len(events)
    return dict(MAE_h=math.fsum(absolute)/n, median_absolute_error_h=statistics.median(absolute),
                RMSE_h=math.sqrt(math.fsum(e*e for e in errors)/n), bias_h=math.fsum(errors)/n,
                MAPE=math.fsum(relative)/n, fraction_within_2h=sum(a<=2 for a in absolute)/n,
                fraction_within20percent=sum(a<=.2 for a in relative)/n, censored=0)


def audit(block, expected):
    events=block['events']; models=sorted(events[0]['predicted_duration_h'])
    assert len(events)==expected and len({r['id'] for r in events})==expected
    assert all(set(r['predicted_duration_h'])==set(models) for r in events)
    assert all(isinstance(r['observed_duration_h'], (int,float)) and math.isfinite(r['observed_duration_h']) and r['observed_duration_h']>0 for r in events)
    assert all(isinstance(v,(int,float)) and math.isfinite(v) and v>0 for r in events for v in r['predicted_duration_h'].values())
    sites=sorted({r['site'] for r in events}); recomputed={}; differences=[]
    for site in sites+['pooled']:
        selected=[r for r in events if site=='pooled' or r['site']==site]
        assert block['metrics'][site]['eligible']==len(selected)
        assert block['metrics'][site]['paired_complete']==len(selected)
        calculated={model:metrics(selected,model) for model in models}
        for model,values in calculated.items():
            values['MAE_skill_vs_training_constant']=1-values['MAE_h']/calculated[BASELINE]['MAE_h']
            for key,saved in block['metrics'][site]['models'][model].items():
                actual=values[key]; delta=abs(saved-actual)
                differences.append(delta)
                assert math.isclose(saved, actual, rel_tol=1e-12, abs_tol=1e-12),(site,model,key,saved,actual)
        recomputed[site]=dict(n=len(selected),models=calculated)
    equal_site={model:dict(MAE_h=statistics.mean(recomputed[s]['models'][model]['MAE_h'] for s in sites),
                          bias_h=statistics.mean(recomputed[s]['models'][model]['bias_h'] for s in sites),
                          RMSE_h=math.sqrt(statistics.mean(recomputed[s]['models'][model]['RMSE_h']**2 for s in sites))) for model in models}
    return dict(complete_events=len(events),models=models,no_model_censoring=True,identical_complete_event_support=True,
                exclusions_retained=len(block['exclusions']),site_weights={s:recomputed[s]['n']/len(events) for s in sites},
                maximum_absolute_saved_metric_difference=max(differences),recomputed=recomputed,equal_site_weighted=equal_site)


def main():
    first=json.loads((SOURCE/'FIELD_RESULTS.json').read_text())
    second=json.loads((SOURCE/'FIELD_RESULTS_2023.json').read_text())
    cal=json.loads((SOURCE/'CALIBRATION_2021.json').read_text())
    assert first['calibration']==cal==second['original_melt_calibration']
    assert cal['training_year']==2021 and first['test_year']==2022
    assert all(r['start_utc'].startswith('2022-') for r in first['events'])
    for year, filename, calkey in ((2022,'protocol.json','protocol_sha256'),(2023,'protocol2023.json',None)):
        digest=hashlib.sha256((SOURCE/filename).read_bytes()).hexdigest()
        if calkey: assert cal[calkey]==digest
        else: assert json.loads((SOURCE/'WEATHER_TRANSFER_FROZEN.json').read_text())['protocol_sha256']==digest
    result=dict(status='PASS: saved-event metrics independently recomputed; no model rerun or retuning',
                input_sha256={name:hashlib.sha256((SOURCE/name).read_bytes()).hexdigest() for name in ('FIELD_RESULTS.json','FIELD_RESULTS_2023.json','CALIBRATION_2021.json','protocol.json','protocol2023.json','WEATHER_TRANSFER_FROZEN.json')},
                year2022=audit(first,105), year2023={name:audit(block,71) for name,block in second['scenarios'].items()},
                calibration_check=dict(identical_saved_melt_calibration=True,training_year=2021,training_intervals=cal['training_intervals'],
                    training_sites=cal['training_sites'],temporal_holdout_years=[2022,2023],spatial_holdout=False,
                    weather_transfer_training_years=[2021,2022],weather_transfer_test_year=2023,
                    limits='Saved calibration and protocol agree; this arithmetic audit does not independently reconstruct fitting from raw observations or establish preregistration.'),
                interpretation=['Pooled metrics weight events equally, not sites equally; adjacent events share endpoints and weather.',
                    '2022 has two site-year histories; 2023 has four, with repeated sites. These are not 176 independent glaciers or validations.',
                    '2023 weather transfer uses 2021+2022 meteorology, not 2023 melt fitting; 2022 results were already inspected before secondary 2023 evaluation.',
                    'All results are retrospective observed-weather atmospheric-component hindcasts, not issued forecasts.',
                    'Missing hydrological heat is unknown, not measured zero. These models do not validate the new coupled hydraulic-memory/debris/stream system.',
                    'SW +/-15 percent are engineering scenarios, not uncertainty intervals. No confidence interval or accuracy guarantee is inferred.',
                    'Saved events retain no 2023 timestamps; 2023 temporal membership relies on protocol and reviewed evaluation code rather than event-level dates.'])
    (DEST/'METRIC_AUDIT.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    lines=['# Independent saved-event metric audit','',result['status']+'.','',
           'Recomputed every reported metric directly from saved observed and predicted event durations using Python standard-library arithmetic. No model was rerun, selected, or retuned. Input hashes and full recalculations are in `METRIC_AUDIT.json`.','',
           '| Test | Complete paired events | Site counts | Largest metric difference |','|---|---:|---|---:|']
    for label,block in [('2022',result['year2022'])]+[('2023 '+k,v) for k,v in result['year2023'].items()]:
        counts=', '.join(f"{s}: {v['n']}" for s,v in block['recomputed'].items() if s!='pooled')
        lines.append(f"| {label} | {block['complete_events']} | {counts} | {block['maximum_absolute_saved_metric_difference']:.3g} |")
    lines += ['', 'All four models have finite positive predictions on every eligible saved event: 105/105 in 2022 and 71/71 in every 2023 scenario. There is no model-specific censoring or complete-case selection within these eligible sets. The original eligibility exclusions remain; this audit does not imply all raw stake rows were eligible.', '',
              'Pooled scores weight each event equally. Site weights and alternative equal-site MAE, bias, and square-root mean site MSE are reported separately in JSON; averaging site RMSE is not the pooled RMSE. Pooled medians likewise cannot be averaged from site medians.', '',
              'Melt calibration dictionaries are identical across saved 2022 and 2023 results and the frozen 2021 calibration (94 intervals). Holdout is by year, not wholly new sites: training includes E/N/S/W, 2022 tests E/W, and 2023 tests E/N/W/W2. Frozen protocol hashes match. This verifies saved artifacts, not independent reconstruction of calibration history or preregistration. Weather transfer is fit on 2021–2022 paired weather; 2023 is a secondary held-out melt year after 2022 scores were inspected.', '']
    lines += ['- '+item for item in result['interpretation']]
    (DEST/'METRIC_AUDIT.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps(dict(status=result['status'],counts={'2022':105,'2023':71},max_differences={'2022':result['year2022']['maximum_absolute_saved_metric_difference'],**{k:v['maximum_absolute_saved_metric_difference'] for k,v in result['year2023'].items()}}),indent=2))


if __name__=='__main__':main()
