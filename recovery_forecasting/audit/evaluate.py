"""Paired operational recovery forecasts with interval/right-censored events."""
import math
from .evidence import digest


def number(x):return isinstance(x,(int,float)) and not isinstance(x,bool) and math.isfinite(x)


def event_score(observation,predicted_s):
    if predicted_s is not None and (not number(predicted_s) or predicted_s<0):
        raise ValueError('event forecast must be nonnegative time or null (no event)')
    if not isinstance(observation,dict):raise ValueError('event observation must be an object')
    kind=observation.get('kind')
    low=observation.get('lower_s')
    if not number(low) or low<0:raise ValueError('event lower bound must be nonnegative')
    if kind=='interval':
        high=observation.get('upper_s')
        if not number(high) or high<low:raise ValueError('invalid observed event interval')
        if predicted_s is None:return dict(kind=kind,missed_observed_event=True,interval_distance_s=None)
        return dict(kind=kind,missed_observed_event=False,
                    interval_distance_s=max(low-predicted_s,0.,predicted_s-high),
                    qualification='distance to observed bracket, not error against a known exact event time')
    if kind=='right_censored':
        return dict(kind=kind,early_prediction_violation=predicted_s is not None and predicted_s<low,
                    qualification='event-free through lower_s; later event time is unobserved and has no MAE')
    raise ValueError('event kind must be interval or right_censored')


def evaluate(record):
    """Classifications refer to a common explicitly declared operational target.

    Case rows are not presumed independent. No confidence intervals or p-values
    are manufactured. Censored event forecasts are not treated as exact times.
    """
    if not isinstance(record,dict):raise ValueError('evaluation must be an object')
    fingerprint=digest(record)
    target=record.get('target',{})
    threshold=target.get('toe_fraction_threshold');window=target.get('persistence_s')
    if not number(threshold) or not 0<threshold<1 or not number(window) or window<=0:
        raise ValueError('explicit toe fraction threshold and positive persistence_s required')
    if target.get('after_ordinary_flow_return') is not True:
        raise ValueError('target must be assessed after ordinary-flow return')
    models=record.get('models');cases=record.get('cases')
    if not isinstance(models,list) or not models or len(set(models))!=len(models) or not all(isinstance(m,str) and m for m in models):
        raise ValueError('unique named models required')
    if not isinstance(cases,list) or not cases:raise ValueError('nonempty case list required')
    seen=set();scores={m:dict(tp=0,tn=0,fp=0,fn=0,abstentions=0,unresolved_labels=0,event_scores=[]) for m in models}
    paired={};labels={}
    for row in cases:
        case=row.get('case_id')
        if not isinstance(case,str) or not case or case in seen:raise ValueError('unique nonempty case_id required')
        seen.add(case)
        observed=row.get('recovered')
        if observed is not None and not isinstance(observed,bool):raise ValueError('recovered must be bool or null')
        predictions=row.get('predictions')
        if not isinstance(predictions,dict) or set(predictions)!=set(models):raise ValueError('every case requires exactly the declared models')
        labels[case]=observed;paired[case]={}
        for m in models:
            pred=predictions[m]
            if not isinstance(pred,dict):raise ValueError('prediction must be an object')
            value=pred.get('recovered')
            if value is not None and not isinstance(value,bool):raise ValueError('predicted recovery must be bool or null')
            if observed is None:scores[m]['unresolved_labels']+=1
            elif value is None:scores[m]['abstentions']+=1
            else:
                key=('tp' if value else 'fn') if observed else ('fp' if value else 'tn')
                scores[m][key]+=1;paired[case][m]=int(value==observed)
            if 'event' in row:
                e=event_score(row['event'],pred.get('event_time_s'));e['case_id']=case
                scores[m]['event_scores'].append(e)
    for m,s in scores.items():
        n=s['tp']+s['tn']+s['fp']+s['fn']
        s['scored_classifications']=n
        s['accuracy_on_nonabstained']=((s['tp']+s['tn'])/n if n else None)
        s['sensitivity']=(s['tp']/(s['tp']+s['fn']) if s['tp']+s['fn'] else None)
        s['specificity']=(s['tn']/(s['tn']+s['fp']) if s['tn']+s['fp'] else None)
        ds=[e['interval_distance_s'] for e in s['event_scores'] if e.get('interval_distance_s') is not None]
        s['mean_distance_to_observed_event_bracket_s']=sum(ds)/len(ds) if ds else None
        s['event_distance_count']=len(ds)
    comparisons=[]
    for i,a in enumerate(models):
        for b in models[i+1:]:
            diffs=[v[a]-v[b] for v in paired.values() if a in v and b in v]
            comparisons.append(dict(model_a=a,model_b=b,paired_cases=len(diffs),
                a_better=sum(d>0 for d in diffs),b_better=sum(d<0 for d in diffs),
                ties=sum(d==0 for d in diffs),accuracy_difference=sum(diffs)/len(diffs) if diffs else None))
    return dict(target=target,case_count=len(cases),scores=scores,paired_comparisons=comparisons,
                input_sha256=fingerprint,
                qualification='Descriptive paired scores only; no independence assumed, no confidence interval, no readiness approval')
