"""One-gain cumulative-dose calibration using training ablation only.

This is an established empirical calibration, not independent measurement of
heat capacity or proof of a novel physical law. No evaluation outcome is input.
"""
import hashlib,json,math

def finite(x,name):
    if type(x) not in (int,float) or not math.isfinite(x):raise ValueError(name+' must be finite numeric')
    return float(x)

def fit_training_gain(intervals,*,train_end_s,evaluation_start_s,dose_quantity):
    """gain=sum(native ablation)/sum(model dose), preserving signed increments.

    Each interval requires start_s,end_s,observed_mmwe,dose. Gaps are explicit
    unsupported windows, not implicitly filled. Intervals must not overlap.
    """
    end=finite(train_end_s,'train end');ev=finite(evaluation_start_s,'evaluation start')
    if ev<end or not isinstance(dose_quantity,str) or not dose_quantity:raise ValueError('invalid split/quantity')
    if not intervals:raise ValueError('training intervals required')
    rows=[];previous=None
    for row in intervals:
        r={k:finite(row[k],k) for k in ('start_s','end_s','observed_mmwe','dose')}
        if r['end_s']<=r['start_s'] or r['end_s']>end or (previous is not None and r['start_s']<previous):
            raise ValueError('invalid, overlapping or evaluation-contaminated training intervals')
        if r['dose']<0:raise ValueError('negative dose unsupported; no silent clipping')
        rows.append(r);previous=r['end_s']
    total_dose=math.fsum(r['dose'] for r in rows);total_y=math.fsum(r['observed_mmwe'] for r in rows)
    if total_dose<=0 or total_y<=0:raise ValueError('positive integrated training dose and ablation required')
    gain=total_y/total_dose
    return dict(gain_mmwe_per_dose=gain,dose_quantity=dose_quantity,training_intervals=len(rows),
                training_total_mmwe=total_y,training_total_dose=total_dose,
                training_sha256=hashlib.sha256(json.dumps(rows,sort_keys=True,separators=(',',':')).encode()).hexdigest(),
                selection='one nonnegative cumulative-dose gain from all declared training intervals',
                parameter_is_independently_measured=False,uses_evaluation_targets=False,
                interpretation='effective predictive gain; not automatically physical conductivity or heat capacity')

def predict_from_frozen_gain(calibration,intervals,*,evaluation_start_s):
    """No observation field is accepted. Emit predictions before target release."""
    gain=finite(calibration['gain_mmwe_per_dose'],'gain');start=finite(evaluation_start_s,'evaluation start')
    if gain<=0:raise ValueError('positive trained gain required')
    result=[];previous=None
    for row in intervals:
        if set(row)!= {'start_s','end_s','dose'}:raise ValueError('prediction inputs must contain times/dose only; no outcomes')
        a,b,q=[finite(row[k],k) for k in ('start_s','end_s','dose')]
        if a<start or b<=a or q<0 or (previous is not None and a<previous):raise ValueError('invalid evaluation interval')
        result.append(dict(start_s=a,end_s=b,predicted_mmwe=gain*q));previous=b
    return result
