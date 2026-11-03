"""Finite-family calibrated-goal diameter. Not a statistical confidence interval.
No glacier data or hidden camera outcomes are imported.
"""
import math

def certificate(models, training_target, tolerance):
    if not math.isfinite(training_target) or training_target <= 0:
        raise ValueError('positive finite training target required')
    if not math.isfinite(tolerance) or tolerance < 0:
        raise ValueError('finite nonnegative absolute tolerance required')
    rows=[]
    for model in models:
        d=float(model['training_dose']); f=float(model['future_dose'])
        if not math.isfinite(d) or d<=0 or not math.isfinite(f) or f<0:
            raise ValueError('positive training and nonnegative future doses required')
        rows.append({'id':model['id'], 'ratio':f/d,
                     'calibrated_prediction':training_target*f/d})
    if not rows: raise ValueError('empty feasible family')
    low=min(rows,key=lambda x:x['calibrated_prediction'])
    high=max(rows,key=lambda x:x['calibrated_prediction'])
    radius=(high['calibrated_prediction']-low['calibrated_prediction'])/2
    return {'family':rows,'lower':low['calibrated_prediction'],
            'upper':high['calibrated_prediction'],
            'minimax_point':(low['calibrated_prediction']+high['calibrated_prediction'])/2,
            'conditional_minimax_absolute_error':radius,
            'abstain':radius>tolerance,
            'witness_ids':[low['id'],high['id']],
            'scope':'Only finite declared feasible model family; no continuous-set or field-coverage guarantee.'}

def passive_rc_doses(tau,train_duration=1.,future_duration=1.):
    """Unit step during TRAIN, zero future: passive first-order low-pass flux.
    State initial0, y'= (u-y)/tau. Output is nonnegative.
    """
    if min(tau,train_duration,future_duration)<=0:raise ValueError('positive durations')
    train=train_duration-tau*(-math.expm1(-train_duration/tau))
    end=-math.expm1(-train_duration/tau)
    future=end*tau*(-math.expm1(-future_duration/tau))
    return train,future

if __name__=='__main__':
    import json
    rows=[]
    for tau in [.1,1.,10.]:
        t,f=passive_rc_doses(tau)
        rows.append({'id':f'passive_tau_{tau}','training_dose':t,'future_dose':f})
    print(json.dumps(certificate(rows,100.,5.),indent=2))
