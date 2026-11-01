"""Operational recovery classification with explicit finite observation scope.

Sampled series use piecewise-linear interpolation. Event brackets remain the
adjacent observation times, not statistical confidence intervals. A final wet
sample never waives the declared persistence requirement.
"""
import hashlib
import json
import math


def finite(x,name):
    if isinstance(x,bool) or not isinstance(x,(int,float)) or not math.isfinite(x):
        raise ValueError(name+' must be finite numeric')
    return float(x)


def classify(days, fractions, *, theta=1e-3, persistence_days=1.,
             actual_removal_day=None, requested_removal_day=None,
             assessment_start_day=None, expected_horizon_day=None):
    t=[finite(x,'time') for x in days];f=[finite(x,'fraction') for x in fractions]
    if len(t)!=len(f) or len(t)<2 or any(a>=b for a,b in zip(t,t[1:])):
        raise ValueError('equal arrays of >=2 strictly increasing times required')
    theta=finite(theta,'threshold');persistence=finite(persistence_days,'persistence')
    if not 0<=theta<=1 or persistence<=0:raise ValueError('invalid operational definition')
    # Tiny numerical negative flows are permitted; the operational threshold
    # handles them. No clipping or presumed statistical uncertainty is applied.
    removal=None if actual_removal_day is None else finite(actual_removal_day,'actual removal')
    requested=None if requested_removal_day is None else finite(requested_removal_day,'requested removal')
    if removal is not None and not t[0]<=removal<=t[-1]:raise ValueError('actual removal outside observed interval')
    if assessment_start_day is None:start=removal if removal is not None else t[0]
    else:start=finite(assessment_start_day,'assessment start')
    if start<t[0] or start>t[-1] or (removal is not None and start<removal):
        raise ValueError('assessment must be observed and not precede actual removal')
    horizon=t[-1] if expected_horizon_day is None else finite(expected_horizon_day,'horizon')
    if horizon<t[-1]:raise ValueError('expected horizon precedes observations')
    above=[x>=theta for x in f]
    transitions=[]
    for i in range(1,len(t)):
        if above[i]!=above[i-1]:
            crossing=t[i-1]+(t[i]-t[i-1])*(theta-f[i-1])/(f[i]-f[i-1])
            transitions.append(dict(kind='wet' if above[i] else 'dry',
                                    estimate_day=crossing,bracket_day=[t[i-1],t[i]]))
    wet=[];wet_start=t[0] if above[0] else None;wet_bracket=[t[0],t[0]] if above[0] else None
    for tr in transitions:
        if tr['kind']=='wet':wet_start=tr['estimate_day'];wet_bracket=tr['bracket_day']
        elif wet_start is not None:
            wet.append(dict(start_day=wet_start,end_day=tr['estimate_day'],
                            start_bracket_day=wet_bracket,end_bracket_day=tr['bracket_day'],
                            right_censored=False));wet_start=None
    if wet_start is not None:
        wet.append(dict(start_day=wet_start,end_day=t[-1],start_bracket_day=wet_bracket,
                        end_bracket_day=[t[-1],t[-1]],right_censored=True))
    losses=[tr for tr in transitions if tr['kind']=='dry']
    loss=losses[0] if losses else None
    prior_loss=removal is not None and loss is not None and loss['estimate_day']<removal
    confirmed_prior_loss=prior_loss and loss['bracket_day'][1]<=removal
    candidates=[];preintervention_rewet=[]
    if prior_loss:
        for interval in wet:
            if interval['start_day']<=loss['estimate_day']:continue
            if interval['start_bracket_day'][1]<removal:
                preintervention_rewet.append(interval)
                continue
            beginning=max(interval['start_day'],start)
            estimated_duration=max(0.,interval['end_day']-beginning)
            # Conservative duration within the reported crossing brackets.
            confirmed_start=max(interval['start_bracket_day'][1],start)
            confirmed_end=interval['end_bracket_day'][0]
            confirmed_duration=max(0.,confirmed_end-confirmed_start)
            if estimated_duration>=persistence:
                candidates.append(dict(start_estimate_day=beginning,
                    onset_bracket_day=interval['start_bracket_day'],
                    onset_order_relative_to_removal=('straddles_removal' if interval['start_bracket_day'][0]<removal else 'at_or_after_removal'),
                    estimated_persistence_days=estimated_duration,
                    bracket_confirmed_persistence_days=confirmed_duration,
                    persistence_confirmed_by_brackets=confirmed_duration>=persistence))
    truncated=t[-1]<horizon
    if not above[0]:status='initially_dry'
    elif prior_loss and not confirmed_prior_loss:status='unresolved_loss_order'
    elif prior_loss and candidates:
        status='genuine_reconnection' if candidates[0]['persistence_confirmed_by_brackets'] else 'unresolved_persistence_bracket'
    elif prior_loss and preintervention_rewet and above[-1]:status='reconnected_before_intervention'
    elif truncated:status='unresolved_truncated'
    elif loss is None:status='maintained_contact'
    elif removal is None:status='lost_contact_no_intervention'
    elif not prior_loss:status='loss_after_intervention'
    elif t[-1]-start<persistence:status='unresolved_short_followup'
    elif above[-1]:status='unresolved_terminal_wet_episode'
    else:status='failed_recovery_by_observed_horizon'
    return dict(classification=status,initially_wet=above[0],loss=loss,
                loss_before_actual_removal=prior_loss,wet_intervals=wet,
                loss_before_actual_removal_confirmed_by_bracket=bool(confirmed_prior_loss),
                sustained_recovery_candidates=candidates,
                preintervention_rewet_intervals=preintervention_rewet,
                actual_removal_day=removal,requested_removal_day=requested,
                removal_delay_days=(None if removal is None or requested is None else removal-requested),
                assessment_start_day=start,observed_end_day=t[-1],expected_horizon_day=horizon,
                truncated=truncated,threshold=theta,persistence_days=persistence,
                scope='sampled operational contact; piecewise-linear interpolation, brackets not confidence intervals')


def history_key(*,initial_state,physics,boundary,forcing_history,interventions=(),event_day):
    """Group identical causal input histories, ignoring names/future events/grid.

    Caller supplies actual physical geometry and closure parameters. A hash is
    not proof of statistical independence or of equal computed trajectories.
    Numerics/mesh should be recorded separately, not counted as new physics.
    Each forcing/intervention item must have a time_day field.
    """
    event=finite(event_day,'event_day')
    def prefix(items):
        result=[]
        for item in items:
            if not isinstance(item,dict) or 'time_day' not in item:raise ValueError('history item needs time_day')
            if finite(item['time_day'],'history time')<event:result.append(item)
        return sorted(result,key=lambda r:r['time_day'])
    payload=dict(initial_state=initial_state,physics=physics,boundary=boundary,
                 forcing_history=prefix(forcing_history),interventions=prefix(interventions))
    raw=json.dumps(payload,sort_keys=True,separators=(',',':'),allow_nan=False).encode()
    return hashlib.sha256(raw).hexdigest()


def classify_run(run, *, theta=1e-3,persistence_days=1.,assessment_start_day=None):
    """Adapter for actual cases.run JSON; requires actual export event records."""
    if not isinstance(run,dict) or not isinstance(run.get('records'),list):raise ValueError('run records required')
    requested=run.get('removal_day');events=run.get('export',[])
    if not isinstance(events,list) or len(events)>1:raise ValueError('expected at most one actual clearing export')
    if requested is not None and not events:raise ValueError('requested removal lacks actual export timestamp')
    actual=None
    if events:
        if not isinstance(events[0],dict) or 't' not in events[0]:raise ValueError('actual export time t in seconds required')
        actual=finite(events[0]['t'],'actual export time')/86400.
    result=classify([r['day'] for r in run['records']],[r['f_toe'] for r in run['records']],
        theta=theta,persistence_days=persistence_days,actual_removal_day=actual,
        requested_removal_day=requested,assessment_start_day=assessment_start_day,
        expected_horizon_day=run.get('days'))
    result['actual_removal_source']='export[0].t / 86400' if events else 'no intervention'
    return result


if __name__=='__main__':
    import argparse
    from pathlib import Path
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('input',type=Path);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--theta',type=float,default=1e-3);p.add_argument('--persistence-days',type=float,default=1.)
    args=p.parse_args()
    if args.input.resolve()==args.output.resolve():p.error('output must not overwrite input')
    try:
        raw=json.loads(args.input.read_text(),parse_constant=lambda x:(_ for _ in ()).throw(ValueError('invalid JSON '+x)))
        answer=classify_run(raw,theta=args.theta,persistence_days=args.persistence_days)
        args.output.write_text(json.dumps(answer,indent=2,allow_nan=False)+'\n')
    except (ValueError,KeyError,TypeError,OSError) as exc:p.exit(2,'classification refused: '+str(exc)+'\n')
