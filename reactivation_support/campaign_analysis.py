"""Descriptive event analysis for Claude's M3 run JSON; no production approval.

Threshold and persistence window must be declared by the experiment designer.
Piecewise-linear event estimates do not resolve unrecorded transients.
"""
from pathlib import Path
import argparse,json,math
import numpy as np
from scipy.integrate import trapezoid


def contact_intervals(times,fractions,threshold):
    t=np.asarray(times,dtype=float);f=np.asarray(fractions,dtype=float)
    if t.ndim!=1 or f.shape!=t.shape or len(t)<2 or not np.all(np.isfinite(t)) or not np.all(np.isfinite(f)) or np.any(np.diff(t)<=0):
        raise ValueError('finite matching series at strictly increasing times required')
    if not math.isfinite(threshold) or not 0<threshold<1:raise ValueError('threshold must lie in (0,1)')
    intervals=[];crossings=[];start=float(t[0]) if f[0]>=threshold else None
    for i in range(len(t)-1):
        w0,w1=f[i]>=threshold,f[i+1]>=threshold
        if w0==w1:continue
        cross=float(t[i]+(t[i+1]-t[i])*(threshold-f[i])/(f[i+1]-f[i]))
        crossings.append(dict(day=cross,direction='gain' if w1 else 'loss',observation_bracket_days=[float(t[i]),float(t[i+1])]))
        if w1:start=cross
        else:
            if start is not None:intervals.append([start,cross])
            start=None
    if start is not None:intervals.append([start,float(t[-1])])
    return intervals,crossings


def summarize(run,*,threshold,persistence_hours,ordinary_from_day=None):
    if not math.isfinite(persistence_hours) or persistence_hours<=0:raise ValueError('positive persistence window required')
    records=run['records'];t=[r['day'] for r in records];f=[r['f_toe'] for r in records]
    intervals,crossings=contact_intervals(t,f,threshold)
    losses=[x for x in crossings if x['direction']=='loss']
    first_loss=losses[0] if losses else None
    ordinary=run.get('removal_day') if ordinary_from_day is None else ordinary_from_day
    ordinary_source='explicit argument' if ordinary_from_day is not None else 'requested removal (constant forcing assumed)'
    if ordinary_from_day is None and ordinary is not None and run.get('export'):
        if len(run['export'])!=1:raise ValueError('multiple exports require an explicit ordinary-flow start')
        actual=run['export'][0].get('t')
        if actual is not None:
            ordinary=float(actual)/86400
            ordinary_source='actual export timestamp (constant forcing assumed)'
    if ordinary is not None and (not math.isfinite(ordinary) or ordinary<t[0] or ordinary>t[-1]):raise ValueError('ordinary-flow start outside observations')
    persistence=persistence_hours/24
    reconnect=None;followup_limited=False
    if ordinary is not None and first_loss is not None:
        for a,b in intervals:
            start=max(a,ordinary)
            if start<first_loss['day'] or start>b:continue
            if b-start>=persistence:
                reconnect=start;break
            if b==t[-1] and b>start:followup_limited=True
    final_active=f[-1]>=threshold
    wet_days=sum(b-a for a,b in intervals)
    row=dict(layout=run.get('layout'),dx=run.get('dx'),morphology_step_s=run.get('dtm'),
             observed_horizon_days=[t[0],t[-1]],requested_horizon_days=run.get('days'),
             threshold=threshold,persistence_hours=persistence_hours,ordinary_from_day=ordinary,
             ordinary_time_source=ordinary_source,requested_removal_day=run.get('removal_day'),
             initially_active=f[0]>=threshold,contact_intervals_days=intervals,
             total_observed_contact_days=wet_days,first_loss=first_loss,
             sustained_reconnection_start_day=reconnect,insufficient_followup=followup_limited,
             active_at_observed_horizon=final_active,max_sampling_interval_hours=float(max(np.diff(t))*24),
             qualification='descriptive interpolation; no independent physical validation or final-data eligibility asserted')
    if ordinary is not None:row['post_ordinary_contact_days']=sum(max(0,b-max(a,ordinary)) for a,b in intervals)
    for key,outkey in [('Q_toe','integrated_net_toe_water_m3'),('Q_byp','integrated_net_bypass_water_m3')]:
        if all(key in r for r in records):
            vals=np.array([r[key] for r in records],dtype=float)
            if not np.all(np.isfinite(vals)):raise ValueError('nonfinite flow series')
            row[outkey]=float(trapezoid(vals,t)*86400)
    for key in ('cum_melt_toe','cum_melt_byp'):
        if all(key in r for r in records):
            vals=np.array([r[key] for r in records],dtype=float)
            if not np.all(np.isfinite(vals)) or np.any(np.diff(vals)<-1e-12):raise ValueError('invalid cumulative melt')
            row[key+'_m3']=float(vals[-1]-vals[0])
    row['hydraulic_convergence_recorded']=all('hyd_converged' in r for r in records)
    row['all_recorded_hydraulic_solves_converged']=all(r['hyd_converged'] is True for r in records) if row['hydraulic_convergence_recorded'] else None
    row['requested_horizon_completed']=math.isclose(t[-1],run['days'],rel_tol=0,abs_tol=1e-9) if 'days' in run else None
    return row


def refinement_diagnostic(values,ratio=2.):
    if len(values)!=3 or any(v is None for v in values):return dict(status='three uncensored values required')
    if not all(math.isfinite(v) for v in values) or not math.isfinite(ratio) or ratio<=1:raise ValueError('invalid refinement data')
    a,b,c=values;d0=a-b;d1=b-c
    if d0*d1<=0:return dict(status='nonmonotone or degenerate; no Richardson extrapolation')
    p=math.log(abs(d0/d1))/math.log(ratio)
    if p<=0:return dict(status='differences not contracting; no asymptotic estimate',observed_order=p)
    estimate=abs(d1)/(ratio**p-1)
    return dict(status='conditional asymptotic diagnostic only',observed_order=p,
                extrapolated_value=c+(c-b)/(ratio**p-1),fine_grid_estimated_error=estimate,
                qualification='not a confidence interval; excludes model/reference mismatch and requires an asymptotic, consistent geometry sequence')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('inputs',nargs='+',type=Path)
    p.add_argument('--threshold',type=float,required=True);p.add_argument('--persistence-hours',type=float,required=True)
    p.add_argument('--ordinary-from-day',type=float);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    out=[]
    for f in args.inputs:
        result=summarize(json.loads(f.read_text()),threshold=args.threshold,persistence_hours=args.persistence_hours,ordinary_from_day=args.ordinary_from_day)
        out.append(dict(source=str(f),summary=result))
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(out,indent=2,allow_nan=False))
