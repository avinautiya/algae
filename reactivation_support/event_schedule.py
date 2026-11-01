"""Exact event-aligned morphology intervals, independent of hydraulic stepping."""
import math


def schedule(horizon_s,nominal_step_s,events_s=()):
    if any(isinstance(x,bool) or not math.isfinite(x) or x<=0 for x in (horizon_s,nominal_step_s)):
        raise ValueError('positive finite horizon and nominal step required')
    events=[]
    for t in events_s:
        if isinstance(t,bool) or not math.isfinite(t) or not 0<=t<=horizon_s:
            raise ValueError('event outside simulation horizon')
        events.append(float(t))
    points={0.,float(horizon_s),*events}
    # Preserve event timestamps if a floating-point regular tick is coincident.
    aliases={}
    for e in events:
        k=round(e/nominal_step_s);tick=k*nominal_step_s
        if abs(e-tick)<=8*math.ulp(max(abs(e),abs(tick),1.)):
            aliases[k]=e
    k=1
    while k*nominal_step_s<horizon_s:
        points.add(aliases.get(k,k*nominal_step_s));k+=1
    times=sorted(points)
    return [(a,b,b-a) for a,b in zip(times[:-1],times[1:])]
