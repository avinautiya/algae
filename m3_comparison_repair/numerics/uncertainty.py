"""Conditional three-grid Richardson/GCI estimate, never a statistical CI."""
import json
import math


def estimate(values, spacings, safety_factor=1.25):
    if len(values)!=3 or len(spacings)!=3:
        raise ValueError('exactly three coarse-to-fine observations required')
    if not all(not isinstance(v,bool) and isinstance(v,(int,float)) and math.isfinite(v) for v in [*values,*spacings,safety_factor]):
        raise ValueError('finite numbers required')
    if not spacings[0]>spacings[1]>spacings[2]>0 or safety_factor<1:
        raise ValueError('strictly decreasing positive spacing and safety >=1 required')
    r=spacings[0]/spacings[1]
    if not math.isclose(r,spacings[1]/spacings[2],rel_tol=1e-10):
        raise ValueError('constant refinement ratio required')
    a,b,c=values;d1=a-b;d2=b-c
    base=dict(values=list(values),spacings=list(spacings),refinement_ratio=r,
              scope='conditional numerical extrapolation; not a confidence interval',
              assumes='same observable and problem, asymptotic single-power leading error; not established by three points')
    if d1==0 or d2==0:
        return dict(**base,status='not_estimable_zero_difference')
    if d1*d2<0:return dict(**base,status='not_estimable_oscillatory')
    if abs(d2)>=abs(d1):return dict(**base,status='not_estimable_noncontracting')
    p=math.log(abs(d1/d2))/math.log(r)
    extrapolated=c+(c-b)/math.expm1(p*math.log(r))
    error=abs(c-extrapolated);gci=safety_factor*error
    return dict(**base,status='conditional_estimate',observed_order=p,
                extrapolated_value=extrapolated,fine_grid_absolute_error_estimate=error,
                safety_factor=safety_factor,gci_absolute=gci,
                gci_relative=(gci/abs(c) if c!=0 else None),
                numerical_estimate_band=[c-gci,c+gci],
                next_halving_prediction=extrapolated+(c-extrapolated)/r**p)


if __name__=='__main__':
    print(json.dumps(estimate([3.165,2.721,2.478],[.1,.05,.025]),indent=2,allow_nan=False))
