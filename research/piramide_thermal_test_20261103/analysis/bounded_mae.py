"""Conservative bounded-error certificates for a locked MAE comparison.

Elementary triangle-inequality bound. No historical theorem claim, probability
model, independence assumption, or automatic source-error certification.
"""
import math
def bounded_loss_difference(residual_a,residual_b,coeff_a,coeff_b,absolute_parameter_bounds):
    """Residuals r_m(e)=r_m(0)+A_m e, |e_j|<=delta_j.

    The reverse triangle inequality gives | |r+u|-|r| |<=|u|. Summing
    both model losses bounds the change of their MAE difference. This remains
    valid with shared endpoint errors, where a corner-only search may fail.
    """
    n=len(residual_a);bounds=list(absolute_parameter_bounds)
    if n<1 or len(residual_b)!=n or len(coeff_a)!=n or len(coeff_b)!=n:raise ValueError('inconsistent residual rows')
    if not bounds or any(type(x) not in (int,float) or not math.isfinite(x) or x<0 for x in bounds):raise ValueError('finite nonnegative conditional bounds required')
    for matrix in (coeff_a,coeff_b):
        if any(len(row)!=len(bounds) or any(not math.isfinite(x) for x in row) for row in matrix):raise ValueError('invalid error coefficients')
    if any(not math.isfinite(x) for x in list(residual_a)+list(residual_b)):raise ValueError('nonfinite residual')
    nominal=math.fsum(abs(a)-abs(b) for a,b in zip(residual_a,residual_b))/n
    radius=math.fsum((abs(a)+abs(b))*d for aa,bb in zip(coeff_a,coeff_b) for a,b,d in zip(aa,bb,bounds))/n
    lo,hi=nominal-radius,nominal+radius
    return dict(nominal_candidate_minus_baseline_mae=nominal,
        certified_conservative_difference_interval=[lo,hi],
        conclusion='candidate_strictly_better_throughout_box' if hi<0 else 'baseline_strictly_better_throughout_box' if lo>0 else 'unresolved_by_conservative_bound',
        exact_extrema_claim=False,confidence_interval=False,source_bounds_independently_verified_by_this_function=False)

def issued_forecast_bound(prediction_a,prediction_b,observations,per_reading_absolute_bound):
    n=len(observations)
    if len(prediction_a)!=n or len(prediction_b)!=n:raise ValueError('forecast support mismatch')
    # One common sequence of n+1 uncertain survey readings; each observed
    # increment uses the same shared middle reading in neighboring periods.
    matrix=[]
    for i in range(n):
        row=[0.]*(n+1);row[i]=1.;row[i+1]=-1.;matrix.append(row)
    return bounded_loss_difference([p-y for p,y in zip(prediction_a,observations)],
        [p-y for p,y in zip(prediction_b,observations)],matrix,matrix,[per_reading_absolute_bound]*(n+1))
