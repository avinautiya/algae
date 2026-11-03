"""Bounded-error ranking for two locally calibrated ablation predictors.

Both predict f*training_ablation. Consecutive calibration/evaluation intervals
share a reading, so their measurement errors cannot be varied independently
under a per-reading uncertainty interpretation. No novelty priority is claimed.
"""
import math
def finite(v):
 try:return type(v) in (int,float) and math.isfinite(v)
 except OverflowError:return False
def rank_bound(train_ablation,test_ablation,f_a,f_b,delta,*,association):
 if not all(finite(v) and v>=0 for v in [train_ablation,test_ablation,f_a,f_b,delta]):raise ValueError('finite nonnegative native-unit values required')
 train_radius=2*delta if association=='each_shared_reading' else delta if association=='each_interval' else None
 if train_radius is None:raise ValueError('explicit labelled uncertainty interpretation required')
 if train_ablation<=train_radius:raise ValueError('positive true calibration ablation not established')
 total=f_a+f_b;midpoint_residual=total*train_ablation-2*test_ablation
 radius=(2*total+4)*delta if association=='each_shared_reading' else (total+2)*delta
 lo=midpoint_residual-radius;hi=midpoint_residual+radius
 # ea²-eb² = train*(fa-fb)*[(fa+fb)*train-2*test].
 if f_a==f_b:status='identical_predictions'
 elif lo<=0<=hi:status='ranking_not_identified_under_bounds'
 elif (f_a-f_b)*(lo if lo>0 else hi)<0:status='a_has_lower_absolute_error_throughout_bounds'
 else:status='b_has_lower_absolute_error_throughout_bounds'
 return dict(nominal_prediction_a=f_a*train_ablation,nominal_prediction_b=f_b*train_ablation,nominal_error_a=abs(f_a*train_ablation-test_ablation),nominal_error_b=abs(f_b*train_ablation-test_ablation),midpoint_residual=midpoint_residual,midpoint_residual_bounds=[lo,hi],status=status,uncertainty_association=association,statistical_confidence_interval=False)
