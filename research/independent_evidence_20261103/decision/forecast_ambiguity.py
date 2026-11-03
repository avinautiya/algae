"""Exact finite-heat observational-equivalence fixture, not a new theorem.

One training delivered-power observation does not identify h and inlet T.
All values are synthetic, not inferred Swiss sensor metadata.
"""
import math

def power(F,A,h,T):
 if any(type(x) not in (float,int) or not math.isfinite(x) or x<=0 for x in [F,A,h,T]):raise ValueError('positive finite arguments required')
 return F*T*(-math.expm1(-h*A/F))

def equivalent_temperature(P,F,A,h):
 if type(P) not in (float,int) or not math.isfinite(P) or P<=0:raise ValueError('positive delivered power required')
 return P/power(F,A,h,1.)

def demonstration():
 training_F=1000.;training_P=100.;A=1.;future_F=100.
 pairs=[]
 for h in [100.,300.]:
  T=equivalent_temperature(training_P,training_F,A,h)
  pairs.append(dict(h_W_m2_K=h,inlet_T_K=T,training_power_W=power(training_F,A,h,T),future_power_W=power(future_F,A,h,T)))
 return dict(scope='synthetic exact finite-heat equivalence; not measured field uncertainty or historical novelty',training_F_W_K=training_F,future_F_W_K=future_F,training_P_W=training_P,area_m2=A,pairs=pairs,future_power_difference_W=abs(pairs[0]['future_power_W']-pairs[1]['future_power_W']),consequence='Identical training power need not imply identical predictions after a flow change. Extra parameters do not supply missing measurements.')

if __name__=='__main__':
 import json
 print(json.dumps(demonstration(),indent=2,allow_nan=False))
