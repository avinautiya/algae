"""Necessary recovery-energy bounds, not an environmental forecast or novelty claim."""
import math
RHOI=917.;LF=334000.;RHOW=1000.;CP=4217.
def finite(x,name,nonnegative=True):
 if isinstance(x,bool) or not isinstance(x,(int,float)) or not math.isfinite(x) or (nonnegative and x<0):raise ValueError(name+' invalid')
 return x

def water_ice_flux(q_m3_s,area_m2,h_W_m2_K,superheat_K):
 for k,v in locals().copy().items():finite(v,k)
 if area_m2<=0:raise ValueError('positive area required')
 if q_m3_s==0:return 0.
 f=RHOW*CP*q_m3_s
 return f*superheat_K*(-math.expm1(-h_W_m2_K*area_m2/f))/area_m2

def screen(*,required_gap_repair_m,horizon_s,pulse_discharge_m3_s,bypass_area_m2,
 toe_h_upper_W_m2_K,toe_temperature_upper_K,bypass_h_lower_W_m2_K,bypass_temperature_lower_K,
 additional_toe_ice_flux_upper_W_m2,additional_bypass_ice_flux_lower_W_m2,
 evidence):
 """All bounds must apply throughout positive-gap interval; equal hydraulic K.

 Additional ports are heat actually delivered to horizontal ice area, not raw
 atmospheric radiation on a vertical cliff or above an insulating deposit.
 Negative external bypass heat is conservatively floored at zero total melt;
 no refreezing or cold-state prediction is asserted. Metadata cannot prove bounds.
 """
 for k,v in list(locals().items()):
  if k!='evidence':finite(v,k,not k.startswith('additional_'))
 if horizon_s<=0 or bypass_area_m2<=0:raise ValueError('positive time/area required')
 if not isinstance(evidence,dict) or not evidence.get('bound_sources') or evidence.get('equal_rating_controls_verified') is not True or evidence.get('bounds_hold_over_window') is not True:raise ValueError('explicit independently justified bound contract required')
 toe_upper=max(0.,toe_h_upper_W_m2_K*toe_temperature_upper_K+additional_toe_ice_flux_upper_W_m2)
 bypass_lower=max(0.,water_ice_flux(pulse_discharge_m3_s/2,bypass_area_m2,bypass_h_lower_W_m2_K,bypass_temperature_lower_K)+additional_bypass_ice_flux_lower_W_m2)
 rate=(toe_upper-bypass_lower)/(RHOI*LF)
 maximum_repair=max(0.,rate)*horizon_s
 return dict(required_gap_repair_m=required_gap_repair_m,maximum_repair_m=maximum_repair,
  repair_rate_upper_m_s=rate,ruled_out_within_window=maximum_repair<required_gap_repair_m,
  minimum_extra_toe_ice_flux_to_avoid_this_exclusion_W_m2=max(0.,required_gap_repair_m*RHOI*LF/horizon_s+bypass_lower-toe_h_upper_W_m2_K*toe_temperature_upper_K),
  evidence=evidence,interpretation='Necessary condition only. Passing means not excluded; it does not predict recovery, sustained access or field skill.')
