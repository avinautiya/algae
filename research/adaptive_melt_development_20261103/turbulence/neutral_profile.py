"""Conditional neutral, zero-displacement constant-flux profile diagnostics.

Upward-positive H/LE. No surfaceTs, roughness, melt or target calibration.
Not eddy covariance or validated turbulent-flux measurements.
"""
import math
HEIGHTS=(.5,1.,2.,2.7);KAPPA=.4;RD=287.05;CP=1005.;LV=2.5e6

def log_slope(values):
 if len(values)!=4 or not all(math.isfinite(v)for v in values):raise ValueError('four finite levels required')
 x=[math.log(z)for z in HEIGHTS];mx=sum(x)/4;my=sum(values)/4;sxx=sum((v-mx)**2 for v in x);weights=[(v-mx)/sxx for v in x]
 slope=sum(a*b for a,b in zip(weights,values));intercept=my-slope*mx;residual=sum((v-intercept-slope*z)**2 for z,v in zip(x,values));sst=sum((v-my)**2 for v in values)
 return slope,weights,None if sst==0 else 1-residual/sst

def specific_humidity(T_C,RH_percent,pressure_Pa):
 if not all(math.isfinite(v)for v in [T_C,RH_percent,pressure_Pa])or not -80<T_C<80 or not 0<=RH_percent<=100 or pressure_Pa<=10000:raise ValueError('unsupported physical input')
 e=RH_percent/100*611*math.exp(17.27*T_C/(237.3+T_C))
 if e>=pressure_Pa:raise ValueError('vapor pressure exceeds total')
 return .622*e/(pressure_Pa-.378*e)

def neutral_flux(wind,T_C,RH_percent,pressure_Pa):
 if len(T_C)!=4 or len(RH_percent)!=4 or any(u<0 for u in wind):raise ValueError('four nonnegative wind levels')
 bu,w,ru=log_slope(wind);bt,_,rt=log_slope(T_C)
 q=[specific_humidity(t,r,pressure_Pa)for t,r in zip(T_C,RH_percent)];bq,_,rq=log_slope(q)
 windradius=sum(abs(x)*e for x,e in zip(w,[.45,.45,.89,.89]));tradius=.35*sum(abs(x)for x in w)
 out={'u_log_slope_m_s':bu,'T_log_slope_K':bt,'q_log_slope_kg_kg':bq,'wind_profile_R2':ru,'temperature_profile_R2':rt,'humidity_profile_R2':rq,'u_slope_conditional_accuracy_radius':windradius,'T_slope_conditional_accuracy_radius':tradius,'wind_slope_sign_resolved_given_stated_accuracy':bu>windradius,'T_slope_sign_resolved_given_stated_accuracy':abs(bt)>tradius,'pressure_Pa':pressure_Pa,'neutral_zero_displacement_assumed':True,'actual_flux_validation':False}
 if bu<=0:return {**out,'status':'unsupported_nonpositive_wind_log_slope','H_up_W_m2':None,'LE_up_W_m2':None}
 ustar=KAPPA*bu;Tmean=sum(T_C)/4+273.15;qmean=sum(q)/4;rho=pressure_Pa/(RD*Tmean*(1+.61*qmean));tv=[(t+273.15)*(1+.61*qi)for t,qi in zip(T_C,q)];btv,_,_=log_slope(tv)
 return {**out,'status':'conditional_neutral_proxy_only','u_star_m_s':ustar,'air_density_kg_m3':rho,'H_up_W_m2':-rho*CP*KAPPA*ustar*bt,'LE_up_W_m2':-rho*LV*KAPPA*ustar*bq,'gradient_Ri_at1m_approx':9.81*btv/(sum(tv)/4*bu*bu)}
