"""Execute declared feature-only fit. No mass fields/imports/network access."""
import json,pathlib,hashlib,csv,datetime as dt,math
from .spectroscopy import harmonic_extract,bounded_fit,jacobian_rank,transfer
P=pathlib.Path(__file__).parent;DATA=P.parent/'data'/'kennicott'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 plan=json.loads((P/'FEATURE_ONLY_PLAN.json').read_text());raw=(DATA/'features'/'KDT_2020b.csv').read_bytes()
 if hashlib.sha256(raw).hexdigest()!=plan['source_sha256']:raise ValueError('Source feature hash mismatch')
 rows=list(csv.DictReader(raw.decode('utf-8-sig').splitlines()));a=dt.datetime.fromisoformat(plan['fit_start_UTC']);b=dt.datetime.fromisoformat(plan['fit_end_UTC']);r=[x for x in rows if a<=dt.datetime.fromisoformat(x['Datetime'])<=b];t=[(dt.datetime.fromisoformat(x['Datetime'])-a).total_seconds() for x in r];u=[float(x['Depth: 5 cm']) for x in r];y=[float(x['Depth: 15 cm']) for x in r]
 result={'plan_sha256':sha(P/'FEATURE_ONLY_PLAN.json'),'code_sha256':sha(P/'spectroscopy.py'),'run_code_sha256':sha(P/'run_feature_fit.py'),'mass_values_used':False,'source_zero_22cm_used':False,'native_rows':len(r),'native_cadence_s':sorted(set(v-w for w,v in zip(t,t[1:]))),'input_or_output_negative_count':sum(v<0 for v in u+y),'scope':'Effective periodic homogeneous Robin response, fixedgeometry; not independent material-contact identification.'}
 try:obs=harmonic_extract(t,u,y)
 except ValueError as e:
  result.update({'status':'NO_FIT_SOURCE_RESOLUTION_GATE','reason':str(e)});(P/'RESULTS_FEATURE_ONLY.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2));return
 result['harmonic_extraction']=obs
 for label,o in [('primary_DC_and_harmonics',obs['observations']),('diagnostic_harmonics_only',obs['observations'][1:])]:
  fit=bounded_fit(o);null=bounded_fit(o,null=True);fit['local_rank']=jacobian_rank(fit['alpha_m2_s'],fit['b_kR_m'],o);fit['conditional_capacity_resistance']=[{'C_J_m3_K':C,'k_W_m_K':fit['alpha_m2_s']*C,'R_m2_K_W':fit['b_kR_m']/(fit['alpha_m2_s']*C),'interpretation':'assumedC; not independently measuredR'} for C in [.8e6,1.5e6,3.5e6]];fit['ratio_residuals']=[]
  for z in o:
   pred=transfer(fit['alpha_m2_s'],fit['b_kR_m'],z['omega_rad_s']);measurement=complex(z['ratio_real'],z['ratio_imag']);fit['ratio_residuals'].append({'omega_rad_s':z['omega_rad_s'],'predicted_real':pred.real,'predicted_imag':pred.imag,'residual_absolute':abs(pred-measurement),'residual_in_declared_sensor_error_scales':abs(pred-measurement)/z['ratio_error_scale']})
  result[label]={'Robin_fit':fit,'b0_null':null,'objective_difference_null_minus_robin':null['weighted_objective']-fit['weighted_objective'],'not_a_likelihood_ratio_test':True}
 result['status']='BOUNDED_FEATURE_FIT_COMPLETE'
 (P/'RESULTS_FEATURE_ONLY.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
