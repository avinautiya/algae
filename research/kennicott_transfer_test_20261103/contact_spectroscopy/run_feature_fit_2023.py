"""First/deepest native thermistors only. No outcome files or loader imports."""
import json,pathlib,hashlib,csv,datetime as dt
from .spectroscopy import harmonic_extract
from .profile_fit import fit_geometry,rank_geometry
P=pathlib.Path(__file__).parent

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 f=P/'FEATURE_ONLY_PLAN_2023.json';plan=json.loads(f.read_text());start=dt.datetime.fromisoformat(plan['fit_start_UTC']);end=dt.datetime.fromisoformat(plan['fit_end_UTC']);results=[]
 for s in plan['profiles']:
  raw=P.parent/'data'/'kennicott'/'features'/f"{s['profile']}.csv";assert sha(raw)==s['source_sha256'];rows=list(csv.DictReader(raw.read_text(encoding='utf-8-sig').splitlines()));r=[x for x in rows if start<=dt.datetime.fromisoformat(x['Datetime'])<=end];t=[(dt.datetime.fromisoformat(x['Datetime'])-start).total_seconds() for x in r];u=[float(x[s['upper_column']]) for x in r];y=[float(x[s['deep_column']]) for x in r];o={'profile':s['profile'],'stake':s['stake'],'rows':len(r),'geometry_scenes':[],'negative_temperatures_count':sum(x<0 for x in u+y)}
  try:fit_data=harmonic_extract(t,u,y)
  except ValueError as e:o.update({'status':'NO_FIT_SOURCE_RESOLUTION_GATE','reason':str(e)});results.append(o);continue
  o['harmonic_extraction']=fit_data
  for name,d in s['geometry_scenes_m'].items():
   L=d-s['upper_z_m'];x=s['deep_z_m']-s['upper_z_m'];g={'scene':name,'d_m':d,'L_m':L,'x_m':x}
   if x>L+1e-12:g.update({'status':'REFUSED_GEOMETRY','reason':'Deep probe below modeled ice boundary; do not relabel/move probe'});o['geometry_scenes'].append(g);continue
   for label,obs in [('primary_DC_and_harmonics',fit_data['observations']),('diagnostic_harmonics_only',fit_data['observations'][1:])]:
    z=fit_geometry(obs,L=L,x=x);n=fit_geometry(obs,L=L,x=x,null=True);z['local_rank']=rank_geometry(z['alpha_m2_s'],z['b_kR_m'],obs,L=L,x=x);g[label]={'Robin_fit':z,'b0_null':n,'not_a_likelihood_ratio_test':True}
   g['status']='FEATURE_FIT_COMPLETE';o['geometry_scenes'].append(g)
  o['status']='FEATURE_FIT_COMPLETE';results.append(o)
 out={'plan_sha256':sha(f),'base_OLS_sha256':sha(P/'spectroscopy.py'),'geometry_fit_sha256':sha(P/'profile_fit.py'),'runner_sha256':sha(P/'run_feature_fit_2023.py'),'2023_mass_values_used':False,'scope':'Fixed-geometry periodic thermal-response fit; no historical priority or field accuracy claim.','profiles':results};(P/'RESULTS_FEATURE_ONLY_2023.json').write_text(json.dumps(out,indent=2)+'\n')
 for p in results:
  print(p['profile'],p['status'],p.get('reason',''))
  for g in p['geometry_scenes']:
   if g['status']!='FEATURE_FIT_COMPLETE':print(' ',g['scene'],g['status']);continue
   q=g['primary_DC_and_harmonics'];h=g['diagnostic_harmonics_only'];z=q['Robin_fit'];print(' ',g['scene'], 'alpha',z['alpha_m2_s'],'b',z['b_kR_m'],'score',z['weighted_objective'],'null',q['b0_null']['weighted_objective'],'harm-only b',h['Robin_fit']['b_kR_m'],'rank',z['local_rank']['singular_values'])
if __name__=='__main__':main()
