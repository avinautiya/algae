"""Temperature-only Robin/extra-clearance discrimination diagnostic."""
import cmath,math,json,pathlib
from .spectroscopy import transfer
P=pathlib.Path(__file__).parent

def shifted_boundary_transfer(alpha,b,omega,L,x):
 if omega==0:return complex((L+b-x)/(L+b),0)
 g=cmath.sqrt(1j*omega/alpha);shift=cmath.atanh(b*g)/g
 return cmath.sinh(g*(L+shift-x))/cmath.sinh(g*(L+shift))

def _row(profile,fit,obs,L,x):
 a,b=fit['alpha_m2_s'],fit['b_kR_m'];rows=[]
 for o in obs:
  w=o['omega_rad_s'];r=transfer(a,b,w,L,x);g=transfer(a,0,w,L+b,x);v=shifted_boundary_transfer(a,b,w,L,x);rows.append({'omega_rad_s':w,'absolute_Robin_minus_geometric_extension':abs(r-g),'difference_in_source_error_scales':abs(r-g)/o['ratio_error_scale'],'exact_complex_shift_identity_error':abs(r-v)})
 return {'profile':profile,'alpha_m2_s':a,'effective_contact_length_m':b,'pure_geometry_alternative_d_extension_m':b,'rows':rows}

def main():
 a=json.loads((P/'RESULTS_FEATURE_ONLY.json').read_text());rows=[_row('KDT_2020b',a['primary_DC_and_harmonics']['Robin_fit'],a['harmonic_extraction']['observations'],.17,.1)]
 b=json.loads((P/'RESULTS_FEATURE_ONLY_2023.json').read_text())
 for p in b['profiles']:
  if p['status']!='FEATURE_FIT_COMPLETE':continue
  g=next(x for x in p['geometry_scenes'] if x['scene']=='nominal_profile');rows.append(_row(p['profile'],g['primary_DC_and_harmonics']['Robin_fit'],p['harmonic_extraction']['observations'],g['L_m'],g['x_m']))
 result={'mass_values_used':False,'scope':'Change in temperature response at same probe position; not a continuous geometry bound or proof equivalent physical melt.','profiles':rows};(P/'GEOMETRY_EQUIVALENCE_DIAGNOSTIC.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
