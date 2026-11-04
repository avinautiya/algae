"""Additional finite-grid compatibility diagnostic; not a confidence region."""
import json,pathlib,math
from .spectroscopy import transfer
P=pathlib.Path(__file__).parent
r=json.loads((P/'RESULTS_FEATURE_ONLY.json').read_text());o=r['harmonic_extraction']['observations'];out={}
for label,obs in [('DC_and_harmonics',o),('harmonics_only',o[1:])]:
 accepted=[]
 for i in range(61):
  alpha=math.exp(math.log(1e-8)+i/60*math.log(500))
  for j in range(61):
   b=.3*j/60;maxres=max(abs(transfer(alpha,b,x['omega_rad_s'])-complex(x['ratio_real'],x['ratio_imag']))/x['ratio_error_scale'] for x in obs)
   if maxres<=1:accepted.append((alpha,b))
 out[label]={'accepted_finite_grid_points':len(accepted),'alpha_grid_range':[min(a for a,b in accepted),max(a for a,b in accepted)] if accepted else None,'b_grid_range_m':[min(b for a,b in accepted),max(b for a,b in accepted)] if accepted else None,'b0_accepted_on_finite_grid':any(b==0 for a,b in accepted),'criterion':'Each complex-ratio residual <= its deterministic propagated source-error scale','not_statistical_confidence_or_continuous_enclosure':True,'finite_window_model_discrepancy_not_included':True}
(P/'ADDITIONAL_GRID_DIAGNOSTIC.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
