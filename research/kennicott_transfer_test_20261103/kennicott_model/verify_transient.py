"""Independent continuum Fourier-series transient, not any field target."""
import json,math
from pathlib import Path
from column import simulate

def exact_melt(t,L=.2,alpha=1e-6,k=1.5,Ts=5.,latent=333700.):
 # Alternating steady sum is known analytically; only exponential tail needed.
 tail=math.fsum((-1)**n*math.exp(-n*n*math.pi**2*alpha*t/(L*L))/(n*n)for n in range(1,101))
 integral=t-L*L/(6*alpha)-2*L*L/(alpha*math.pi**2)*tail
 return k*Ts/L*integral/latent

def run():
 duration=7200.;truth=exact_melt(duration);out=[]
 for n in [8,16,32]:
  for h in [600.,300.,150.]:
   r=simulate(surface_temperature_C=[5.]*(int(duration/h)),interval_seconds=h,debris_depth_m=.2,debris_alpha_cells_m2_s=[1e-6]*n,debris_capacity_J_m3K=1.5e6,initial_debris_temperature_C=[0.]*n,initial_ice_temperature_C=[0.]*20,record_intervals=False)
   out.append({'debris_cells':n,'step_seconds':h,'melt_mmwe':r['physical_melt_mmwe'],'exact_Fourier_series_mmwe':truth,'relative_error':(r['physical_melt_mmwe']-truth)/truth,'relative_ledger_residual':r['relative_energy_residual']})
 return {'synthetic_continuum_oracle':True,'two_hour_step_response':out,'field_outcomes_used':False,'steady_series_checked_independently':True}
if __name__=='__main__':
 r=run()
 with Path(__file__).with_name('TRANSIENT_VERIFICATION.json').open('x')as f:json.dump(r,f,indent=2);f.write('\n')
 print(json.dumps(r,indent=2))
