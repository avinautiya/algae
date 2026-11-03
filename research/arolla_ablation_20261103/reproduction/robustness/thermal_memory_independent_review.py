"""Independent analytic review and declared-synthetic periodic demonstration.

Imports the additive model but changes none of its files. This is not a field
validation or an empirical heat-capacity estimate.
"""
import cmath,hashlib,json,math,sys
from pathlib import Path
MODEL=Path('/private/tmp/melt-validation-20261103/model')
sys.path.insert(0,str(MODEL))
from thermal_memory.solver import Column,SurfaceBoundary,simulate
P={'status':'declared_synthetic','source':'review synthetic analytic fixture','location':'not a glacier site','method':'prescribed independently of all melt observations','uses_evaluation_melt':False}
S={key:dict(P) for key in ['conductivity','heat_capacity','depth','initial_profile']}
S['heat_capacity']['quantity']='bulk_debris_volumetric_heat_capacity'
S['basal_ice_state']=dict(P,state='temperate_melting_interface',temperature_C=0.)
def solve(depth,n,step,hours,surface,initial,k=.96,C=1.5e6):
 return simulate(Column(k,C,depth,tuple(initial),S),SurfaceBoundary(tuple(surface),dict(P,quantity='debris_surface_temperature')),horizon_s=hours*3600,max_horizon_s=hours*3600,dt_s=step,synthetic=True)
def eigen_check():
 depth=.18;k=.96;C=1.5e6;n=32;step=30;hours=2
 dz=depth/n;alpha=k/C;decay=4*alpha*math.sin(math.pi/(2*n))**2/dz**2
 initial=[math.sin(math.pi*(j+.5)/n) for j in range(n)]
 r=solve(depth,n,step,hours,[0]*hours,initial,k,C)
 factor=(1+step*decay)**(-(hours*3600)//step)
 exact_discrete=[v*factor for v in initial]
 err=max(abs(a-b) for a,b in zip(r['final_temperature_C'],exact_discrete))
 return dict(max_temperature_error_against_exact_discrete_eigenmode_C=err,continuous_eigenvalue_s_inverse=alpha*(math.pi/depth)**2,discrete_eigenvalue_s_inverse=decay,energy_relative_residual=r['discrete_energy_relative_residual'])
def steady_single_cell():
 r=solve(.18,1,60,1,[10],[5])
 expected=.96*10/.18
 return dict(expected_flux_W_m2=expected,max_flux_error_W_m2=max(abs(h['q_basal_W_m2']-expected) for h in r['history']),max_step_ledger_J_m2=r['max_absolute_step_residual_J_m2'])
def harmonic(data,key):
 # Uniform full-period samples. A+iB is the cosine/sine coefficient pair.
 rows=[r for r in data if r['time_s']>6*86400]
 omega=2*math.pi/86400
 a=2*sum(r[key]*math.cos(omega*r['time_s']) for r in rows)/len(rows)
 b=2*sum(r[key]*math.sin(omega*r['time_s']) for r in rows)/len(rows)
 return complex(a,b)
def demo(depth):
 n=40;step=60;hours=8*24
 surface=[10+8*math.sin(2*math.pi*(i+.5)/24) for i in range(hours)]
 initial=[10*(1-(i+.5)/n) for i in range(n)]
 r=solve(depth,n,step,hours,surface,initial)
 out=harmonic(r['history'],'q_basal_W_m2');inp=harmonic(r['history'],'surface_temperature_C')
 # A+iB convention conjugates the e^(iwt) physical transfer function.
 numerical=out/inp;omega=2*math.pi/86400;k=.96;C=1.5e6
 lam=cmath.sqrt(1j*omega*C/k);analytic=k*lam/cmath.sinh(lam*depth)
 lag=cmath.phase(numerical)/omega/3600
 return dict(depth_m=depth,conductivity_W_mK=k,volumetric_heat_capacity_J_m3K=C,
  status='synthetic dimensional demonstration; no field fit or novelty claim',
  numerical_transfer_W_m2_per_K=[numerical.real,numerical.imag],
  continuous_analytic_transfer_W_m2_per_K=[analytic.real,analytic.imag],
  numerical_basal_harmonic_amplitude_W_m2=abs(out),quasistatic_amplitude_W_m2=abs(inp)*k/depth,
  numerical_lag_modulo_24h=lag%24,analytic_lag_modulo_24h=(-cmath.phase(analytic)/omega/3600)%24,
  phase_lag_is_modulo_one_forcing_period=True,
  relative_transfer_error=abs(numerical-analytic.conjugate())/abs(analytic),
  total_incoming_heat_J_m2=r['surface_heat_J_m2'],total_basal_heat_J_m2=r['basal_heat_J_m2'],
  final_minus_initial_storage_J_m2=r['final_sensible_storage_J_m2']-r['initial_sensible_storage_J_m2'],
  energy_relative_residual=r['discrete_energy_relative_residual'])
def continuum_eigen_convergence():
 depth=.18;alpha=.96/1.5e6;hours=1;step=1
 rows=[]
 for n in [8,16,32]:
  initial=[math.sin(math.pi*(j+.5)/n) for j in range(n)]
  r=solve(depth,n,step,hours,[0]*hours,initial)
  exact=[v*math.exp(-alpha*(math.pi/depth)**2*hours*3600) for v in initial]
  rows.append(dict(cells=n,dt_s=step,max_temperature_error_C=max(abs(a-b) for a,b in zip(r['final_temperature_C'],exact))))
 return rows
if __name__=='__main__':
 p=Path(__file__).parent/'THERMAL_MEMORY_INDEPENDENT_RESULTS_V3.json'
 if p.exists():raise ValueError('refuse overwrite')
 result=dict(model_solver_sha256=hashlib.sha256((MODEL/'thermal_memory/solver.py').read_bytes()).hexdigest(),single_cell=steady_single_cell(),eigenmode=eigen_check(),continuum_eigen_convergence=continuum_eigen_convergence(),demos=[demo(.18),demo(.6)],field_validation=False)
 p.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
