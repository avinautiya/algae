"""Post-study independent early-time numerical audit; leaves frozen files intact."""
import sys,json,time,hashlib
from pathlib import Path
ROOT=Path('/private/tmp/cliff-upgrade-20261101/package');sys.path.insert(0,str(ROOT))
from thermal.resolved_slab import SlabParameters,SlabState,advance
from thermal.robin_reference import RobinSlab
OUT=Path(__file__).resolve().parent
cases=json.loads((ROOT/'mechanism/wet_boundary_cases.json').read_text())
rows=[]
for case in cases:
 p=SlabParameters(case['thickness_m'],case['conductivity_W_m_K'],case['volumetric_heat_capacity_J_m3_K'],case['wet_area_m2'])
 K=case['ice_conductance_W_K'];C=p.diffusion_time_s*K;t=.1*p.diffusion_time_s
 ref=RobinSlab.from_contact(C,K,case['water_capacity_W_K'],case['film_conductance_W_K'],case['inlet_superheat_K'],case['absorbed_source_W'],terms=512).response(t)['ice_heat_J']
 for n in [64,128]:
  for steps in [4000,16000]:
   state=SlabState((0.,)*n);ice=0.;res=0.;start=time.perf_counter()
   for j in range(steps//10):
    r=advance(state,p,p.diffusion_time_s/steps,case['inlet_superheat_K'],case['water_capacity_W_K'],debris_water_conductance_W_K=case['film_conductance_W_K'],atmospheric_power_W=case['absorbed_source_W'])
    state=r['state'];ice+=r['covered_ice_heat_J'];res=max(res,abs(r['global_energy_residual_J']))
   rows.append(dict(id=case['id'],cells=n,steps_per_tau=steps,analytic_J=ref,FV_J=ice,relative_error=(ice-ref)/ref,max_discrete_residual_J=res,runtime_s=time.perf_counter()-start))
 (OUT/'RESULTS.json').write_text(json.dumps(dict(scope='post-study numerical audit; not frozen production results',rows=rows),indent=2))
 print(case['id'],flush=True)
