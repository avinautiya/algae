"""Pre-TRAIN numerical gate using source forcing, never melt outcomes."""
import json,time
from pathlib import Path
from experiment import prepare_plan,forward,PRIMARY

def run():
 plan=prepare_plan();results=[]
 for stake in ['P_AWS_1','P11','P12','P13']:
  for method in [PRIMARY,'physical_homogeneous_local']:
   rows=[]
   for n,h in[(16,300.),(16,150.),(32,150.)]:
    start=time.monotonic();r=forward(plan,stake,method,1.5e6,debris_cells=n,ice_cells=20,step_seconds=h)
    rows.append({'debris_cells':n,'ice_cells':20,'dt_seconds':h,'melt_mmwe':r['physical_melt_mmwe'],'relative_energy_residual':r['relative_energy_residual'],'wall_seconds':time.monotonic()-start})
   timechange=(rows[1]['melt_mmwe']-rows[0]['melt_mmwe'])/rows[1]['melt_mmwe'];spacechange=(rows[2]['melt_mmwe']-rows[1]['melt_mmwe'])/rows[2]['melt_mmwe'];overall=(rows[2]['melt_mmwe']-rows[0]['melt_mmwe'])/rows[2]['melt_mmwe']
   results.append({'stake':stake,'method':method,'runs':rows,'dt_refinement_relative_change':timechange,'debris_mesh_refinement_relative_change':spacechange,'combined_relative_change':overall,'one_percent_seasonal_gate_passed':max(abs(timechange),abs(spacechange),abs(overall))<.01})
 return {'source_forced_pretrain_verification':True,'capacity_reference_J_m3K':1.5e6,'layer_mapping':'exactcellharmonicintersection','all_comparisons_pass':all(x['one_percent_seasonal_gate_passed']for x in results),'results':results,'native_melt_targets_read':False,'ice_mesh_refined_here':False,'engineering_gate_not_confidence_interval':True}
if __name__=='__main__':
 r=run()
 with Path(__file__).with_name('SEASONAL_REFINEMENT_FINAL.json').open('x')as f:json.dump(r,f,indent=2,allow_nan=False);f.write('\n')
 print(json.dumps(r,indent=2))
