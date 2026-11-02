"""Small physical-time v5 smoke/observation check, not a glacier campaign."""
import argparse,hashlib,json,sys,time
from pathlib import Path
import numpy as np
if __package__:
 from .matched_roi import capture,observe,save_bundle,load_bundle,compare
 from ..phase.swe2d_v5 import Sim,SOLVER_VERSION
else:
 from matched_roi import capture,observe,save_bundle,load_bundle,compare
 sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'phase'))
 from swe2d_v5 import Sim,SOLVER_VERSION

CONFIGS=[dict(id='coarse',dx=.1,step_s=.05),dict(id='fine',dx=.05,step_s=.05),dict(id='fine_dt_half',dx=.05,step_s=.025)]
END_S=2.
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
 if (a.out/'SPEC.json').exists():raise ValueError('new output directory required; no overwrite')
 source=Path(__file__).resolve().parents[1]/'phase/swe2d_v5.py';sha=hashlib.sha256(source.read_bytes()).hexdigest()
 base_configuration=dict(Lpool=.2,Lr=.4,Ltail=.2,Wl=.2,Q=.02,Tw=.5,ramp_ratio=0.)
 source_root=source.parent;source_hashes={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in [source,source_root/'swe2d_v4_pinned.py',Path(__file__),Path(__file__).with_name('matched_roi.py')]}
 (a.out/'SPEC.json').write_text(json.dumps(dict(configurations=CONFIGS,base_configuration=base_configuration,physical_duration_s=END_S,source_sha256=source_hashes,purpose='short numerical verification only; no environmental prediction or novelty claim'),indent=2)+'\n')
 records=[]
 for cfg in CONFIGS:
  start=time.perf_counter();sim=Sim(dx=cfg['dx'],**base_configuration)
  ref=capture(sim,roi_id='fixed_toe_reach',bounds_m=(.2,.6,0.,.2),frame='idealised_cartesian_grid_origin0_m',provenance=dict(definition='fixed rectangular toe-reach cells',source='explicit verification configuration, not measured terrain'))
  sim.regions={'fixed_toe_reach':ref.mask.copy()};sim.region_vol={};history=[]
  for _ in range(round(END_S/cfg['step_s'])):history.append(sim.step_physical(cfg['step_s']))
  obs=observe(sim,ref,physical_time_s=END_S,heat_implied_volume_m3=float(sim.region_vol['fixed_toe_reach']),absolute_ledger_tolerance_m3=1e-12)
  led=sim.combined_ledgers()
  if not led['physical_split_only'] or abs(led['combined_ice_water_mass_residual_kg'])>1e-8 or abs(led['physical_split_thermal_latent_residual_J'])>1e-6:raise RuntimeError('short physical verification ledger failed')
  if abs(sim.physical_time_s-END_S)>1e-12:raise RuntimeError('physical time mismatch')
  meta=dict(solver_version=SOLVER_VERSION,solver_source_sha256=sha,source_sha256=source_hashes,configuration=dict(**base_configuration,**cfg,physical_duration_s=END_S),evidence_class='conditional_numerical_run')
  path=a.out/(cfg['id']+'.json');save_bundle(path,sim,ref,obs,run_metadata=meta);load_bundle(path)
  summary=dict(id=cfg['id'],observation=obs,ledgers=led,wall_s=time.perf_counter()-start);records.append(summary)
  print(cfg['id'],obs['geometric_melt_m3'],led['combined_ice_water_mass_residual_kg'],flush=True)
 (a.out/'SUMMARY.json').write_text(json.dumps(dict(records=records,coarse_vs_fine=compare(records[0]['observation'],records[1]['observation'],relative_reporting_min_volume_m3=1e-10),fine_step_vs_halfstep=compare(records[1]['observation'],records[2]['observation'],relative_reporting_min_volume_m3=1e-10),interpretation='Three short transient verification configurations; no field accuracy, event-timing convergence or scientific discovery established.'),indent=2,allow_nan=False)+'\n')
if __name__=='__main__':main()
