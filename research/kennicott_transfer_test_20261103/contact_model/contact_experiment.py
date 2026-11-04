"""Additive TRAIN-informed basal-contact alternative; TEST targets stay sealed.

Uses frozen old temperature-only feature plan. Physical C fixed independently;
one R is fitted on2020 only. Old capacity-refusal results are not rewritten.
"""
import hashlib,json,math,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'kennicott_model'),str(ROOT/'algorithm')]
import experiment as base
from column_contact import simulate
PRIMARY='contact_layered_local'
CONTACT_METHODS=[PRIMARY,'contact_homogeneous_local','contact_fixed2020_layered','contact_weather_boundary_layered']
SPEC={'primary_candidate':PRIMARY,'development_timing':'Newphysicalalternative developedAFTER2020TRAIN360kg/m2 known;2023meltmagnitudes remainsealed','capacity_reference_J_m3K':1.5e6,'capacity_sensitivity_J_m3K':[.8e6,3.5e6],'capacity_fitted':False,'R_engineering_bracket_m2K_W':[0.,.5],'R_measured':False,'R_equivalent_dry_air_gap_m_at_k_point024':[0.,.012],'contact_physics':'Zeroheatcapacitypassiveresistance;debrisfaceT−icefaceT=qR; noenergy removedatcontact','calibration':'ONEglobalRperphysicalmethodandcapacityscene fittedonly2020TRAIN; nooutputrescaling','old_primary_refusal_preserved':True,'profile_parameters':'Useunchangedoldtemperature-onlyinference/geometry/forcing; no2023meltfit','reference_ice_k_W_mK':2.1,'reference_ice_C_J_m3K':917.*2100.,'latent_J_kg':333700.,'fixed_interface':True,'historical_physics_novelty_claim':False,'field_superiority_established':False}

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load_plan():return base.load(ROOT/'kennicott_model/FEATURE_ONLY_PLAN_FINAL.json')

def forward(plan,stake_id,method,R_m2K_W,capacity_J_m3K=1.5e6,*,proxy_error_multiplier=0.,initial_ice_T_C=0.,debris_cells=16,ice_cells=20,step_seconds=300.,probe_absolute_depths_m=()):
 if method not in CONTACT_METHODS or not math.isfinite(R_m2K_W)or not 0<=R_m2K_W<=.5:raise ValueError('frozencontactmethod/nonnegativeengineeringRrequired')
 site=plan['sites'][stake_id];b=base.boundary_for(plan,stake_id);fitrow=plan['profiles'][site['inverse_profile']];layered=fitrow['alpha_layered_m2_s'];middle=fitrow['depths_m'][1]
 if method=='contact_fixed2020_layered':layered=plan['profiles']['KDT_2020a']['alpha_layered_m2_s']
 if method=='contact_homogeneous_local':layered=[fitrow['alpha_homogeneous_m2_s']]*2
 upper=site['forcing_plane_m'];length=site['stake_thickness_m']-upper
 if length<=0:raise ValueError('stakeaboveforcingplane')
 alpha=base.layer_cell_alpha(upper,site['stake_thickness_m'],middle,layered,debris_cells)
 raw=b['proxy_boundary_C']if method=='contact_weather_boundary_layered'else b['boundary_C'];offset=proxy_error_multiplier*2*b['metadata']['blocked_RMSE_C']
 hourly=[v+(offset if method=='contact_weather_boundary_layered'or flag else 0.)for v,flag in zip(raw,b['is_proxy'])]
 if 3600%step_seconds:raise ValueError('stepmustdividehour')
 per=int(3600/step_seconds);surface=[a+(z-a)*(j+1)/per for a,z in zip(hourly[:-1],hourly[1:])for j in range(per)];init=[hourly[0]*(1-(j+.5)/debris_cells)for j in range(debris_cells)]
 r=simulate(surface_temperature_C=surface,interval_seconds=step_seconds,debris_depth_m=length,debris_alpha_cells_m2_s=alpha,debris_capacity_J_m3K=capacity_J_m3K,initial_debris_temperature_C=init,initial_ice_temperature_C=[initial_ice_T_C]*ice_cells,ice_depth_m=2.,record_intervals=True,basal_contact_resistance_m2K_W=R_m2K_W,probe_depths_m=[z-upper for z in probe_absolute_depths_m])
 trace=r.pop('trace');cumulative=[0.];total=0.;maxjump=0.;hourly_probes=[]
 for j,row in enumerate(trace):
  total+=row['melt_mmwe'];maxjump=max(maxjump,abs(row['interface_debris_side_T_C']-row['interface_T_C']))
  if(j+1)%per==0:
   cumulative.append(total);hourly_probes.append(row['debris_probe_temperature_C'])
 r.update({'method':method,'stake_id':stake_id,'capacity_J_m3K':capacity_J_m3K,'R_m2K_W':R_m2K_W,'layered_alpha_m2_s':layered,'layer_interface_depth_m':middle,'forcing_plane_m':upper,'stake_ice_depth_m':site['stake_thickness_m'],'no_below_ice_thermistor_forcing':True,'hourly_cumulative_melt_mmwe':cumulative,'clock_labels_UTC':b['times_utc'],'maximum_contact_temperature_jump_C':maxjump,'probe_absolute_depths_m':list(probe_absolute_depths_m),'hourly_probe_times_UTC':b['times_utc'][1:]if probe_absolute_depths_m else[],'hourly_probe_temperature_C':hourly_probes if probe_absolute_depths_m else[],'proxy_error_multiplier':proxy_error_multiplier,'initial_ice_T_C':initial_ice_T_C,'temperature_features_only':True,'posthoc_physical_ledger_rescaling':False})
 return r

def R_probe(plan,method,capacity_J_m3K=1.5e6):
 rows=[]
 for R in[0.,.05,.1,.2,.35,.5]:
  r=forward(plan,'P_AWS_1',method,R,capacity_J_m3K);rows.append({'R_m2K_W':R,'melt_mmwe':r['physical_melt_mmwe'],'relative_energy_residual':r['relative_energy_residual']})
 return {'method':method,'capacity_J_m3K':capacity_J_m3K,'sampled_R_response':rows,'sampled_strictly_decreasing':all(a['melt_mmwe']>b['melt_mmwe']for a,b in zip(rows[:-1],rows[1:])),'global_monotonicity_theorem_established':False,'TEST_targets_read':False}

def calibrate_R(plan,method,train_mass_loss_kg_m2,capacity_J_m3K=1.5e6,*,relative_tolerance=1e-5,max_iterations=24):
 if not math.isfinite(train_mass_loss_kg_m2)or train_mass_loss_kg_m2<=0:raise ValueError('positiveexplicitTRAINrequired')
 lo,hi=SPEC['R_engineering_bracket_m2K_W'];a=forward(plan,'P_AWS_1',method,lo,capacity_J_m3K)['physical_melt_mmwe'];b=forward(plan,'P_AWS_1',method,hi,capacity_J_m3K)['physical_melt_mmwe'];history=[{'R':lo,'melt':a},{'R':hi,'melt':b}]
 if not b<=train_mass_loss_kg_m2<=a:return {'status':'TRAIN_outside_fixed_R_bracket','method':method,'capacity_J_m3K':capacity_J_m3K,'TRAIN_mass_loss_kg_m2':train_mass_loss_kg_m2,'history':history,'output_rescaling':False,'TEST_targets_read':False}
 if train_mass_loss_kg_m2 in(a,b):return {'status':'calibrated','R_m2K_W':lo if train_mass_loss_kg_m2==a else hi,'method':method,'capacity_J_m3K':capacity_J_m3K,'physical_TRAIN_melt_mmwe':train_mass_loss_kg_m2,'history':history,'TEST_targets_read':False}
 for _ in range(max_iterations):
  R=(lo+hi)/2;r=forward(plan,'P_AWS_1',method,R,capacity_J_m3K);y=r['physical_melt_mmwe'];history.append({'R':R,'melt':y})
  if not b<=y<=a:return {'status':'nonmonotone_R_response_refusal','history':history,'TEST_targets_read':False}
  if abs(y-train_mass_loss_kg_m2)<=relative_tolerance*train_mass_loss_kg_m2:return {'status':'calibrated','R_m2K_W':R,'method':method,'capacity_J_m3K':capacity_J_m3K,'physical_TRAIN_melt_mmwe':y,'TRAIN_mass_loss_kg_m2':train_mass_loss_kg_m2,'history':history,'physical_energy_ledger_rescaled':False,'source_mass_density_reconverted':False,'TEST_targets_read':False}
  if y>train_mass_loss_kg_m2:lo,a=R,y
  else:hi,b=R,y
 return {'status':'calibration_iteration_limit_refusal','history':history,'TEST_targets_read':False}

if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=['probe','calibrate']);p.add_argument('--output',required=True);p.add_argument('--train-kg-m2',type=float);a=p.parse_args();plan=load_plan()
 if a.mode=='probe':r={'spec':SPEC,'probes':[R_probe(plan,m)for m in CONTACT_METHODS],'TEST_targets_read':False}
 else:
  if a.train_kg_m2 is None:p.error('TRAINvalue explicitlyrequired')
  r={'spec':SPEC,'calibrations':[calibrate_R(plan,m,a.train_kg_m2)for m in CONTACT_METHODS],'primary_capacity_sensitivity_calibrations':[calibrate_R(plan,PRIMARY,a.train_kg_m2,C)for C in[.8e6,3.5e6]],'TEST_targets_read':False}
 with Path(a.output).open('x')as f:json.dump(r,f,indent=2,allow_nan=False);f.write('\n')
 print(a.output)
