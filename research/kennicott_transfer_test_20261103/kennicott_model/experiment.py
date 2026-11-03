"""Source-feature-only plan and bounded physical-capacity calibration.

No outcome reader exists here. TRAIN value may be supplied only after external
method publication; 2023 values must remain sealed until predictions publish.
"""
import csv,datetime as dt,hashlib,json,math,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'algorithm')]
from three_sensor import fit,integrated_rows,UnsupportedModel
from column import simulate
DATA=ROOT/'data/kennicott'
PRIMARY='physical_layered_local'
THERMAL=[PRIMARY,'physical_homogeneous_local','physical_fixed2020_layered','physical_weather_boundary_layered']
INDICES=['air_pdd','buried_temp_exposure','conventional_gradient','depth_corrected_air_pdd','constant_rate']
EMPIRICAL_THERMAL=['homogeneous_empirical_scale']
MODEL_IDS=THERMAL+EMPIRICAL_THERMAL+INDICES
SPEC={'primary':PRIMARY,'model_ids':MODEL_IDS,'temperature_parameter_fit':'deepest3 REAL measured sensors, first14days after5dayinstallationreset','training_material_profile':'KDT_2020a transferredtoKDT_2020b by identicalsummarycoordinates/cotime/grainclass; notsamepitproof','synthetic_zero_2020b_depth22cm_not_used_in_inverse':True,'balance_window_seconds':21600.,'minimum_singular_ratio':.001,'surface_forcing':'firstREAL measured sensorplane, observedafterreset plusAIR2proxybefore reset; hourlylinearreconstruction','forward_step_seconds':300.,'debris_cells':16,'within_cell_layer_rule':'exactharmonicintersection preservesseriesresistance','ice_cells':20,'ice_depth_m':2.,'ice_k_W_mK':2.1,'ice_C_J_m3K':917.*2100.,'ice_initial_T_C':0.,'latent_heat_J_kg':333700.,'ice_model_density_kg_m3':917.,'source_mass_density_kg_m3':900.,'fixed_interface_no_ice_advection':True,'capacity_engineering_bracket_J_m3K':[.8e6,3.5e6],'capacity_probe_J_m3K':[.8e6,1.5e6,3.5e6],'capacity_bounds_are_measured_coverage':False,'global_capacity_calibration':'ONE2020TRAINmassloss only, no outputrescaling ofphysicalmodels','physical_baseline_calibration':'homogeneous/weather variants eachreceiveoneglobal2020capacity parameter inSAMEbracket; fixed2020sharesprimary2020fit','homogeneous_empirical_comparator':'Samehomogeneousalpha/boundary/geometry withfixedC1.5MJrawphysicalresponse, thenone2020TRAINoutputgain; rescaledoutputNOTphysicalenergyledger','mass_target_density_reconversion':False,'temperature_conditioned_hindcast_not_advance_forecast':True}

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):return json.loads(Path(p).read_text())
def profile_fit(name):
 path=DATA/'features'/f'{name}.csv';manifest=load(DATA/'FEATURE_EXPORT_MANIFEST.json');record=next(x for x in manifest['files']if x['local_path']=='features/'+name+'.csv')
 if sha(path)!=record['sha256']:raise ValueError('sourcefeaturehashmismatch')
 with path.open()as f:rows=list(csv.DictReader(f))
 cols=[c for c in rows[0]if c.startswith('Depth:')][-3:];z=[float(c.split(':')[1].split('cm')[0].strip())/100 for c in cols]
 first=dt.datetime.fromisoformat(rows[0]['Datetime']);start=first+dt.timedelta(days=5);end=start+dt.timedelta(days=14);r=[row for row in rows if start<=dt.datetime.fromisoformat(row['Datetime'])<=end]
 times=[(dt.datetime.fromisoformat(row['Datetime'])-start).total_seconds()for row in r];T=[[float(row[c])for c in cols]for row in r]
 if len(times)!=14*24*12+1:raise UnsupportedModel('complete5minute14dayfeaturefitrequired')
 if any(all(row[j]==0 for row in T)for j in range(3)):raise UnsupportedModel('constantzeroendpointisnotindependentthermistor')
 f=fit(z,times,T,balance_window_seconds=21600,minimum_singular_ratio=.001,maximum_gap_seconds=300,sensor_absolute_bound_C=.2)
 if not f['supported_point_estimate']:raise UnsupportedModel('temperatureonlylayeredidentifiabilitygatefailed')
 A,y,_=integrated_rows(z,times,T,21600,300);x=[sum(row)for row in A];den=math.fsum(v*v for v in x)
 if den<=0:raise UnsupportedModel('homogeneousinverseunidentified')
 hom=math.fsum(a*b for a,b in zip(x,y))/den
 if hom<=0:raise UnsupportedModel('homogeneousinversegivesnonpositivealpha')
 return {'profile':name,'source_sha256':sha(path),'depths_m':z,'columns':cols,'fit_start_UTC':start.isoformat()+'Z','fit_end_UTC':end.isoformat()+'Z','alpha_layered_m2_s':f['alpha_m2_s'],'alpha_homogeneous_m2_s':hom,'fit_diagnostics':f}

def prepare_plan():
 binding=load(DATA/'SOURCE_BINDING.json');profiles={name:profile_fit(name)for name in ['KDT_2020a','KDT_2023a','KDT_2023b','KDT_2023c']};sites={}
 for pair in binding['pairs']:
  name=pair['profile_id'];bpath=ROOT/'kennicott_forcing/nominal_boundaries'/f'{name}__Debris_T2.json';b=load(bpath)
  sites[pair['stake_id']]={'source_pair':pair,'boundary_path':str(bpath.relative_to(ROOT)),'boundary_sha256':sha(bpath),'inverse_profile':'KDT_2020a'if name=='KDT_2020b'else name,'stake_thickness_m':pair['stake_debris_thickness_m'],'forcing_plane_m':b['metadata']['upper_boundary_plane_m'],'boundary_metadata':b['metadata']}
 return {'spec':SPEC,'profiles':profiles,'sites':sites,'source_binding_sha256':sha(DATA/'SOURCE_BINDING.json'),'melt_targets_read':False,'historical_novelty_claim':False}

def boundary_for(plan,stake_id):
 r=plan['sites'][stake_id];p=ROOT/r['boundary_path']
 if sha(p)!=r['boundary_sha256']:raise ValueError('frozenboundarychanged')
 return load(p)

def layer_cell_alpha(upper,lower,middle,alpha_pair,cells):
 if cells<1 or lower<=upper or min(alpha_pair)<=0:raise ValueError('positivegeometry/material')
 dz=(lower-upper)/cells;result=[]
 for j in range(cells):
  a=upper+j*dz;b=a+dz;left=max(min(b,middle)-a,0.);left=min(left,dz);right=dz-left
  result.append(dz/(left/alpha_pair[0]+right/alpha_pair[1]))
 return result

def forward(plan,stake_id,method,capacity_J_m3K,*,proxy_error_multiplier=0.,initial_ice_T_C=0.,debris_cells=16,ice_cells=20,step_seconds=300.):
 if method not in THERMAL:raise ValueError('thermalmethodrequired')
 site=plan['sites'][stake_id];b=boundary_for(plan,stake_id);fitrow=plan['profiles'][site['inverse_profile']];layered=fitrow['alpha_layered_m2_s'];middle=fitrow['depths_m'][1]
 if method=='physical_fixed2020_layered':layered=plan['profiles']['KDT_2020a']['alpha_layered_m2_s']
 if method=='physical_homogeneous_local':layered=[fitrow['alpha_homogeneous_m2_s']]*2
 upper=site['forcing_plane_m'];length=site['stake_thickness_m']-upper
 if length<=0:raise UnsupportedModel('stakebelowforcingplane')
 alpha=layer_cell_alpha(upper,site['stake_thickness_m'],middle,layered,debris_cells)
 raw=b['proxy_boundary_C']if method=='physical_weather_boundary_layered'else b['boundary_C'];offset=proxy_error_multiplier*2*b['metadata']['blocked_RMSE_C']
 # Only predeclared reconstructed epochs shift in the primary scene.
 hourly=[v+(offset if method=='physical_weather_boundary_layered'or flag else 0.)for v,flag in zip(raw,b['is_proxy'])]
 if 3600%step_seconds:raise ValueError('stepmustdividehour')
 per=int(3600/step_seconds);surface=[a+(z-a)*(j+1)/per for a,z in zip(hourly[:-1],hourly[1:])for j in range(per)]
 init=[hourly[0]*(1-(j+.5)/debris_cells)for j in range(debris_cells)]
 r=simulate(surface_temperature_C=surface,interval_seconds=step_seconds,debris_depth_m=length,debris_alpha_cells_m2_s=alpha,debris_capacity_J_m3K=capacity_J_m3K,initial_debris_temperature_C=init,initial_ice_temperature_C=[initial_ice_T_C]*ice_cells,ice_depth_m=2.,record_intervals=True)
 cumulative=[0.];total=0.
 for j,row in enumerate(r.pop('trace')):
  total+=row['melt_mmwe']
  if(j+1)%per==0:cumulative.append(total)
 r.update({'method':method,'stake_id':stake_id,'capacity_J_m3K':capacity_J_m3K,'layered_alpha_m2_s':layered,'layer_interface_depth_m':middle,'forcing_plane_m':upper,'stake_ice_depth_m':site['stake_thickness_m'],'no_below_ice_thermistor_forcing':True,'hourly_cumulative_melt_mmwe':cumulative,'clock_labels_UTC':b['times_utc'],'proxy_error_multiplier':proxy_error_multiplier,'initial_ice_T_C':initial_ice_T_C,'temperature_features_only':True})
 return r

def positive_linear(a,b,duration=3600.):
 if min(a,b)>=0:return duration*(a+b)/2
 if max(a,b)<=0:return 0.
 return duration*max(a,b)**2/(2*abs(b-a))

def index_doses(plan,stake):
 b=boundary_for(plan,stake);upper=b['boundary_C'];air=b['air_C'];length=plan['sites'][stake]['stake_thickness_m']-plan['sites'][stake]['forcing_plane_m'];buried=math.fsum(positive_linear(a,z)for a,z in zip(upper[:-1],upper[1:]));return {'air_pdd':math.fsum(positive_linear(a,z)for a,z in zip(air[:-1],air[1:])),'buried_temp_exposure':buried,'conventional_gradient':buried/length,'depth_corrected_air_pdd':math.fsum(positive_linear(a,z)for a,z in zip(air[:-1],air[1:]))/length,'constant_rate':float(len(upper)-1)*3600}

def empirical_homogeneous_reference(plan,stake_id):
 r=forward(plan,stake_id,'physical_homogeneous_local',1.5e6)
 return {'method':'homogeneous_empirical_scale','reference_dose':r['physical_melt_mmwe'],'reference_capacity_J_m3K':1.5e6,'reference_physical_ledger':{k:r[k]for k in ['initial_storage_J_m2','final_storage_J_m2','surface_heat_J_m2','deep_heat_out_J_m2','latent_heat_J_m2','energy_residual_J_m2','relative_energy_residual']},'one_global_TRAIN_output_gain_required':True,'rescaled_output_energy_conservation_claim':False}

def capacity_probe(plan,method):
 values=[]
 for C in SPEC['capacity_probe_J_m3K']:
  r=forward(plan,'P_AWS_1',method,C);values.append({'capacity_J_m3K':C,'melt_mmwe':r['physical_melt_mmwe'],'relative_energy_residual':r['relative_energy_residual']})
 return {'method':method,'probe':values,'sampled_response_strictly_increasing':all(a['melt_mmwe']<b['melt_mmwe']for a,b in zip(values[:-1],values[1:])),'all_capacity_monotonicity_theorem_established':False,'train_target_read':False}

def calibrate_capacity(plan,method,train_mass_loss_kg_m2,*,relative_tolerance=1e-5,max_iterations=24):
 """Requires external method freeze+TRAINrelease; NOT a target-reader."""
 if not math.isfinite(train_mass_loss_kg_m2)or train_mass_loss_kg_m2<=0:raise UnsupportedModel('positiveexplicitTRAINmasslossrequired')
 lo,hi=SPEC['capacity_engineering_bracket_J_m3K'];a=forward(plan,'P_AWS_1',method,lo)['physical_melt_mmwe'];b=forward(plan,'P_AWS_1',method,hi)['physical_melt_mmwe'];history=[{'C':lo,'melt':a},{'C':hi,'melt':b}]
 if not a<=train_mass_loss_kg_m2<=b:raise UnsupportedModel('TRAINtargetoutsideengineeredphysicalcapacitybracket; nooutputrescaling')
 if train_mass_loss_kg_m2 in (a,b):return {'capacity_J_m3K':lo if train_mass_loss_kg_m2==a else hi,'physical_TRAIN_melt_mmwe':train_mass_loss_kg_m2,'TRAIN_mass_loss_kg_m2':train_mass_loss_kg_m2,'history':history,'physical_energy_ledger_rescaled':False,'source_mass_density_reconverted':False}
 for _ in range(max_iterations):
  mid=(lo+hi)/2;r=forward(plan,'P_AWS_1',method,mid);y=r['physical_melt_mmwe'];history.append({'C':mid,'melt':y})
  if not a<=y<=b:raise UnsupportedModel('nonmonotonecapacityresponserootnotcertified')
  if abs(y-train_mass_loss_kg_m2)<=relative_tolerance*train_mass_loss_kg_m2:return {'capacity_J_m3K':mid,'physical_TRAIN_melt_mmwe':y,'TRAIN_mass_loss_kg_m2':train_mass_loss_kg_m2,'history':history,'physical_energy_ledger_rescaled':False,'source_mass_density_reconverted':False}
  if y<train_mass_loss_kg_m2:lo,a=mid,y
  else:hi,b=mid,y
 raise UnsupportedModel('physicalcapacitycalibrationdidnotmeetfrozenaccuracy')

if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=['prepare','probe','reference']);p.add_argument('--output',required=True);a=p.parse_args();plan=prepare_plan()
 if a.mode=='prepare':result=plan
 elif a.mode=='probe':result={'physical_capacity_probes':[capacity_probe(plan,m)for m in [PRIMARY,'physical_homogeneous_local','physical_weather_boundary_layered']],'no_melt_targets_read':True}
 else:
  result={'spec':SPEC,'no_melt_targets_read':True,'reference_capacity_J_m3K':1.5e6,'sites':{}}
  for stake in plan['sites']:
   thermal={m:forward(plan,stake,m,1.5e6)for m in THERMAL}
   result['sites'][stake]={'physical_reference':thermal,'homogeneous_empirical_reference_dose':thermal['physical_homogeneous_local']['physical_melt_mmwe'],'index_doses':index_doses(plan,stake)}
 with Path(a.output).open('x')as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
 print(a.output)
