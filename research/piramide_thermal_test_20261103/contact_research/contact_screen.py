"""Exploratory TRAIN-only contact plausibility screen, not a new fitted model.
Uses already released Tower1 TRAIN calibration and independent temperature data.
No Tower2 camera or future melt outcome is read. Steady assumptions are explicit.
"""
import argparse,cmath,csv,datetime as dt,hashlib,json,math
from pathlib import Path

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def transfer(omega,k,C,depth,R):
 if not all(type(v) in (float,int) and math.isfinite(v) for v in [omega,k,C,depth,R]) or min(k,C,depth)<=0 or R<0 or omega<0:raise ValueError('finite physical parameters required')
 if omega==0:return complex(1/(depth/k+R))
 lam=cmath.sqrt(1j*omega*C/k)
 return k*lam/(cmath.sinh(lam*depth)+R*k*lam*cmath.cosh(lam*depth))
def same_sensor_resistance(R,clearance,k):return R+clearance/k

def screen(process,normalized,predictions):
 n=json.loads(Path(normalized).read_text());p=json.loads(Path(predictions).read_text());name=p['primary_candidate']
 row=next(r for r in n['results'] if f"{'forward' if r['method']=='training_temperature_forward_fit' else 'derivative_OLS'}__{r['geometry']['id']}__{r['hold']}__{r['capacity_id']}"==name)
 calibration=p['calibrations'][name];gain=calibration['postprocessing_gain_relative_to_physical_component'];k=row['debris_conductivity_W_mK'];C=row['physicalcapacity_J_m3K'];geom=row['geometry'];depth=geom['ice_depth_m']-geom['sensor_depths_m'][0]
 a,b=map(dt.datetime.fromisoformat,[n['start_local'],n['split_local']]);totals={f'DT{i}':0. for i in range(1,6)};covered=0.;minimum=float('inf')
 with Path(process).open() as stream:rows=list(csv.DictReader(stream))
 for first,last in zip(rows,rows[1:]):
  t0=dt.datetime.fromisoformat(first['SOURCE_LOCAL_LABEL']);t1=dt.datetime.fromisoformat(last['SOURCE_LOCAL_LABEL']);seconds=max(0.,(min(b,t1)-max(a,t0)).total_seconds())
  if not seconds:continue
  if not 0<(t1-t0).total_seconds()<=1800:raise ValueError('unsupported TRAIN temperature gap')
  covered+=seconds
  for key in totals:
   value=float(first[key])
   if not math.isfinite(value):raise ValueError('nonfinite TRAIN temperature')
   totals[key]+=value*seconds
  minimum=min(minimum,float(first['DT1']))
 if covered!=(b-a).total_seconds():raise ValueError('complete TRAIN support required')
 means={key:value/covered for key,value in totals.items()};L=333700.;q=p['calibration']['observed_ice_ablation_mmwe']*L/covered
 Rd=depth/k;R_equiv=Rd*(1/gain-1);R_sensor_zero_clearance=means['DT5']/q;R_sensor_upper=(means['DT5']+.4)/q
 clearance=geom['ice_depth_m']-geom['sensor_depths_m'][-1]
 spectra=[]
 for hours in [6,12,24,48,168]:
  w=2*math.pi/(hours*3600);h0=transfer(w,k,C,depth,0);hr=transfer(w,k,C,depth,R_equiv)
  spectra.append(dict(period_hours=hours,R0_transfer_magnitude=abs(h0),contact_transfer_magnitude=abs(hr),contact_to_R0_magnitude_ratio=abs(hr/h0),contact_minus_R0_phase_rad=cmath.phase(hr/h0)))
 return dict(scope='Exploratory steady-temperate TRAIN-only plausibility/identifiability calculation; not a transient field falsification or new model',source_sha256=dict(process=sha(process),normalized=sha(normalized),released_train_predictions=sha(predictions)),train_start=a.isoformat(),train_end=b.isoformat(),temperature_means_C=means,minimum_TRAIN_upper_temperature_C=minimum,observed_TRAIN_mean_latent_flux_W_m2=q,source_model=dict(k_W_mK=k,C_J_m3K=C,computational_debris_depth_m=depth,empirical_gain_relative_raw=gain),steady_bulk_R_m2K_W=Rd,equivalent_R_contact_for_entire_gain_m2K_W=R_equiv,required_mean_contact_jump_K=q*R_equiv,steady_Tdeep_based_R_upper_if_ice0_and_clearance0_m2K_W=R_sensor_zero_clearance,steady_Tdeep_based_R_upper_with_plus04C_sensor_bound_m2K_W=R_sensor_upper,source_assumed_sensor_clearance_m=clearance,steady_R_contact_estimate_after_clearance_m2K_W=R_sensor_zero_clearance-clearance/k,illustrative_pure_air_gap_m_if_kair0025=R_equiv*.025,illustrative_radiation_conductance_eps09_at273K_W_m2K=4*.9*5.670374419e-8*273.15**3,spectra=spectra,critical_assumptions=['time-averaged local steady conduction','known co-located sensor geometry','ice side approximately0C','no important net storage/latent-pore/advection contribution'],contact_mechanism_demonstrated=False,novelty_demonstrated=False,Tower2_camera_read=False)
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--root',type=Path,required=True);a.add_argument('--out',type=Path,required=True);x=a.parse_args()
 if x.out.exists():raise ValueError('new output required')
 result=screen(x.root/'data/decoded/debtemp_Tower1.csv',x.root/'model/frozen_piramide_v3/NORMALIZED_HEAT_PREDICTIONS.json',x.root/'model/frozen_piramide_v3/melt_comparison/PREDICTIONS.json');x.out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps({k:result[k] for k in ['equivalent_R_contact_for_entire_gain_m2K_W','required_mean_contact_jump_K','temperature_means_C','steady_Tdeep_based_R_upper_with_plus04C_sensor_bound_m2K_W']},indent=2))
