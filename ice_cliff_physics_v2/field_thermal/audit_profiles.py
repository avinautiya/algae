"""Audit public thermal geometry/time support without fitting or scoring models."""
import csv,datetime as dt,hashlib,json,math,re,statistics
from pathlib import Path
ROOT=Path(__file__).resolve().parent
SOURCE=ROOT/'source'

def num(s):
 try:v=float(s)
 except (ValueError,TypeError):return None
 return v if math.isfinite(v) else None

def parse(s):return dt.datetime.fromisoformat(s)

def readcsv(p):
 with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))

def main():
 folder=SOURCE/'Debris_Temperature_Profiles';summary=readcsv(folder/'debris_thermistors_profiles_summary.csv');audit=[];matches=[]
 mass=readcsv(SOURCE/'Point_Mass_Balance/point_mass_balance.csv')
 for meta in summary:
  name=meta['Profile'];p=folder/'Profile_Temperature_Records'/f'{name}.csv';rows=readcsv(p)
  cols=[k for k in rows[0] if k.startswith('Depth:')];depths=[float(re.search(r'Depth:\s*([\d.]+)',c).group(1)) for c in cols]
  depthfile=[float(r['Depth']) for r in readcsv(folder/'Profile_Sensor_Depths'/f'{name}_Depth.csv')]
  times=[parse(r['Datetime']) for r in rows];intervals=[(b-a).total_seconds() for a,b in zip(times,times[1:])]
  if len(times)!=len(set(times)) or any(v<=0 for v in intervals):raise ValueError(name+' unordered/duplicate time')
  values=[[num(r[col]) for col in cols] for r in rows]
  fields=[]
  for j,(c,z) in enumerate(zip(cols,depths)):
   valid=[r[j] for r in values if r[j] is not None]
   fields.append(dict(column=c,depth_cm=z,is_observed_probe=not(name=='KDT_2020b' and z==22.),observation_status='author_inserted_unsensed_zero_C_boundary' if name=='KDT_2020b' and z==22. else 'measured_thermistor_per_source_sensor_count',valid_count=len(valid),missing_count=len(rows)-len(valid),negative_count=sum(v<0 for v in valid),min_C=min(valid) if valid else None,max_C=max(valid) if valid else None))
  nominal=statistics.median(intervals);longest=run=0;longest_start=start=None
  for i,(t,vs) in enumerate(zip(times,values)):
   complete=all(v is not None for v in vs)
   if not complete:run=0;start=None;continue
   if i and run and (t-times[i-1]).total_seconds()<=1.01*nominal:run+=1
   else:run=1;start=t
   if run>longest:longest=run;longest_start=start
  shifts=sorted(set((parse(r['Datetime'])-parse(r['Datetime_ADT'])).total_seconds()/3600 for r in rows if r.get('Datetime_ADT')))
  d=dict(profile=name,source_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),rows=len(rows),utc_start=times[0].isoformat(),utc_end=times[-1].isoformat(),span_days=(times[-1]-times[0]).total_seconds()/86400,
   cadence_seconds=dict(min=min(intervals),median=nominal,max=max(intervals),unique=sorted(set(intervals))),
   depth_columns_cm=depths,depth_file_cm=depthfile,depth_header_matches_file=depths==depthfile,summary_N=meta['N'],summary_debris_depth_cm=float(meta['Debris depth (cm)']),
   shallowest_depth_cm=min(depths),deepest_depth_cm=max(depths),strictly_interior_probe_depths_cm=depths[1:-1],
   temperature_accuracy_C=.2,utc_minus_ADT_hours=shifts,sensors=fields,
   longest_all_probes_finite_run=dict(rows=longest,start_utc=longest_start.isoformat() if longest_start else None,duration_days=max(0,(longest-1))*nominal/86400),
   all_probes_finite_rows=sum(all(v is not None for v in vs) for vs in values),all_probes_nonnegative_finite_rows=sum(all(v is not None and v>=0 for v in vs) for vs in values),
   author_linked_ablation_stake=meta['Ablation Stake'] or None,
   coordinate=dict(lat=float(meta['lat']),lon=float(meta['lon']),elevation_m=float(meta['ele'])),
   published_depth_averaged_k_not_independent_prior=float(meta['Depth-averaged k']))
  observed_depths=[v['depth_cm'] for v in fields if v['is_observed_probe']]
  d.update(observed_depth_columns_cm=observed_depths,actual_measured_sensor_count=len(observed_depths),primary_measured_endpoint_validation_eligible=len(observed_depths)>=3,primary_interior_target_depths_cm=observed_depths[1:-1])
  audit.append(d)
  if meta['Ablation Stake']:
   for row in mass:
    if row['site']!=meta['Ablation Stake']:continue
    a=dt.datetime.strptime(row['date_start'],'%m/%d/%Y') if row['date_start']!='unknown' else None
    b=dt.datetime.strptime(row['date_end'],'%m/%d/%Y')
    # Date-only overlap is descriptive, not a conversion to exact UTC visit hours.
    if a is not None and a.date()<=times[-1].date() and b.date()>=times[0].date():
     match=dict(profile=name,stake=row['site'],date_start=row['date_start'],date_end=row['date_end'],mass_change_kg_m2=float(row['mass_change_kg_m2']),stake_debris_thickness_cm=float(row['debris_thick_mm'])/10,profile_debris_thickness_cm=d['summary_debris_depth_cm'],temperature_start_utc=d['utc_start'],temperature_end_utc=d['utc_end'],mass_date_clock='date-only, timezone not established here',source_link='explicit summary Ablation Stake ID, not nearest-neighbor assignment',raw_point_row=row)
     if row['lat_start'] and row['lon_start']:
      # Metadata says degrees WEST for stake longitudes, unlike signed profile longitude.
      lat=float(row['lat_start']);lon=-abs(float(row['lon_start']));dl=math.radians(lat-d['coordinate']['lat']);do=math.radians(lon-d['coordinate']['lon']);av=math.radians((lat+d['coordinate']['lat'])/2)
      match['initial_coordinate_distance_m_approx']=6371000*math.sqrt(dl*dl+(math.cos(av)*do)**2)
     matches.append(match)
 result=dict(scope='data feasibility and support audit only; no model fitting or predictive scores',profiles=audit,explicit_ablation_matches=matches,identifiability='Temperature diffusion can identify effective alpha=k/C only. Scaling k and C together preserves temperature but changes heat flux; melt flux needs independent material/ablation constraints.')
 (ROOT/'PROFILE_FEASIBILITY.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
 print(json.dumps([dict(profile=x['profile'],rows=x['rows'],start=x['utc_start'],end=x['utc_end'],depths=x['depth_columns_cm'],missing=[s['missing_count'] for s in x['sensors']],longest_days=x['longest_all_probes_finite_run']['duration_days']) for x in audit],indent=2))
if __name__=='__main__':main()
