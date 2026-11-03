"""Feature-only source contract; target numerical values never imported."""
import pathlib,csv,json,datetime,math,hashlib
P=pathlib.Path(__file__).parent
PAIRS=[('KDT_2020b','P_AWS_1','train','2020-07-25','2020-09-06',.22,.21,[.05,.15,.22],'Debris_T3'),('KDT_2023a','P11','evaluation','2023-06-08','2023-09-02',.24,.24,[.02,.07,.12,.17,.22],'Debris_T2'),('KDT_2023b','P12','evaluation','2023-06-08','2023-09-02',.32,.29,[.09,.165,.24,.315],'Debris_T2'),('KDT_2023c','P13','evaluation','2023-06-08','2023-09-02',.26,.255,[.02,.095,.17,.245],'Debris_T2')]
def parse(s):return datetime.datetime.fromisoformat(s)
def finite(s):
 try:return math.isfinite(float(s))
 except (ValueError,TypeError):return False
def main():
 out={'scope':'Source/feature contract only. Candidate proposed gap scenes not experimentally selected or target-tested.','ice_density_kg_m3':900,'ice_density_provenance':'Petersen2025 source assumption for reported mass change','target_quantity':'signed mass_change_kg_m2 equals mmwater equivalent; negative means loss','target_conversion':'positive ablation = -signed mass change only for declared no-other-mass-change/snow-accounting scenario','source_stake_dates':'month/day/year date strings; no native visit hour or timezone given','nominal_clock_assumption':'sourceADT(UTC-8) noon on each visit date; NOTmeasured visit times','thermal_accuracy_absolute_C':.2,'depth_accuracy_m':None,'published_k_or_pooled_capacity_allowed':False,'pairs':[],'weather':[]}
 for prof,stake,role,d0,d1,dp,ds,z,air in PAIRS:
  raw=(P/'features'/f'{prof}.csv').read_bytes();r=list(csv.DictReader(raw.decode('utf-8-sig').splitlines()));t0=parse(d0)+datetime.timedelta(hours=20);t1=parse(d1)+datetime.timedelta(hours=20);first=parse(r[0]['Datetime']);reset_end=first+datetime.timedelta(days=5)
  out['pairs'].append({'profile_id':prof,'stake_id':stake,'role':role,'target_start_date':d0,'target_end_date':d1,'profile_path':'features/'+prof+'.csv','profile_sha256':hashlib.sha256(raw).hexdigest(),'sensor_depths_m':z,'profile_debris_thickness_m':dp,'stake_debris_thickness_m':ds,'forward_upper_boundary_plane_m':z[0],'forward_domain_length_stake_m':ds-z[0],'preferred_native_air_station_by_proximity':air,'nominal_start_UTC':t0.isoformat(),'nominal_end_UTC':t1.isoformat(),'first_profile_UTC':first.isoformat(),'nominal_initial_no_profile_hours':(first-t0).total_seconds()/3600,'proposed_installation_exclusion_hours_for_parameter_fit_only':120,'proposed_profile_fit_start_UTC':reset_end.isoformat(),'proposed_proxy_fit_window_days':14,'source_native_deepest_in_stake_domain':z[-1]<=ds,'warnings':(['last22cmcolumn identicallyzero; not independently observed thirdsensor'] if prof=='KDT_2020b' else [])+(['deepest31.5cm exceeds29cmstake thickness; use onlyinprofile inference, neverinside stakecolumn'] if prof=='KDT_2023b' else [])})
 for name,cols in [('Debris_T2',['T_1_75m']),('Debris_T3',['T_1_75m']),('Amazon_Creek',['T_3m','SWin','LWin'])]:
  raw=(P/'features'/f'{name}.csv').read_bytes();rows=list(csv.DictReader(raw.decode('utf-8-sig').splitlines()))
  for year,d0,d1 in [(2020,'2020-07-25T20:00:00','2020-09-06T20:00:00'),(2023,'2023-06-08T20:00:00','2023-09-02T20:00:00')]:
   rr=[r for r in rows if parse(d0)<=parse(r['Date_Time'])<=parse(d1)]
   out['weather'].append({'file':'features/'+name+'.csv','year':year,'nominal_clock':True,'rows_inclusive':len(rr),'expected_rows_inclusive':int((parse(d1)-parse(d0)).total_seconds()/3600)+1,'first':rr[0]['Date_Time'] if rr else None,'last':rr[-1]['Date_Time'] if rr else None,'missing_counts':{c:sum(not finite(r[c]) for r in rr) for c in cols}})
 (P/'SOURCE_BINDING.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'pairs':out['pairs'],'weather':out['weather']},indent=2))
if __name__=='__main__':main()
