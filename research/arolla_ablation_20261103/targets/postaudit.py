"""Independent raw-source audit; does not import or recalibrate production models."""
import csv,datetime as d,hashlib,io,json,math,pathlib,xml.etree.ElementTree as E,zipfile
P=pathlib.Path('/private/tmp/melt-validation-20261103');M=P/'model';T=P/'targets';tz=d.timezone(d.timedelta(hours=2));ns={'x':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
r=json.loads((M/'results/frozen_ablation_comparison/RESULTS.json').read_text());freeze=json.loads((M/'FREEZE.json').read_text())
assert r['preexecution_commit']=='bc7916b852b9667820ef68e7442d4759e06ff342'
assert all(hashlib.sha256((M/n).read_bytes()).hexdigest()==h for n,h in freeze['source_sha256'].items())
assert hashlib.sha256((T/'AROLLA_PARSED.json').read_bytes()).hexdigest()==freeze['data_sha256']
with zipfile.ZipFile(T/'1_ARO_2010.zip') as z:
 with zipfile.ZipFile(io.BytesIO(z.read('ARO_2010/ARO_2010_VAL_stakes.xlsx'))) as w:
  cells={c.attrib['r']:c.find('x:v',ns).text for c in E.fromstring(w.read('xl/worksheets/sheet1.xml')).findall('.//x:c',ns) if c.find('x:v',ns)is not None}
 weather=list(csv.DictReader(io.StringIO(z.read('ARO_2010/ARO_2010_60.csv').decode())))
 rho=float(cells['J4']);assert rho==900
 targets={}
 for stake,row,rawrow,start in [('H1D',16,4,'D'),('G1D',17,5,'D'),('F1D',18,6,'E'),('E1D',19,7,'E')]:
  c=[start,'F','G'];vals=[float(cells[x+str(row)]) for x in c];reads=[float(cells[x+str(rawrow)]) for x in c]
  a,b=vals[1]-vals[0],vals[2]-vals[1]
  assert abs(a-(reads[1]-reads[0])*.01*rho)<1e-9 and abs(b-(reads[2]-reads[1])*.01*rho)<1e-9
  dates=[(d.datetime(1899,12,30,tzinfo=tz)+d.timedelta(days=float(cells[x+'2']))).date().isoformat() for x in c]
  targets[stake]=(a,b,dates)
 forcing=[]
 for x in weather:
  hour=int(x['HOUR']);t=d.datetime(int(x['YEAR']),int(x['MONTH']),int(x['DAY']),tzinfo=tz)+d.timedelta(hours=hour,minutes=int(x['MINUTE']))
  forcing.append((t,x))
assert all(b[0]-a[0]==d.timedelta(hours=1) for a,b in zip(forcing,forcing[1:]))
def dose(a,b,offset,kind):
 total=seconds=0.
 for label,x in forcing:
  hi=label+d.timedelta(hours=offset);lo=hi-d.timedelta(hours=1);s=max(0.,(min(hi,b)-max(lo,a)).total_seconds())
  if s:
   seconds+=s
   if kind=='constant_rate':v=1/86400
   elif kind=='positive_degree_day':v=max(float(x['T']),0)/86400
   else:v=max(float(x['SWIN'])-float(x['SWOUT']),0)/1e6
   assert math.isfinite(v);total+=v*s
 assert abs(seconds-(b-a).total_seconds())<1e-7
 return total
pub=list(csv.reader((M/'arolla_outputs/DEBCF/ARO/DEBCF_ARO_EXP1_EXP2.csv').read_text().splitlines()));labels=[d.datetime.strptime(x,'%d-%b-%Y %H:%M:%S').replace(tzinfo=tz) for x in (M/'arolla_outputs/Dates_Rows/Dates_hourly_ARO.csv').read_text().splitlines()]
assert len(pub)==len(labels)==1686
maximum=0.;seen=set()
for x in r['rows']:
 stake=x['stake'];a,b,dates=targets[stake];hours=x['visit_hours'];offset=x['mean_label_end_offset_h'];key=(stake,x['model'],offset,tuple(hours));assert key not in seen;seen.add(key)
 times=[d.datetime.combine(d.date.fromisoformat(date),d.time(),tz)+d.timedelta(hours=h) for date,h in zip(dates,hours)]
 assert x['heldout_observed_mmwe']==b and x['test_interval']==[q.isoformat() for q in times[1:]]
 if x['model'].startswith('published'):
  assert stake=='H1D' and x['local_calibration_applied'] is False and x['pure_ice_component_match_established'] is False
  expected=seconds=0.
  for label,y in zip(labels,pub):
   hi=label+d.timedelta(hours=offset);lo=hi-d.timedelta(hours=1);s=max(0.,(min(hi,times[2])-max(lo,times[1])).total_seconds());seconds+=s;expected+=float(y[0])*1000*s/3600
  assert seconds==(times[2]-times[1]).total_seconds()
 else:
  assert x['calibration_observed_mmwe']==a and x['train_interval']==[q.isoformat() for q in times[:2]]
  train=dose(times[0],times[1],offset,x['model']);test=dose(times[1],times[2],offset,x['model']);coef=a/train;expected=coef*test
  assert math.isclose(coef,x['fitted_training_scale'],rel_tol=2e-13)
  assert math.isclose(train,x['train_dose'],rel_tol=2e-13) and math.isclose(test,x['test_dose'],rel_tol=2e-13)
 delta=abs(expected-x['predicted_mmwe']);maximum=max(maximum,delta);assert delta<1e-8
 assert abs(x['signed_error_mmwe']-(expected-b))<1e-8 and abs(x['absolute_error_mmwe']-abs(expected-b))<1e-8
 assert x['uncertainty']['confidence_interval'] is False
assert len(seen)==1053
nom=[x for x in r['rows'] if x['nominal']];mae={m:sum(x['absolute_error_mmwe'] for x in nom if x['model']==m)/sum(x['model']==m for x in nom) for m in r['model_summary']}
for m,v in mae.items():assert abs(v-r['model_summary'][m]['available_stake_nominal_mae_mmwe'])<1e-10
report={'all_1053_rows_independently_reproduced_from_raw_sources':True,'maximum_prediction_difference_mmwe':maximum,'frozen_hashes_match':True,'native_workbook_and_reading_increment_agree':True,'first_interval_only_coefficients_verified':True,'full_clock_scenario_support_verified':True,'nominal_rows':nom,'all_available_stake_nominal_mae_mmwe':mae,'scope':'Retrospective chronological native ice-ablation holdout, one shared weather event. Not coupled cliff validation or novel/professionally superior model. Professional output quantity equivalence unresolved.','original_results_sha256':hashlib.sha256((M/'results/frozen_ablation_comparison/RESULTS.json').read_bytes()).hexdigest()}
(T/'POSTSCORE_REVIEW.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({'rows':len(seen),'max_prediction_difference':maximum,'mae':mae,'result_sha256':report['original_results_sha256']}))
