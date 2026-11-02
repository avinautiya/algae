"""Inspect public original records without executing author processing code."""
import csv,datetime as dt,json,math
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def main():
 out=[];series={}
 for name,depth in [('20457129.csv',5),('20457130.csv',15)]:
  with (ROOT/'original_2022'/name).open(encoding='utf-8-sig',newline='') as f:next(f);rows=list(csv.DictReader(f))
  tk=next(k for k in rows[0] if k.startswith('Date Time'));vk=next(k for k in rows[0] if k.startswith('Temp,'))
  data={dt.datetime.strptime(r[tk],'%m/%d/%y %H:%M:%S')+dt.timedelta(hours=8):float(r[vk]) for r in rows};ts=sorted(data);series[depth]=data
  out.append(dict(file=name,assumed_july_depth_cm=depth,raw_clock_header=tk,rows=len(rows),start_UTC=ts[0].isoformat(),end_UTC=ts[-1].isoformat(),max_gap_seconds=max((b-a).total_seconds() for a,b in zip(ts,ts[1:])),depth_status='July/August depth configuration exists; exact transition date from June depths unverified'))
 with (ROOT/'source/Debris_Temperature_Profiles/Profile_Temperature_Records/KDT_2020b.csv').open(encoding='utf-8-sig') as f:rows=list(csv.DictReader(f))
 compare={}
 for depth,data in series.items():
  errors=[];missing=0
  for r in rows:
   t=dt.datetime.fromisoformat(r['Datetime'])
   if t not in data:missing+=1;continue
   errors.append(abs(float(r[f'Depth: {depth} cm'])-data[t]))
  compare[str(depth)]=dict(joined=len(errors),missing=missing,max_abs_difference_C=max(errors),mean_abs_difference_C=math.fsum(errors)/len(errors))
 zero=all(float(r['Depth: 22 cm'])==0 for r in rows)
 result=dict(scope='original-source data/clock comparison, not model accuracy',profiles=out,overlap_processed_comparison=compare,synthetic_boundary=dict(profile='KDT_2020b',column='Depth: 22 cm',summary_real_sensors=2,all_values_exact_zero=zero,confirmed_author_profile_info='20457129 at5cm,20457130 at15cm,Base at22cm',author_code_behavior='Base constructed as shallow temperature times zero; documented as a non-sensor ice-interface boundary',source_code_url='https://github.com/ericivanpetersen/supraglacial_debris_heat_flux_jgr/blob/4863e0bec8e496452421f97d616b03925451f215/plot_debris_temp_paper.py#L56-L86'))
 (ROOT/'ORIGINAL_2022_SOURCE_AUDIT.json').write_text(json.dumps(result,indent=2)+'\n')
 print('Source comparison complete; no predictive models fitted or scored')
if __name__=='__main__':main()
