"""Source-resolved historical baseline check; refuses to claim an unmatched score."""
import argparse,csv,datetime as dt,hashlib,io,json,math,re,zipfile
from pathlib import Path
import xml.etree.ElementTree as ET
NS={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
def xlsx_rows(path):
 with zipfile.ZipFile(path) as z:
  ss=[''.join(e.itertext()) for e in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('m:si',NS)] if 'xl/sharedStrings.xml' in z.namelist() else []
  root=ET.fromstring(z.read('xl/worksheets/sheet1.xml'));out=[]
  for row in root.findall('.//m:row',NS):
   record={}
   for c in row.findall('m:c',NS):
    v=c.find('m:v',NS);x=v.text if v is not None else ''
    record[c.attrib['r'].rstrip('0123456789')]=ss[int(x)] if c.attrib.get('t')=='s' and x else x
   out.append(record)
  return out

def percent_error(model,observed):
 if any(type(x) not in (float,int) or not math.isfinite(x) for x in [model,observed]) or model<0 or observed<=0:raise ValueError('finite nonnegative model and positive observation required')
 return 100*(model-observed)/observed

def calculate(protocol,source_zip,outputs):
 run=next(r for r in xlsx_rows(protocol/'MMIC_model_run_dates.xlsx') if r.get('A')=='CN-5470_20152017')
 start=dt.datetime(*[int(run[k]) for k in ['F','G','H','I']]);end=dt.datetime(*[int(run[k]) for k in ['J','K','L','M']])
 expected=int(run['C'])-int(run['B'])+1
 model_path=outputs/'DEBCF/CN/DEBCF_CN_EXP1_EXP2.csv'
 matrix=[[float(x) for x in row] for row in csv.reader(model_path.read_text().splitlines())]
 if len(matrix)!=expected or any(len(r)!=13 or not all(math.isfinite(x) and x>=0 for x in r) for r in matrix):raise ValueError('actual-thickness matrix size/values inconsistent with protocol')
 dates_path=outputs/'Dates_Rows/Dates_hourly_CN.csv';date_text=dates_path.read_text().splitlines()
 # The exact archive labels are recorded; these bounds follow the explicit run
 # spreadsheet, not a fit-selected trimming of the output matrix.
 if len(date_text)!=expected:raise ValueError('output date count mismatch')
 month={k:i for i,k in enumerate(['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'],1)}
 def output_stamp(text):
  m=re.fullmatch(r'(\d{2})-([A-Za-z]{3})-(\d{4}) (\d{2}):(\d{2}):(\d{2})',text.strip())
  if m:
   day,mon,year,h,minute,sec=m.groups();return dt.datetime(int(year),month[mon],int(day),int(h),int(minute),int(sec))
  t=dt.datetime.fromisoformat(text)
  if t.tzinfo is not None:raise ValueError('archive timestamps must match source-naive local convention')
  return t
 stamps=[output_stamp(x) for x in date_text]
 if stamps!=[start+dt.timedelta(hours=i) for i in range(expected)] or stamps[-1]!=end:raise ValueError('output timestamps do not match source-prescribed hourly grid')
 with zipfile.ZipFile(source_zip) as z:
  forcing=z.read('CN-5470_20152017/CN-5470_20152017_60.csv')
  forcing_rows=list(csv.reader(io.StringIO(forcing.decode()),delimiter=';'))[1:]
  def stamp(row):return dt.datetime(*map(int,row[:3]))+dt.timedelta(hours=int(row[3]),minutes=int(row[4]))
  if stamp(forcing_rows[int(run['B'])-1])!=start or stamp(forcing_rows[int(run['C'])-1])!=end:raise ValueError('forcing indices do not match published run dates')
  raw=z.read('CN-5470_20152017/CN-5470_20152017_VAL_dailysurfaceheight.csv')
  rows=list(csv.reader(io.StringIO(raw.decode()),delimiter=';'))[2:]
  selected=[r for r in rows if start.date()<=dt.date(*map(int,r[:3]))<=end.date()]
  dates=[dt.date(*map(int,r[:3])) for r in selected]
  if len(set(dates))!=len(dates) or dates!=[start.date()+dt.timedelta(days=i) for i in range((end.date()-start.date()).days+1)]:raise ValueError('daily observation support incomplete')
  values=[float(r[4]) for r in selected]
  if not all(math.isfinite(v) and v>=0 for v in values):raise ValueError('daily ABL missing/invalid')
 model=sum(r[0] for r in matrix)
 # Protocol FINAL_RUNS explicitly prescribes m w eq and ice density900kg/m3.
 # The source ABL column is thresholded lowering(m), not assumed pure ice melt.
 lowering=sum(values);observed=lowering*.9;score=percent_error(model,observed)
 table=next(r for r in xlsx_rows(protocol/'published_table3.xlsx') if r.get('B')=='DEBCF')
 published=float(table['D'].replace('−','-'));match=abs(score-published)<=.05
 return dict(model='DEBCF',site='Changri Nup',actual_thickness_column=1,hourly_rows=expected,run_local_start=start.isoformat(),run_local_end=end.isoformat(),forcing_begin_index=int(run['B']),forcing_end_index=int(run['C']),daily_calendar_observations=len(values),model_sum_m_water_equivalent=model,thresholded_surface_lowering_m=lowering,protocol_ice_density_kg_m3=900,calendar_ABL_sum_as_m_water_equivalent=observed,conditional_final_cumulative_error_percent=score,published_Table3_percent=published,published_rounding_tolerance_percentage_points=.05,absolute_discrepancy_percentage_points=abs(score-published),published_score_reproduced=match,independent_validation=False,observation_alignment_certified=False,interpretation='Conditional raw-calendar ABL reproduction using source-prescribed units/density. Exact processed observation and daily interval aggregation remain unresolved; no alternative conversion is chosen by fit.')

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--protocol-dir',type=Path,required=True);p.add_argument('--source-data-zip',type=Path,required=True);p.add_argument('--selected-outputs',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
 if a.out.exists():raise ValueError('new output path required')
 r=calculate(a.protocol_dir,a.source_data_zip,a.selected_outputs)
 r['input_sha256']={str(x):hashlib.sha256(x.read_bytes()).hexdigest() for x in [a.protocol_dir/'MMIC_model_run_dates.xlsx',a.protocol_dir/'published_table3.xlsx',a.source_data_zip,a.selected_outputs/'DEBCF/CN/DEBCF_CN_EXP1_EXP2.csv',a.selected_outputs/'Dates_Rows/Dates_hourly_CN.csv']}
 r['executed_source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
 a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(r,indent=2,allow_nan=False)+'\n');print(json.dumps(r,indent=2))
 raise SystemExit(0 if r['published_score_reproduced'] and r['observation_alignment_certified'] else 2)
if __name__=='__main__':main()
