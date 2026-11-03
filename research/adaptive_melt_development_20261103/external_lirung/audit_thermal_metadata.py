"""Metadata-only thermistor audit. No ablation file request; no values exported."""
import urllib.request,pathlib,hashlib,csv,io,json,datetime,math,collections
URL='https://zenodo.org/api/records/4446387/files/Debris_temperature_profile_Lirung_Glacier_CHAND.csv/content'
def main():
 b=urllib.request.urlopen(URL,timeout=30).read()
 assert len(b)==50324 and hashlib.md5(b).hexdigest()=='65d67be6eaf8e47aa3a492a17e0ac3e7'
 r=list(csv.reader(io.StringIO(b.decode('cp1252'))));groups=[]
 for i in range(0,32,4):
  ts=[];missing=[0]*3
  for row in r[7:]:
   if len(row)<=i or not row[i].strip():continue
   try:t=datetime.datetime.strptime(row[i].strip(),'%m/%d/%Y %H:%M')
   except ValueError:continue
   ts.append(t)
   for j in range(3):
    try:ok=math.isfinite(float(row[i+j+1]))
    except (ValueError,IndexError):ok=False
    missing[j]+=not ok
  groups.append({'group':i//4+1,'site_header':r[5][i:i+4], 'depth_headers':r[6][i+1:i+4], 'n_clock_rows':len(ts),'start':min(ts).isoformat() if ts else None,'end':max(ts).isoformat() if ts else None,'cadence_seconds_counts':dict(collections.Counter(str(int((v-u).total_seconds())) for u,v in zip(ts,ts[1:]))),'missing_per_temperature_column':missing})
 out={'url':URL,'bytes':len(b),'md5':hashlib.md5(b).hexdigest(),'sha256':hashlib.sha256(b).hexdigest(),'encoding':'cp1252','groups':groups,'native_air_file_available_in_record':False,'outcome_file_requested':False,'thermal_values_exported':False,'clock':'Literal source local date labels; timezone undocumented; do not infer UTC.'}
 p=pathlib.Path(__file__).with_name('THERMAL_COVERAGE_ONLY.json');p.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
