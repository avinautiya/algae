"""Thermal and weather file coverage/missingness only; no outcomes or temperature values exported."""
from audit_paired_metadata import fetch,P
import csv,io,json,hashlib,collections,math,datetime
FILES=['Debris_Temperature_Profiles/Profile_Temperature_Records/KDT_2023a.csv','Debris_Temperature_Profiles/Profile_Temperature_Records/KDT_2023b.csv','Debris_Temperature_Profiles/Profile_Temperature_Records/KDT_2023c.csv','AWS/Debris_T2.csv','AWS/Debris_T3.csv']
def main():
 out=[]
 for s in FILES:
  b,e=fetch(s);rs=list(csv.DictReader(io.StringIO(b.decode('utf-8-sig'))));headers=list(rs[0]);key=next(k for k in headers if k.lower().replace('_','') in ['datetime','datetimeutc']);ts=[];miss={k:0 for k in headers if 'date' not in k.lower()}
  for row in rs:
   try:t=datetime.datetime.fromisoformat(row[key])
   except ValueError:
    try:t=datetime.datetime.strptime(row[key],'%m/%d/%Y %H:%M')
    except ValueError:continue
   ts.append(t)
   for k in miss:
    try:ok=math.isfinite(float(row[k]))
    except ValueError:ok=False
    miss[k]+=not ok
  item={'member':s,'crc32':e['crc32'],'sha256':hashlib.sha256(b).hexdigest(),'headers':headers,'clock_key':key,'rows':len(rs),'parseable_clock_rows':len(ts),'first':min(ts).isoformat() if ts else None,'last':max(ts).isoformat() if ts else None,'cadence_seconds_counts':dict(collections.Counter(str((v-u).total_seconds()) for u,v in zip(ts,ts[1:]))),'missing_numeric_counts':miss,'outcome_data_accessed':False,'feature_values_exported':False};out.append(item);print(json.dumps(item,indent=2),flush=True)
 (P/'FEATURE_CLOCK_AUDIT.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
