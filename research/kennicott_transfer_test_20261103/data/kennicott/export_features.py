"""Explicitly allowed feature-only member export. Never request point-mass-balance bodies."""
from audit_paired_metadata import fetch,P
import pathlib,csv,io,json,hashlib,datetime,collections,math
ALLOW=['Debris_Temperature_Profiles/Profile_Temperature_Records/KDT_2020a.csv','Debris_Temperature_Profiles/Profile_Temperature_Records/KDT_2020b.csv','Debris_Temperature_Profiles/Profile_Temperature_Records/KDT_2023a.csv','Debris_Temperature_Profiles/Profile_Temperature_Records/KDT_2023b.csv','Debris_Temperature_Profiles/Profile_Temperature_Records/KDT_2023c.csv','AWS/Debris_T2.csv','AWS/Debris_T3.csv','AWS/Amazon_Creek.csv']
def main():
 dest=P/'features';dest.mkdir(exist_ok=True);manifest=[]
 for member in ALLOW:
  name=pathlib.Path(member).name;path=dest/name
  if path.exists():
   b=path.read_bytes();meta=json.loads((P/'ZIP_METADATA_ONLY.json').read_text());e=next(x for x in meta['entries'] if x['name'].endswith('/'+member) and not x['name'].startswith('__MAC'));import zlib;assert len(b)==e['uncompressed_bytes'] and zlib.crc32(b)&0xffffffff==e['crc32']
  else:b,e=fetch(member);path.write_bytes(b)
  rows=list(csv.DictReader(io.StringIO(b.decode('utf-8-sig'))));headers=list(rows[0]);clock=next(k for k in headers if k.lower().replace('_','')=='datetime');times=[]
  for row in rows:
   try:t=datetime.datetime.fromisoformat(row[clock])
   except ValueError:t=datetime.datetime.strptime(row[clock],'%m/%d/%Y %H:%M')
   times.append(t)
  x={'member':member,'local_path':'features/'+name,'sha256':hashlib.sha256(b).hexdigest(),'bytes':len(b),'crc32':e['crc32'],'headers':headers,'rows':len(rows),'clock_field':clock,'clock':'UTC per source metadata','first':times[0].isoformat(),'last':times[-1].isoformat(),'source_outcome_dependency':False}
  manifest.append(x);print(json.dumps({k:v for k,v in x.items() if k not in ['headers']},indent=2),flush=True)
 (P/'FEATURE_EXPORT_MANIFEST.json').write_text(json.dumps({'feature_values_authorized':True,'outcome_values_accessed_or_exported':False,'source_zip_metadata':'ZIP_METADATA_ONLY.json','files':manifest},indent=2)+'\n')
if __name__=='__main__':main()
