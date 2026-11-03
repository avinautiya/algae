"""Read only source-binding/date/QC fields; never convert or export mass outcomes."""
import json,pathlib,struct,zlib,csv,io,hashlib
from list_remote_zip import get_range
P=pathlib.Path(__file__).parent
def fetch(suffix):
 m=json.loads((P/'ZIP_METADATA_ONLY.json').read_text());e=next(x for x in m['entries'] if x['name'].endswith('/'+suffix) and not x['name'].startswith('__MAC'));o=e['local_offset'];f=struct.unpack('<4s5H3L2H',get_range(m['url'],o,o+29));a=o+30+f[-2]+f[-1];b=get_range(m['url'],a,a+e['compressed_bytes']-1);b=zlib.decompress(b,-15) if e['compression']==8 else b;assert len(b)==e['uncompressed_bytes'] and zlib.crc32(b)&0xffffffff==e['crc32'];return b,e

def main():
 b,e=fetch('Debris_Temperature_Profiles/debris_thermistors_profiles_summary.csv');rs=list(csv.DictReader(io.StringIO(b.decode('utf-8-sig'))));summary={'source_member':e['name'],'source_sha256':hashlib.sha256(b).hexdigest(),'headers':list(rs[0]),'rows':[{k:v for k,v in r.items() if not any(x in k.lower() for x in ['conduct','diffus','heat','kmelt','averaged k']) and k.lower()!='k'} for r in rs]};(P/'PROFILE_SUMMARY_METADATA_ONLY.json').write_text(json.dumps(summary,indent=2)+'\n');print('PROFILE_SCHEMA',summary)
 b,e=fetch('Point_Mass_Balance/point_mass_balance.csv');rs=list(csv.DictReader(io.StringIO(b.decode('utf-8-sig'))));allow={'site','date_start','date_end','period_length','debris_thick_mm','lat_start','lon_start','elevation_start','surface_type_start','seasonal_snow_depth_start_mm','seasonal_snow_density_start_kg_m3','lat_end','lon_end','elevation_end','surface_type_end','seasonal_snow_depth_end_mm','seasonal_snow_density_end_kg_m3','summer_snow_depth_end_mm','summer_snow_density_end_kg_m3','institution'}
 out={'source_member':e['name'],'source_sha256':hashlib.sha256(b).hexdigest(),'headers':list(rs[0]),'row_count':len(rs),'raw_file_persisted':False,'mass_field_converted_or_exported':False,'rows':[{k:v for k,v in r.items() if k in allow} for r in rs if r.get('site') in {'P11','P12','P13','P_AWS_1'}]};(P/'PAIRED_STAKE_METADATA_ONLY.json').write_text(json.dumps(out,indent=2)+'\n');print('PAIRED_ROWS',json.dumps(out,indent=2))
if __name__=='__main__':main()
