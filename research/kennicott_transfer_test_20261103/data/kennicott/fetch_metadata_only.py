"""Allowlisted ZIP text metadata only. Rejects any target/CSV body request."""
from list_remote_zip import get_range
import pathlib,json,struct,zlib,hashlib
P=pathlib.Path(__file__).parent
ALLOW=['Debris_Temperature_Profiles/Metadata for Debris_Temperature_Profiles.txt','Point_Mass_Balance/point_mass_balance_metadata.txt','AWS/Debris_T_metadata.txt','AWS/Debris_EB_2020_metadata.txt','AWS/Amazon_Creek_metadata.txt']
def main():
 m=json.loads((P/'ZIP_METADATA_ONLY.json').read_text());manifest=[]
 for suffix in ALLOW:
  e=next(x for x in m['entries'] if x['name'].endswith('/'+suffix) and not x['name'].startswith('__MAC'))
  o=e['local_offset'];h=get_range(m['url'],o,o+29);f=struct.unpack('<4s5H3L2H',h);assert f[0]==b'PK\x03\x04';start=o+30+f[-2]+f[-1]
  b=get_range(m['url'],start,start+e['compressed_bytes']-1);b=zlib.decompress(b,-15) if e['compression']==8 else b
  assert len(b)==e['uncompressed_bytes'] and zlib.crc32(b)&0xffffffff==e['crc32']
  name=pathlib.Path(suffix).name;(P/name).write_bytes(b);manifest.append({'member':e['name'],'file':name,'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest(),'crc32':e['crc32']});print('\nFILE',name);print(b.decode('utf-8-sig',errors='replace'))
 (P/'METADATA_FILES_MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
if __name__=='__main__':main()
