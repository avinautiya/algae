"""Fetch only selected small ZIPmembers; verify CRC32/size, never execute source."""
import hashlib,json,struct,urllib.request,zlib
from pathlib import Path
ROOT=Path(__file__).parent
INDEX=ROOT/'ARCHIVE_INDEX.json'
INDEX_SHA256='149537e63f200115c655de7aaa779483ade2e252d5c1587cca9272adb8ad34f1'
URL='https://zenodo.org/api/records/11581073/files/Debris_properties_estimations.zip/content'
NAMES=['run_data_Piramide_Tower1_26092023_01-Nov-2022-24-Apr-2023.mat','run_data_Piramide_Tower2_26092023_01-Nov-2022-24-Apr-2023.mat','PIR22-23_60-kd_ERA5.csv','z0Calc.m','TC_Piramide_Calibration.m','ShortwaveFluxesOS.m','Surface_Temperature_Debris.m','MAIN_FRAME.m']
def fetch(a,b):
 req=urllib.request.Request(URL,headers={'Range':f'bytes={a}-{b}','User-Agent':'glacier-source-audit/1.0'})
 with urllib.request.urlopen(req,timeout=60)as r:
  if r.status!=206 or not r.headers.get('Content-Range','').startswith(f'bytes {a}-{b}/'):raise ValueError('exact ranged response required')
  return r.read()
def main():
 if hashlib.sha256(INDEX.read_bytes()).hexdigest()!=INDEX_SHA256:raise ValueError('sourcearchiveindexhashchanged')
 manifest=[]
 for name in NAMES:
  rec=next(r for r in json.loads(INDEX.read_text())if r['name'].endswith('/'+name));head=fetch(rec['offset'],rec['offset']+29);fields=struct.unpack('<4s5H3I2H',head)
  if fields[0]!=b'PK\x03\x04':raise ValueError('local ZIPsignature')
  begin=rec['offset']+30+fields[-2]+fields[-1];compressed=fetch(begin,begin+rec['compressed_size']-1)
  raw=zlib.decompress(compressed,-15)if rec['method']==8 else compressed
  if len(raw)!=rec['size']or zlib.crc32(raw)&0xffffffff!=rec['crc32']:raise ValueError('member integrity failed')
  path=ROOT/name
  if path.exists()and path.read_bytes()!=raw:raise ValueError('refuse changed source overwrite')
  path.write_bytes(raw);manifest.append({'archive_record':'10.5281/zenodo.11581073','member':rec['name'],'size':len(raw),'crc32':rec['crc32'],'sha256':hashlib.sha256(raw).hexdigest(),'path':str(path),'whole_archive_md5_verified':False,'executed_source':False})
 with (ROOT/'REPRODUCED_DOWNLOAD_MANIFEST.json').open('x')as out:out.write(json.dumps(manifest,indent=2)+'\n')
 print(json.dumps(manifest,indent=2))
if __name__=='__main__':main()
