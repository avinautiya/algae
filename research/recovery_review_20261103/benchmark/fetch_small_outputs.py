"""Bounded HTTP range fetch of three published members, not a full ZIP download."""
import hashlib,json,struct,urllib.request,zlib
from pathlib import Path
ROOT=Path(__file__).resolve().parent
URL='https://zenodo.org/api/records/15754456/files/DCG-MIP_Model_Outputs.zip/content'
WANTED=['Dates_Rows/Dates_daily_CN.csv','Dates_Rows/Dates_hourly_CN.csv','DEBCF/CN/DEBCF_CN_EXP1_EXP2.csv']
def get(start,n):
 if n>1000000:raise ValueError('request exceeds bounded member limit')
 r=urllib.request.urlopen(urllib.request.Request(URL,headers={'Range':f'bytes={start}-{start+n-1}'}),timeout=40)
 if r.status!=206 or r.headers.get('Content-Range')!=f'bytes {start}-{start+n-1}/1222386989':r.close();raise ValueError('range response mismatch')
 b=r.read(n+1);r.close()
 if len(b)!=n:raise ValueError('range length mismatch')
 return b
def main():
 entries=json.loads((ROOT/'ZIP_MEMBERS.json').read_text());out=[]
 for name in WANTED:
  e=next(x for x in entries if x['name']==name)
  if e['size']>1000000 or e['method']!=8:raise ValueError('unsupported member')
  h=struct.unpack('<4s5H3I2H',get(e['offset'],30))
  if h[0]!=b'PK\x03\x04':raise ValueError('bad local header')
  raw=get(e['offset']+30+h[9]+h[10],e['compressed']);data=zlib.decompress(raw,-15)
  if len(data)!=e['size'] or zlib.crc32(data)!=e['crc']:raise ValueError('member CRC/length mismatch')
  dest=ROOT/'selected_outputs'/name;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(data)
  out.append(dict(member=name,sha256=hashlib.sha256(data).hexdigest(),size=len(data),crc32_verified=True))
 (ROOT/'SELECTED_OUTPUTS_MANIFEST.json').write_text(json.dumps(dict(record=15754456,archive_md5_not_fully_verified=True,method='bounded HTTPS byte ranges plus member CRC32; no full archive MD5 claim',members=out),indent=2)+'\n')
 print(json.dumps(out,indent=2))
if __name__=='__main__':main()
