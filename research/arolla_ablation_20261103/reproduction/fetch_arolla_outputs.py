"""Fetch exact Arolla members by bounded byte ranges; never download1.22GB ZIP."""
import hashlib,json,struct,urllib.request,zlib
from pathlib import Path
ROOT=Path(__file__).resolve().parent
URL='https://zenodo.org/api/records/15754456/files/DCG-MIP_Model_Outputs.zip/content'
TOTAL=1222386989
def get(start,n):
 if type(start) is not int or type(n) is not int or start<0 or not 0<n<=1000000:raise ValueError('invalid range')
 req=urllib.request.Request(URL,headers={'Range':f'bytes={start}-{start+n-1}'})
 with urllib.request.urlopen(req,timeout=45) as r:
  if r.status!=206 or r.headers.get('Content-Range')!=f'bytes {start}-{start+n-1}/{TOTAL}':raise ValueError('range response mismatch')
  data=r.read(n+1)
 if len(data)!=n:raise ValueError('range length mismatch')
 return data
def main():
 entries=json.loads((ROOT/'AROLLA_MEMBERS.json').read_text());records=[]
 for e in entries:
  if e['method']!=8 or e['size']>1000000:raise ValueError('unsupported member')
  h=struct.unpack('<4s5H3I2H',get(e['offset'],30))
  if h[0]!=b'PK\x03\x04':raise ValueError('bad member header')
  compressed=get(e['offset']+30+h[9]+h[10],e['compressed']);raw=zlib.decompress(compressed,-15)
  if len(raw)!=e['size'] or zlib.crc32(raw)!=e['crc']:raise ValueError('member integrity failure')
  p=ROOT/'arolla_outputs'/e['name'];p.parent.mkdir(parents=True,exist_ok=True)
  if p.exists():
   if p.read_bytes()!=raw:raise ValueError('existing member differs')
  else:p.write_bytes(raw)
  records.append(dict(member=e['name'],sha256=hashlib.sha256(raw).hexdigest(),size=len(raw),crc32=e['crc'],crc32_verified=True))
 result=dict(record=15754456,url=URL,verification='bounded HTTPS ranges and central-directory/member CRC32 and SHA256; full archive MD5 NOT verified',members=records,observations_loaded=False,prediction_scores_computed=False)
 (ROOT/'AROLLA_OUTPUT_MANIFEST.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
