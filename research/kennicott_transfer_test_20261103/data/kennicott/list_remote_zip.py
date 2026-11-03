"""Read ZIP central-directory metadata only; never request outcome file bodies."""
import urllib.request,struct,json,pathlib
URL='https://zenodo.org/records/16389806/files/Petersen_etal_2025_Kennicott_Glacier_Data.zip?download=1'
def get_range(url,begin,end):
 req=urllib.request.Request(url,headers={'Range':f'bytes={begin}-{end}'})
 with urllib.request.urlopen(req,timeout=45) as r:
  if r.status!=206:raise RuntimeError(f'Range refused status{r.status}; no responsebody read')
  return r.read(end-begin+1)
def main():
 req=urllib.request.Request(URL,method='HEAD')
 with urllib.request.urlopen(req,timeout=30) as r:size=int(r.headers['Content-Length'])
 tail=get_range(URL,max(0,size-65536),size-1);i=tail.rfind(b'PK\x05\x06');assert i>=0
 fields=struct.unpack_from('<4s4H2LH',tail,i);n,cdsize,cdoffset=fields[4:7]
 cd=get_range(URL,cdoffset,cdoffset+cdsize-1);entries=[];i=0
 while i<len(cd):
  f=struct.unpack_from('<4s6H3L5H2L',cd,i);assert f[0]==b'PK\x01\x02'
  nl,el,cl=f[10:13];name=cd[i+46:i+46+nl].decode('utf-8',errors='replace');entries.append({'name':name,'crc32':f[7],'compressed_bytes':f[8],'uncompressed_bytes':f[9],'compression':f[4],'local_offset':f[16]});i+=46+nl+el+cl
 out={'url':URL,'archive_bytes':size,'central_directory_bytes':cdsize,'entries':entries,'outcome_file_bodies_read':False};pathlib.Path(__file__).with_name('ZIP_METADATA_ONLY.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
