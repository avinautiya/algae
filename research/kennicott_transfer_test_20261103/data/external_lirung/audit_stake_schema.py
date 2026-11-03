"""Explicitly authorized outcome-file download; metadata-only export, no melt magnitudes."""
import pathlib,urllib.request,hashlib,csv,io,re,json
URL='https://zenodo.org/api/records/4446387/files/Ablation_stake_data_Lirung_Glacier_CHAND.csv/content'
def main():
 b=urllib.request.urlopen(URL,timeout=30).read();assert len(b)==2026;assert hashlib.md5(b).hexdigest()=='973dd111fb0d340235b3792058cfe339'
 p=pathlib.Path(__file__).parent;(p/'Ablation_stake_data_Lirung_Glacier_CHAND.csv').write_bytes(b)
 try:s=b.decode('utf-8-sig');enc='utf-8-sig'
 except UnicodeDecodeError:s=b.decode('cp1252');enc='cp1252'
 rows=list(csv.reader(io.StringIO(s)));out={'authorization':'Parent relayed explicit user authorization to download Lirung ablation CSV after earlier review rejection. This script exports schema only; numerical melt magnitudes stay masked.','url':URL,'bytes':len(b),'md5':hashlib.md5(b).hexdigest(),'sha256':hashlib.sha256(b).hexdigest(),'encoding':enc,'row_count':len(rows),'schema_masked':[]}
 for row in rows:
  out['schema_masked'].append(['<NUMERIC SEALED>' if re.fullmatch(r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?',x.strip()) else x for x in row])
 (p/'STAKE_SCHEMA_MASKED.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
