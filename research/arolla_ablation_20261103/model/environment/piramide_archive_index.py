import struct,json
from pathlib import Path
base=Path(__file__).parent
b=(base/'piramide_zip_tail.bin').read_bytes(); i=b.rfind(b'PK\x05\x06')
e=struct.unpack_from('<4s4H2LH',b,i); n=e[4]; size=e[5]; off=e[6]
p=off-(62677571-len(b)); rows=[]
for _ in range(n):
 x=struct.unpack_from('<4s6H3L5H2L',b,p)
 name=b[p+46:p+46+x[10]].decode('utf-8')
 rows.append({'name':name,'method':x[4],'crc32':x[7],'compressed_size':x[8],'size':x[9],'offset':x[16]})
 p+=46+x[10]+x[11]+x[12]
(base/'PIRAMIDE_ARCHIVE_INDEX.json').write_text(json.dumps(rows,indent=2)+'\n')
for r in rows:
 if r['name'].endswith(('Readme.txt','DATA_Piramide.mat','LithologyData.mat','SM_data.mat','Ablation_Tower1.mat','RunThermalConductivity_allSites.m','Calibrate_deb_par_PIR_k_z0.m')): print(r)
