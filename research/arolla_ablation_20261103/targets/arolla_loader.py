"""Verified Arolla forcing and source-provided stake ice-ablation targets.

The loader authorizes comparison design, not accuracy or a new model claim.
Visit hours remain unknown. SNOW is never silently zero-filled.
"""
import argparse,csv,datetime as dt,hashlib,io,json,math,xml.etree.ElementTree as ET,zipfile
from pathlib import Path
NS={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
LOCAL=dt.timezone(dt.timedelta(hours=2))
def scalar(s):
 if str(s).strip().lower() in {'','nan','na','n/a'}:return None
 v=float(s)
 if not math.isfinite(v):raise ValueError('nonfinite numeric source')
 return v
def source_time(r):
 h=int(r['HOUR']);m=int(r['MINUTE'])
 if h==24 and m!=0:raise ValueError('24hour source label requires minute0')
 return (dt.datetime(int(r['YEAR']),int(r['MONTH']),int(r['DAY']),0 if h==24 else h,m,tzinfo=LOCAL)+(dt.timedelta(days=1) if h==24 else dt.timedelta())).isoformat()
def workbook(raw):
 with zipfile.ZipFile(io.BytesIO(raw)) as z:
  ss=[ ''.join(q.itertext()) for q in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('s:si',NS)]
  cells={}
  for c in ET.fromstring(z.read('xl/worksheets/sheet1.xml')).findall('.//s:c',NS):
   v=c.find('s:v',NS)
   if v is not None:cells[c.attrib['r']]=ss[int(v.text)] if c.attrib.get('t')=='s' else v.text
 return cells
def load(archive,record):
 raw=Path(archive).read_bytes();meta=json.loads(Path(record).read_text());f=next(x for x in meta['files'] if x['key']=='1_ARO_2010.zip')
 if hashlib.md5(raw).hexdigest()!=f['checksum'].split(':')[1]:raise ValueError('archive checksum mismatch')
 with zipfile.ZipFile(io.BytesIO(raw)) as z:
  fn='ARO_2010/ARO_2010_60.csv';sn='ARO_2010/ARO_2010_VAL_stakes.xlsx';mn='ARO_2010/ARO_2010_AWS_metadata.docx'
  rows=list(csv.DictReader(io.StringIO(z.read(fn).decode())));forcing=[];last=None
  for r in rows:
   t=source_time(r)
   if last is not None and t<=last:raise ValueError('nonmonotone or duplicate forcing')
   last=t;forcing.append(dict(source_label_local=t,features={k:scalar(v) for k,v in r.items() if k not in ['YEAR','MONTH','DAY','HOUR','MINUTE']}))
  c=workbook(z.read(sn));rho=float(c['J4'])
  if rho!=900 or c['C16']!='H1D':raise ValueError('unexpected workbook schema/density')
  dates={col:(dt.date(1899,12,30)+dt.timedelta(days=int(float(c[col+'2'])))).isoformat() for col in ['D','E','F','G']}
  stakes=[]
  for row,rawrow in [(16,4),(17,5),(18,6),(19,7)]:
   sid=c['C'+str(row)];points=[]
   for col in ['D','E','F','G']:
    if col+str(row) not in c:continue
    mm=float(c[col+str(row)]);cm=float(c[col+str(rawrow)])
    points.append(dict(date_local=dates[col],date_precision='day',source_melt_mm_water_equivalent=mm,original_stake_reading_cm=cm,visit_hour_unknown=True))
   intervals=[]
   for a,b in zip(points,points[1:]):intervals.append(dict(start_date=a['date_local'],end_date=b['date_local'],observed_ice_ablation_mm_water_equivalent=b['source_melt_mm_water_equivalent']-a['source_melt_mm_water_equivalent'],visit_time_range='each visit anywhere in its local calendar day; propagate jointly',source_unit='source workbook MELT [mm w.e.], not a heat-flux conversion'))
   stakes.append(dict(stake=sid,closest_to_AWS=sid=='H1D',source_coordinates=dict(X=float(c['A'+str(rawrow)]),Y=float(c['B'+str(rawrow)]),crs_unspecified=True),points=points,intervals=intervals))
  hashes={n:hashlib.sha256(z.read(n)).hexdigest() for n in [fn,sn,mn,'ARO_2010/ARO_2010_debris.csv']}
 return dict(site='Haut Glacier d Arolla 2010',source_doi='10.5281/zenodo.3047649',archive_sha256=hashlib.sha256(raw).hexdigest(),archive_md5=hashlib.md5(raw).hexdigest(),member_sha256=hashes,clock='local UTC+02:00 explicitly declared by source',forcing=forcing,stake_targets=stakes,source_ice_density_kg_m3=rho,source_stake_uncertainty_m=.1,uncertainty_semantics='source says uncertainty0.1m; probability/sigma/per-reading association unspecified, do not invent confidence intervals',debris=dict(primary_H1D_thickness_m=.06,thickness_relative_uncertainty=.1,other_stake_thickness_not_provided=True,conductivity_W_m_K=.96,conductivity_provenance='supplied value; its DOCX row has no provenance note, so independent local measurement or a particular borrowed source is not established',roughness_provenance='roughness z0 row explicitly borrows Miage properties; do not silently assign that note to conductivity'),snow=dict(post_July4_native_depth_unknown=True,metadata_events=['2010-07-24','2010-07-30','2010-08-06','2010-08-15','2010-08-30'],August15_cover_duration_days=[3,4],other_events_melt_within_hours=24,silent_missing_to_zero_forbidden=True),scope='stake-inferred sub-debris ice ablation with bounded visit timing; NOT coupled cliff/lake/stream or glacier-scale validation',validation_outcome_not_computed=True)
def main():
 p=argparse.ArgumentParser();p.add_argument('archive',type=Path);p.add_argument('record',type=Path);p.add_argument('out',type=Path);a=p.parse_args()
 if a.out.exists():raise ValueError('new output required')
 d=load(a.archive,a.record);a.out.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n');print(json.dumps(dict(forcing_rows=len(d['forcing']),stakes=len(d['stake_targets']),ice_intervals=sum(len(s['intervals']) for s in d['stake_targets']),validation_outcome_not_computed=True)))
if __name__=='__main__':main()
