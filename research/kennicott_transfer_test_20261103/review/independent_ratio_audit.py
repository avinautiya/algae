"""Independent 2D vertex-enumeration audit, no field data or native targets."""
import sys,itertools,json,hashlib
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'continuous_goal'))
from ratio_enclosure import enclosure

def vertices(A,b,bounds):
 rows=list(np.asarray(A,float));rhs=list(np.asarray(b,float))
 for j,(lo,hi)in enumerate(bounds):
  e=np.zeros(2);e[j]=1;rows.extend([e,-e]);rhs.extend([hi,-lo])
 rows=np.asarray(rows);rhs=np.asarray(rhs);vs=[]
 for i,j in itertools.combinations(range(len(rows)),2):
  M=rows[[i,j]]
  if abs(np.linalg.det(M))<1e-12:continue
  x=np.linalg.solve(M,rhs[[i,j]])
  if np.max(rows@x-rhs)<=1e-9:vs.append(x)
 return vs

def run():
 rng=np.random.default_rng(421);cases=[]
 for i in range(100):
  bounds=[(.1,float(rng.uniform(1,5))),(.1,float(rng.uniform(1,5)))];A=rng.uniform(-1,1,(4,2));center=np.mean(bounds,axis=1);b=A@center+rng.uniform(.1,1,4)
  c=rng.uniform(.1,3,2);f=rng.uniform(0,3,2);Y=100.;vs=vertices(A,b,bounds);expected=[Y*(f@x)/(c@x)for x in vs]
  got=enclosure(A,b,bounds,c,f,Y)
  err=max(abs(got['lower_mmwe']-min(expected)),abs(got['upper_mmwe']-max(expected)))
  assert err<1e-7,(i,err)
  cases.append({'case':i,'vertices':len(vs),'max_endpoint_difference_mmwe':err})
 # Two-dimensional parameter set is broad but all predictions exactly300.
 proportional=enclosure([],[],[(.1,100),(.1,100)],[1,2],[3,6],100,tolerance=1e-7)
 assert proportional['radius_mmwe']<1e-7
 # A narrow TRAIN heat does not establish future heat when independent edge varies.
 nonident=enclosure([[1,0],[-1,0]],[1,-1],[(.1,2),(.1,100)],[1,0],[0,1],100)
 assert abs(nonident['lower_mmwe']-10)<1e-7 and abs(nonident['upper_mmwe']-10000)<1e-7
 source=ROOT/'continuous_goal/ratio_enclosure.py'
 out={'scope':'Independent mathematical implementation audit only; no field-coverage or novelty assertion','source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'seed':421,'100_independent_polygon_vertex_checks':cases,'max_difference_mmwe':max(x['max_endpoint_difference_mmwe']for x in cases),'identified_goal_not_parameters':proportional,'narrow_train_not_future':nonident,'verdict':'PASS for sampled well-scaled finite 2D polytopes; no rigorous outward-rounding or arbitrary dimension numerical guarantee'}
 (Path(__file__).parent/'INDEPENDENT_RATIO_AUDIT.json').write_text(json.dumps(out,indent=2)+'\n');print('100 independent vertex checks PASS; max error',out['max_difference_mmwe'])
if __name__=='__main__':run()
