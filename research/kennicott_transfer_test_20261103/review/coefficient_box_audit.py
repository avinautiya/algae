"""Independent fixed-polytope coefficient-box check, no native data."""
import sys,itertools,json
from pathlib import Path
import numpy as np
from independent_ratio_audit import vertices,ROOT
from ratio_enclosure import coefficient_box_enclosure
rng=np.random.default_rng(5021);checks=[]
for k in range(30):
 B=[(.2,2.),(.1,3.)];A=rng.uniform(-1,1,(3,2));b=A@np.array([1.1,1.55])+rng.uniform(.2,1.,3);vs=vertices(A,b,B)
 cl=rng.uniform(.2,1.,2);cu=cl+rng.uniform(.01,.3,2);fl=rng.uniform(.2,1.,2);fu=fl+rng.uniform(.01,.3,2)
 brute=[]
 for bits in itertools.product((0,1),repeat=4):
  c=np.array([cl[i]if bits[i]==0 else cu[i]for i in range(2)]);f=np.array([fl[i]if bits[i+2]==0 else fu[i]for i in range(2)])
  brute.extend(100*f@v/(c@v)for v in vs)
 r=coefficient_box_enclosure(A,b,B,cl,cu,fl,fu,100)
 e=max(abs(min(brute)-r['lower_mmwe']),abs(max(brute)-r['upper_mmwe']));assert e<1e-7
 checks.append({'case':k,'vertex_count':len(vs),'independent_all16_coefficient_corners_max_difference_mmwe':e})
p=Path(__file__).parent/'COEFFICIENT_BOX_AUDIT.json';p.write_text(json.dumps({'scope':'30 random positive-coefficient nonnegative-parameter polygon checks; no nonlinear-phase or true-sensor-error model guarantee','checks':checks},indent=2)+'\n');print('30 coefficient-box independent all-corner audits PASS; max',max(x['independent_all16_coefficient_corners_max_difference_mmwe']for x in checks))
