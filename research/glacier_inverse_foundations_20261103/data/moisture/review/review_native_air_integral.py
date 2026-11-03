"""Independent native-air integral review; no camera or melt targets."""
import sys,csv,datetime as dt,json,math,hashlib
from pathlib import Path
import numpy as np
ROOT=Path(__file__).parents[3];sys.path.insert(0,str(ROOT/'moisture'))
from anomaly_index import load_site,dose,air_exposure

def independent_integral(source,start,end):
    with Path(source).open() as file:rows=list(csv.DictReader(file))
    times=np.array([(dt.datetime.fromisoformat(r['SOURCE_LOCAL_LABEL'])-start).total_seconds() for r in rows]);values=np.array([float(r['TA_3']) for r in rows])
    last=(end-start).total_seconds();knots=np.r_[0.,times[(times>0)&(times<last)],last]
    temps=np.interp(knots,times,values);roots=[]
    for a,b,x,y in zip(knots[:-1],knots[1:],temps[:-1],temps[1:]):
        if x*y<0:roots.append(a-x*(b-a)/(y-x))
    split=np.sort(np.r_[knots,roots]);positive=np.maximum(np.interp(split,times,values),0.)
    return math.fsum((positive[:-1]+positive[1:])*.5*np.diff(split))

reports=[]
for tower in ['Tower1','Tower2']:
 f=load_site(tower);source=ROOT/'data/moisture'/f'{tower}_NATIVE_AIR_TA3.csv';direct=[]
 for a,b,s in zip(f['epochs'][:-1],f['epochs'][1:],f['slices']):
    independent=independent_integral(source,dt.datetime.fromisoformat(a),dt.datetime.fromisoformat(b));actual=math.fsum(f['air'][s])
    assert math.isclose(actual,independent,rel_tol=2e-14,abs_tol=1e-7)
    direct.append({'start':a,'end':b,'production_K_seconds':actual,'independent_K_seconds':independent,'absolute_difference_K_seconds':abs(actual-independent)})
 plain=dose(f,'air_pdd',0.)
 for method in ['anomaly_exp_air','shifted_ratio_air','calendar_exp_air']:np.testing.assert_allclose(dose(f,method,0.),plain,rtol=0,atol=0)
 np.testing.assert_allclose(dose(f,'centered_softplus_air',0.),plain*np.log(2),rtol=2e-15)
 np.testing.assert_allclose(dose(f,'moisture_only_exp',0.),dose(f,'constant_rate'),rtol=0,atol=0)
 reports.append({'tower':tower,'max_native_air_gap_hours':f['source_max_air_gap_seconds']/3600,'integral_comparisons':direct,'zero_parameter_closures_pass':True})
out={'scope':'Nativeairintegration/calculus only; no camera outcomes','source_anomaly_index_sha256':hashlib.sha256((ROOT/'moisture/anomaly_index.py').read_bytes()).hexdigest(),'review':reports,'independent_operator':'Splitnativepiecewiselinearairatendpointsandallzero-crossings, trapezoidmax(T,0)exactoneachpiece','no_moisture_interpolation_claim':'ProductionholdsSM3fixedatleftnativehalfhourvalue for each exactair-exposurecell; this explicitquadratureassumption is notcontinuous measuredmoisture','passed':True}
(Path(__file__).parent/'NATIVE_AIR_INTEGRAL_REVIEW.json').write_text(json.dumps(out,indent=2)+'\n')
print('Bothsites all10periodnativepositiveairintegrals independentlyverified; zero-parameterclosurespass.')
