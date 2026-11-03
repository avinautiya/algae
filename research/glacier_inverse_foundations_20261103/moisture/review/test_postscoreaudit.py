import datetime as dt,math,unittest
from postscoreaudit import positive_linear_integral,native_air_area,reconstruct,FIXED_PARAMS
class IndependentDoseTests(unittest.TestCase):
 def test_positive_linear_crossings(self):
  self.assertEqual(positive_linear_integral(-2,2,10),5)
  self.assertEqual(positive_linear_integral(2,-2,10),5)
  self.assertEqual(positive_linear_integral(-2,-1,10),0)
  self.assertEqual(positive_linear_integral(2,4,10),30)
 def test_native_air_split_independent(self):
  o=dt.datetime(2022,1,1);air=[(o, -2),(o+dt.timedelta(seconds=3600),2)]
  whole=native_air_area(o,o+dt.timedelta(seconds=3600),air)
  split=sum(native_air_area(o+dt.timedelta(seconds=j),o+dt.timedelta(seconds=j+600),air)for j in range(0,3600,600))
  self.assertAlmostEqual(whole,1800);self.assertAlmostEqual(whole,split)
 def test_independent_dose_constant_and_offset(self):
  o=dt.datetime(2022,1,1);labels=[o+dt.timedelta(seconds=1800*i)for i in range(11)]
  features=[{'SOURCE_LOCAL_LABEL':t.isoformat(),'SM3_raw_m3_m3':str(.1+.001*i)}for i,t in enumerate(labels)]
  air=[(t,2.)for t in labels];temps={t.isoformat():{'DT1':'3','DT4':'2','DT5':'.9'}for t in labels}
  epochs=[labels[i].isoformat()for i in [0,2,4,6,8,10]]
  a=reconstruct(features,air,temps,epochs,FIXED_PARAMS);b=reconstruct(features,air,temps,epochs,FIXED_PARAMS,offset=.2)
  self.assertEqual(a['air_pdd'],[7200.]*5);self.assertEqual(a['constant_rate'],[3600.]*5)
  for m in ['anomaly_exp_air','centered_softplus_air','moisture_only_exp']:
   for i in range(1,5):self.assertAlmostEqual(a[m][i]/a[m][0],b[m][i]/b[m][0],places=13)
 def test_no_extrapolation_or_long_gap(self):
  o=dt.datetime(2022,1,1)
  with self.assertRaises(ValueError):native_air_area(o-dt.timedelta(seconds=1),o,[(o,1),(o+dt.timedelta(seconds=3600),1)])
  with self.assertRaises(ValueError):native_air_area(o,o+dt.timedelta(seconds=3600),[(o,1),(o+dt.timedelta(seconds=70000),1)])
if __name__=='__main__':unittest.main()

class RecalibratedSharedReadings(unittest.TestCase):
 def test_all_shared_six_endpoint_corners(self):
  import itertools
  from postscoreaudit import recalibrated_pairwise_lipschitz_per_mmice
  c=[.5,1,2,1.5];b=[1,1,1,1];truth=[10,11,12,14];train=10.;delta=.2
  nominal=sum(abs(train*x-y)-abs(train*z-y)for x,z,y in zip(c,b,truth))/4
  lip=recalibrated_pairwise_lipschitz_per_mmice(c,b)
  for e in itertools.product((-delta,delta),repeat=6):
   yt=train+.917*(e[0]-e[1]);y=[v+.917*(e[j+1]-e[j+2])for j,v in enumerate(truth)]
   loss=sum(abs(yt*x-v)-abs(yt*z-v)for x,z,v in zip(c,b,y))/4
   self.assertLessEqual(abs(loss-nominal),lip*delta+1e-12)
