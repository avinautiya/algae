import csv,datetime as dt,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
from run_benchmark import execute_one

def fixture():
 t=np.arange(8065)*300.;w=2*np.pi*t/86400
 top=8+3*np.sin(w);bottom=.2+.1*np.sin(w-.7);inside=3+1.2*np.sin(w-.8)
 v=np.column_stack([top,inside,bottom]);z=np.array([.02,.1,.2]);origin=dt.datetime(2020,1,1,tzinfo=dt.timezone.utc)
 dates=[origin+dt.timedelta(seconds=float(s)) for s in t]
 return t,dates,z,v
class LeakageTests(unittest.TestCase):
 def test_test_interior_does_not_change_parameters_or_predictions(self):
  original=fixture();changed=tuple(original[:3])+(original[3].copy(),);changed[3][4032:,1]+=100
  split=dict(profile_id='SYNTHETIC_LEAK_TEST',source_sha256='0'*64)
  with tempfile.TemporaryDirectory() as tmp:
   a=Path(tmp)/'a';b=Path(tmp)/'b'
   with patch('run_benchmark.load_profile',return_value=original):ra=execute_one(split,a)
   with patch('run_benchmark.load_profile',return_value=changed):rb=execute_one(split,b)
   for k in ['alpha_m2_s','alpha_grid_training_MSE_C2','filter_parameters','ridge_penalty','ridge_inner_validation_MSE_C2']:
    self.assertEqual(ra[k],rb[k])
   with (a/ra['prediction_csv']).open() as f:aa=list(csv.DictReader(f))
   with (b/rb['prediction_csv']).open() as f:bb=list(csv.DictReader(f))
   for x,y in zip(aa,bb):
    for model in ['diffusion','quasisteady','lag_filter','ridge']:self.assertEqual(x[model+'_C'],y[model+'_C'])
 def test_future_boundary_does_not_change_fitted_parameters(self):
  original=fixture();changed=tuple(original[:3])+(original[3].copy(),);changed[3][4032:,[0,2]]+=10
  split=dict(profile_id='SYNTHETIC_BOUNDARY_TEST',source_sha256='0'*64)
  with tempfile.TemporaryDirectory() as tmp:
   with patch('run_benchmark.load_profile',return_value=original):ra=execute_one(split,Path(tmp)/'a')
   with patch('run_benchmark.load_profile',return_value=changed):rb=execute_one(split,Path(tmp)/'b')
   for k in ['alpha_m2_s','alpha_grid_training_MSE_C2','filter_parameters','ridge_penalty','ridge_inner_validation_MSE_C2']:self.assertEqual(ra[k],rb[k])
if __name__=='__main__':unittest.main()
