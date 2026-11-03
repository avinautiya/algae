import datetime as dt,unittest
from arolla_melt import TZ,bins,doses,prediction,uncertainty,native_targets,stamp,published_bins,published_total,load_published,summarize,MODELS,PUBLISHED_NAME

def source(n=96):
 first=dt.datetime(2000,1,1,tzinfo=TZ)
 return dict(forcing=[dict(source_label_local=(first+dt.timedelta(hours=i+1)).isoformat(),features=dict(T=8.,SWIN=100.,SWOUT=50.,SNOW=None)) for i in range(n)])

class Tests(unittest.TestCase):
 def test_available_stake_mean_includes_primary(self):
  rows=[]
  for model in MODELS:
   for stake,error in zip(['H1D','G1D','F1D','E1D'],[10.,30.,40.,80.]):rows.append(dict(model=model,stake=stake,nominal=True,predicted_mmwe=100.+error,heldout_observed_mmwe=100.,absolute_error_mmwe=error))
  rows.append(dict(model=PUBLISHED_NAME,stake='H1D',nominal=True,predicted_mmwe=102.,heldout_observed_mmwe=100.,absolute_error_mmwe=2.))
  result=summarize(rows)
  self.assertEqual(result['constant_rate']['available_stake_nominal_mae_mmwe'],40.)
  self.assertEqual(result['constant_rate']['nominal_stakes'],4)
  self.assertEqual(result[PUBLISHED_NAME]['available_stake_nominal_mae_mmwe'],2.)
  self.assertEqual(result[PUBLISHED_NAME]['nominal_stakes'],1)
 def test_published_units_and_partial_steps(self):
  start=dt.datetime(2000,1,1,tzinfo=TZ);r=[dict(label_local=(start+dt.timedelta(hours=i+1)).isoformat(),model_melt_mwe_per_step=.001) for i in range(3)]
  s=published_bins(r,0.);self.assertAlmostEqual(published_total(s,start+dt.timedelta(minutes=30),start+dt.timedelta(hours=2)),1.5)
 def test_published_clock_scenarios_explicit(self):
  start=dt.datetime(2000,1,1,tzinfo=TZ);r=[dict(label_local=(start+dt.timedelta(hours=i+1)).isoformat(),model_melt_mwe_per_step=.001) for i in range(4)]
  for off in [0.,.5,1.]:self.assertAlmostEqual(published_total(published_bins(r,off),start+dt.timedelta(hours=1),start+dt.timedelta(hours=3)),2.)
 def test_actual_published_members_verify_without_target_score(self):
  r=load_published();self.assertEqual(len(r),1686);self.assertEqual(r[0]['label_local'],'2010-07-04T12:00:00+02:00')
 def test_exact_model_doses(self):
  d=source();a=stamp('2000-01-01',0);b=stamp('2000-01-02',0);x=doses(bins(d,0.),a,b)
  self.assertEqual(x['constant_rate'],1.);self.assertAlmostEqual(x['positive_degree_day'],8.);self.assertAlmostEqual(x['net_shortwave_dose'],4.32)
 def test_unknown_snow_is_not_replaced(self):
  d=source();self.assertTrue(all(r['features']['SNOW'] is None for r in d['forcing']));doses(bins(d,0.),stamp('2000-01-01',0),stamp('2000-01-02',0));self.assertTrue(all(r['features']['SNOW'] is None for r in d['forcing']))
 def test_missing_forcing_refuses(self):
  d=source();d['forcing'][2]['features']['T']=None
  with self.assertRaises(ValueError):doses(bins(d,0.),stamp('2000-01-01',0),stamp('2000-01-02',0))
 def test_forcing_gap_or_clock_refuses(self):
  d=source();del d['forcing'][2]
  with self.assertRaises(ValueError):bins(d,0.)
  with self.assertRaises(ValueError):bins(source(),.25)
 def test_no_extrapolation(self):
  with self.assertRaises(ValueError):doses(bins(source(),0.),stamp('1999-12-31',0),stamp('2000-01-02',0))
 def test_training_only_scale(self):
  pred,f,coef=prediction(100.,20.,10.);self.assertEqual((pred,f,coef),(50.,.5,5.))
  # There is intentionally no test-target argument in the prediction API.
  with self.assertRaises(ValueError):prediction(100.,0.,10.)
 def test_native_cumulative_measurement_operator(self):
  s=dict(points=[dict(source_melt_mm_water_equivalent=v) for v in [0,100,160]],intervals=[dict(observed_ice_ablation_mm_water_equivalent=v) for v in [100,60]])
  self.assertEqual(native_targets(s),(100,60));s['intervals'][1]['observed_ice_ablation_mm_water_equivalent']=999
  with self.assertRaises(ValueError):native_targets(s)
 def test_shared_middle_reading_bounds(self):
  r=uncertainty(1.,.5,.1,900.);self.assertEqual(r['if_per_reading_absolute_bound_residual_bound_mmwe'],270.);self.assertEqual(r['if_per_interval_absolute_bound_residual_bound_mmwe'],135.);self.assertFalse(r['confidence_interval'])
 def test_shared_clock_ratio(self):
  s=bins(source(),0.);a=stamp('2000-01-01',12);b=stamp('2000-01-02',12);c=stamp('2000-01-03',12)
  x=doses(s,a,b)['positive_degree_day'];y=doses(s,b,c)['positive_degree_day'];self.assertEqual(y/x,1.)
  later=b+dt.timedelta(hours=6);xx=doses(s,a,later)['positive_degree_day'];yy=doses(s,later,c)['positive_degree_day'];self.assertLess(yy/xx,y/x)
if __name__=='__main__':unittest.main()
