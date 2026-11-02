"""Reproduce bounded v5 checks; does not launch a production campaign."""
import hashlib,io,json,sys,time,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from swe2d_v5 import Sim
from test_swe2d_v5 import flat_cell_domain

def main():
 stream=io.StringIO();suite=unittest.defaultTestLoader.discover(str(ROOT),pattern='test_swe2d_v5.py')
 result=unittest.TextTestRunner(stream=stream,verbosity=2).run(suite)
 c=flat_cell_domain();before=dict(volume_m3=c.liquid_volume(),sensible_J=c.sensible_heat());closed=c.step_physical(2.)
 o=Sim(dx=.25,Lpool=.5,Lr=.5,Ltail=.25,Wl=.25,Q=.01,Tw=.5,ramp_ratio=0)
 for _ in range(10):o.step_physical(.2)
 sha=hashlib.sha256((ROOT/'swe2d_v5.py').read_bytes()).hexdigest()
 record=dict(tests_run=result.testsRun,passed=result.wasSuccessful(),test_log=stream.getvalue(),
   closed_cell_before=before,closed_cell_step=closed,open_boundary_ledgers=o.combined_ledgers(),
   python=sys.version,source_sha256=sha,
   scope='actual inheritedv4hydraulic kernel plus phase source; no field validation or contact-event campaign')
 (ROOT/'VERIFICATION.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
 if not result.wasSuccessful():print(stream.getvalue());raise SystemExit(1)
 rows=[]
 for dt in [1.,.5,.25,.125]:
  s=Sim(dx=.25,Lpool=.5,Lr=.5,Ltail=.25,Wl=.25,Q=.01,Tw=2.,heff=5000.,ramp_ratio=0)
  start=time.perf_counter()
  for _ in range(round(20/dt)):s.step_physical(dt)
  led=s.combined_ledgers();rows.append(dict(dt_s=dt,duration_s=20.,melt_m3=led['ice_volume_lost_m3'],
    toe_discharge_m3_s=s.Qsec(s.toe_cols()),water_volume_m3=s.liquid_volume(),
    mass_residual_kg=led['combined_ice_water_mass_residual_kg'],
    thermal_residual_J=led['physical_split_thermal_latent_residual_J'],runtime_s=time.perf_counter()-start))
 fine=rows[-1]
 for row in rows:row['relative_melt_difference_to_finest']=row['melt_m3']/fine['melt_m3']-1
 study=dict(scope='short splitting sensitivity on suppliedidealised input; not field or contact-event convergence',
   source_sha256=sha,case=dict(dx_m=.25,Q_m3_s=.01,Tw_K=2.,heff_W_m2_K=5000.,duration_s=20.),rows=rows)
 (ROOT/'SHORT_TIME_REFINEMENT.json').write_text(json.dumps(study,indent=2,allow_nan=False)+'\n')
 print(json.dumps(dict(passed=True,tests=result.testsRun,source_sha256=sha),indent=2))
if __name__=='__main__':main()
