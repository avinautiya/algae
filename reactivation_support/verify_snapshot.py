"""Reproduce the pinned source defect and benchmark the independent reference."""
from pathlib import Path
import json,time,statistics
import numpy as np
from .snapshot.swe2d import Sim,RHO_I,RHO_W,LF,CP
from .snapshot.verify_swe import gvf_ref
from .hydraulic_reference import partition,Lane
from .integration_adapter import conservative_lower_bed


def mass_fixture():
    rows=[]
    for fixed in (False,True):
        s=Sim(dx=.05,Q=0.,closed_out=True)
        s.last_qacc=np.where(s.ice,100.,0.)
        V0=float(s.h.sum()*s.dx**2);E0=float(RHO_W*CP*s.hT.sum()*s.dx**2)
        if fixed:vol,_,_=conservative_lower_bed(s,3600.)
        else:vol,_,_=s.lower_bed(3600.)
        observed=float(s.h.sum()*s.dx**2)-V0;expected=vol*RHO_I/RHO_W
        rows.append(dict(adapter=fixed,ice_melt_volume_m3=vol,expected_water_added_m3=expected,
                         actual_water_added_m3=observed,water_source_error_m3=observed-expected,
                         sensible_energy_change_J=float(RHO_W*CP*s.hT.sum()*s.dx**2)-E0,
                         heat_contract='latent supplied externally; geometry update must not remove sensible heat twice'))
    return rows


def benchmark_reference():
    lanes=[(-.08,.5),(0.,.5)]
    new=[Lane(z,w) for z,w in lanes]
    # Warm imports and numerical code before serial timing.
    old=gvf_ref(.07,lanes);latest=partition(.07,new)
    old_times=[];new_times=[]
    for _ in range(5):
        start=time.perf_counter();gvf_ref(.07,lanes);old_times.append(time.perf_counter()-start)
        start=time.perf_counter();partition(.07,new);new_times.append(time.perf_counter()-start)
    return dict(old_stage=old[0],new_stage=latest['stage'],
                old_toe_fraction=float(old[1][0]),new_toe_fraction=latest['fractions'][0],
                old_median_s=statistics.median(old_times),new_median_s=statistics.median(new_times),
                speedup=statistics.median(old_times)/statistics.median(new_times),
                limitation='same prismatic reference assumptions; not applicable as an exact strip-layout truth')

if __name__=='__main__':
    result=dict(mass_fixture=mass_fixture(),reference_benchmark=benchmark_reference())
    p=Path(__file__).parent/'records/snapshot_reproduction.json';p.parent.mkdir(exist_ok=True)
    p.write_text(json.dumps(result,indent=2,allow_nan=False));print(json.dumps(result,indent=2))
