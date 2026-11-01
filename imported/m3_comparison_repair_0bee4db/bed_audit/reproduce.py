"""Reproduce pinned lower_bed defects without running hydraulics or a campaign."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
from .operator import lower_ice

ROOT=Path(__file__).resolve().parents[1]


def pinned_module():
    import importlib.util
    path=ROOT/'snapshot'/'swe2d.py'
    spec=importlib.util.spec_from_file_location('_m3_pinned_swe2d',path)
    module=importlib.util.module_from_spec(spec)
    sys.modules[spec.name]=module;spec.loader.exec_module(module)
    return module


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'bed_audit'/'source_reproduction.json')
    args=parser.parse_args()
    mod=pinned_module();results={}
    for dx in (.1,.05,.025):
        s=mod.Sim(dx=dx);before=s.zb.copy();s.last_qacc=np.zeros_like(before)
        declared=s.lower_bed(0.)[0]
        changed=(s.zb!=before)&(~s.ice)
        zero=dict(changed_nonice_cells=int(changed.sum()),
                  maximum_nonice_bed_raise_m=float(np.max(np.where(changed,s.zb-before,0.))),
                  changed_interior_nonice_volume_m3=float(np.where(changed,s.zb-before,0.)[1:-1,1:-1].sum()*dx**2),
                  declared_melt_m3=declared)
        # Isolated uniform accepted heat flux, NOT a campaign field/result.
        s=mod.Sim(dx=dx);before=s.zb.copy();s.last_qacc=np.where(s.ice,100.,0.)
        declared=s.lower_bed(21600.)[0]
        actual=float(np.where(s.ice,before-s.zb,0.).sum()*dx**2)
        unconstrained=lower_ice(before,s.ice,s.last_qacc,21600.,dx**2)
        results[str(dx)]=dict(zero_time_update=zero,
            illustrative_uniform_100W_m2_6h=dict(declared_melt_m3=declared,actual_ice_lowering_volume_m3=actual,
                declared_minus_actual_m3=declared-actual,
                capped_ice_cells=int(np.sum(s.ice&(s.zb==s.pool_floor))),
                corrected_operator_volume_m3=unconstrained['ice_volume_m3'],
                corrected_nonice_unchanged=unconstrained['nonice_unchanged']))
    record=dict(pinned_commit='f3602e35e1595bcdbee4208a354a34944f1ae974',
                source_sha256=hashlib.sha256((ROOT/'snapshot'/'swe2d.py').read_bytes()).hexdigest(),
                evidence='direct pinned-source operator reproduction; uniform heat fixture is synthetic, not measured campaign data',
                results=results,
                consequences=['Zero-time update changes non-ice geometry.',
                              'Counted latent volume can exceed actual ice lowering where floor clipping activates.',
                              'Actual campaign impact requires original field/state exports or a versioned rerun; this fixture does not quantify it.'])
    args.output.write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
    print(json.dumps(record,indent=2,allow_nan=False))


if __name__=='__main__':main()
