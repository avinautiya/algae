"""Verification-only geometry isolation; NOT a replacement production layout.

Changing the inlet/bump changes the physical problem. Results cannot justify
changing a production geometry merely to pass its old reference tolerance.
"""
from pathlib import Path
import json,math,time
import numpy as np
from .snapshot import cases
from .snapshot.verify_swe import weir_ref
from .partition_diagnostic import section,face_discharge


def one(layout,*,smooth_inlet,smooth_deposit,dx=.05):
    s=cases.build(layout,dx=dx,manning=0.,Lpool=3.,tail_drop=0.)
    if smooth_inlet:
        for i in range(1,s.i0):
            x=(i-.5)*dx;f=min(1.,max(0.,(x-.5)/2.5))
            s.zb[i,s.toe_cols()]=-1.5+(-.08+1.5)*f
            s.zb[i,s.byp_cols()]=-1.5+1.5*f
            s.zb[i,s.nw+1]=-1.5+1.5*f
    if smooth_deposit and layout=='strip':
        # Smooth transcritical bump: peak and area differ from rectangular strip.
        s.solid[:]=False;s.zd[:]=0
        for i in range(s.i0,s.i1+1):
            x=(i-s.i0+.5)*dx
            s.solid[i,s.toe_cols()]=True
            s.zd[i,s.toe_cols()]=.2*math.sin(math.pi*x/s.Lr)**2
    s.reset_state(.15)
    start=time.perf_counter();r=s.hydraulics(tmax=1200.,cool=False)
    d=s.diagnostics();eta,f=weir_ref(.07,[cases.crest(s,'toe'),cases.crest(s,'byp')],[.5,.5])
    mid=s.i0+(s.i1-s.i0)//2;pool=section(s,1,slice(1,s.ny+1))
    return dict(layout=layout,smooth_inlet=smooth_inlet,smooth_deposit=smooth_deposit,dx=dx,
                numerical_face_fraction=face_discharge(s,mid,s.toe_cols())/.07,
                cell_fraction=d['f_toe'],ideal_reference_fraction=float(f[0]),
                pool=pool,reference_stage=float(eta),convergence=r,wall_s=time.perf_counter()-start,
                note='different verification geometries; not independent physical validation of abrupt-step production case')

if __name__=='__main__':
    out=[]
    for args in [('bare',False,False),('bare',True,False),('strip',True,False),('strip',True,True)]:
        row=one(args[0],smooth_inlet=args[1],smooth_deposit=args[2]);out.append(row);print(row,flush=True)
        (Path(__file__).parent/'records/smooth_control_diagnostic.json').write_text(json.dumps(out,indent=2))
