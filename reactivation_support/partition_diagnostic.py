"""Read-only diagnostic on a pinned copy; does not modify Claude's solver."""
from pathlib import Path
import argparse,json,time
import numpy as np
from .snapshot import cases
from .snapshot.swe2d import G,RHO_W,CP,HDRY,_hll


def face_discharge(s,i,ys):
    zs=s.zs();Q=0.
    for j in range(ys.start,ys.stop):
        a=i-1;b=i;zf=max(zs[a,j],zs[b,j])
        hl=max(0.,s.h[a,j]+zs[a,j]-zf);hr=max(0.,s.h[b,j]+zs[b,j]-zf)
        ul=s.hu[a,j]/s.h[a,j] if s.h[a,j]>HDRY else 0.
        vl=s.hv[a,j]/s.h[a,j] if s.h[a,j]>HDRY else 0.
        ur=s.hu[b,j]/s.h[b,j] if s.h[b,j]>HDRY else 0.
        vr=s.hv[b,j]/s.h[b,j] if s.h[b,j]>HDRY else 0.
        Q+=_hll(hl,ul,vl,hr,ur,vr,G)[0]*s.dx
    return float(Q)

def section(s,i,ys):
    h=s.h[i,ys];hu=s.hu[i,ys];hv=s.hv[i,ys];z=s.zs()[i,ys]
    u=np.divide(hu,h,out=np.zeros_like(h),where=h>HDRY)
    v=np.divide(hv,h,out=np.zeros_like(h),where=h>HDRY)
    head=z+h+(u*u+v*v)/(2*G)
    q=hu*s.dx;positive=np.maximum(q,0);total=positive.sum()
    wet=h>HDRY
    return dict(i=int(i),Q=float(q.sum()),numerical_face_Q=face_discharge(s,i,ys),stage_min=float(np.min((h+z)[wet])) if wet.any() else None,
                stage_max=float(np.max((h+z)[wet])) if wet.any() else None,
                positive_flux_mean_total_head=float(positive@head/total) if total else None,
                froude_max=float(np.max(np.divide(abs(u),np.sqrt(G*h),out=np.zeros_like(h),where=wet))),
                reverse_volume_flux=float(-np.minimum(q,0).sum()))

def one(layout,dx=.05,manning=.02,**geometry):
    start=time.perf_counter();s=cases.build(layout,dx=dx,manning=manning,**geometry)
    convergence=s.hydraulics(tmax=500.,cool=False)
    d=s.diagnostics();allcols=slice(1,s.ny+1)
    pool=section(s,max(1,s.i0-3),allcols)
    profiles=[]
    for i in range(max(1,s.i0-3),min(s.nx,s.i1+3)+1):
        profiles.append(dict(toe=section(s,i,s.toe_cols()),bypass=section(s,i,s.byp_cols())))
    critical=[]
    for label,ys in [('toe',s.toe_cols()),('bypass',s.byp_cols())]:
        q=s.Qsec(ys);w=s.nw*dx;crest=cases.crest(s,label if label=='toe' else 'byp')
        required=crest+1.5*((q/w)**2/G)**(1/3)
        critical.append(dict(route=label,measured_Q=q,geometric_crest=crest,
                             ideal_critical_energy_head=required,
                             pool_minus_ideal=pool['positive_flux_mean_total_head']-required))
    return dict(layout=layout,dx=dx,manning=manning,geometry=geometry,convergence=convergence,diagnostics=d,
                pool=pool,critical=critical,profiles=profiles,midsection_numerical_face_toe_fraction=face_discharge(s,s.i0+(s.i1-s.i0)//2,s.toe_cols())/s.Q,wall_s=time.perf_counter()-start,
                warning='Cell-centred total-head diagnostic, not an exact discrete mechanical-energy budget; reference applicability remains to be checked.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--dx',type=float,default=.05);args=p.parse_args()
    rows=[]
    for layout,n in [('bare',.02),('strip',.02),('bare',0.),('strip',0.)]:
        r=one(layout,args.dx,n);rows.append(r)
        print(layout,n,r['diagnostics']['f_toe'],r['convergence'],r['critical'],flush=True)
        args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(rows,indent=2,allow_nan=False))
