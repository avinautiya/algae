from __future__ import annotations
import csv,json,sys,hashlib,os
from pathlib import Path
from datetime import datetime,timedelta,timezone
import numpy as np
from scipy.signal import lfilter
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'solver'))
from dirichlet_slab import PreparedSlab
LF=334000.
PROVENANCE={'top_boundary':'measured_probe','bottom_boundary':'assumed_ice_Tm0','source':'frozen ablation-transfer thickness scenario; basal0 is phase assumption not sensor'}
MODES=int(os.environ.get('ABLATION_MODES','64'))
IDS=['2020b','2023a','2023b','2023c']; MASS=[360,1035,873,688]
PROFILE=[.22,.24,.32,.26]; STAKE=[.21,.24,.29,.255]
def stamp(s):return datetime.strptime(s,'%Y-%m-%d %H:%M:%S').replace(tzinfo=timezone.utc).timestamp()
def read(i):
 p=ROOT/'source/Debris_Temperature_Profiles/Profile_Temperature_Records'/('KDT_'+i+'.csv')
 with p.open() as f:r=csv.reader(f);h=next(r);rows=list(r)
 z=np.array([float(a.split(':')[1].split('cm')[0])/100 for a in h[1:-1]])
 t=np.array([stamp(r[0]) for r in rows]);T=np.array([[float(v) for v in r[1:-1]] for r in rows])
 if i=='2020b':z=z[:2];T=T[:,:2]
 return t,z,T

def node(t,top,L,a,initial):
 dt=t[1]-t[0]; rate=4*a/L**2; decay=np.exp(-rate*dt);phi=-np.expm1(-rate*dt)/(rate*dt)
 # exact response to linearly varying equilibrium top/2
 eq=top/2; y=np.empty(len(t));y[0]=initial
 # y[j]=decay*y[j-1]+(1-decay)*eq[j-1]+(1-phi)*(eq[j]-eq[j-1])
 force=(1-decay)*eq[:-1]+(1-phi)*np.diff(eq)
 y[1:]=lfilter([1],[1,-decay],force,zi=[decay*initial])[0]
 return y

def make(i,d,modes=MODES,domain_provenance=PROVENANCE):
 t,z,T=read(i); inside=(z>d*0+z[0])&(z<d); zz=z[inside]; YY=T[:,inside]
 initz=np.r_[z[0],zz,d];initT=np.r_[T[0,0],YY[0],0.]
 train=t<=t[0]+14*86400; tt=t[train]-t[0];score=(tt[::12]>=86400)&(tt[::12]<14*86400);target=YY[train][::12][score]
 prep=PreparedSlab(tt,T[train,0],np.zeros(train.sum()),zz,initz,initT,modes=modes,domain_provenance=PROVENANCE)
 grid=np.logspace(-8,-5,61);loss=[];nl=[];x=(zz-z[0])/(d-z[0]);xn=.5
 for a in grid:
  loss.append(float(np.mean((prep.predict(a)[::12][score]-target)**2)))
  n=node(tt,T[train,0],d-z[0],a,np.interp((z[0]+d)/2,initz,initT))
  pred=np.where(x[None,:]<=.5,T[train,0,None]*(1-2*x[None,:])+n[:,None]*(2*x[None,:]),n[:,None]*2*(1-x[None,:]))
  nl.append(float(np.mean((pred[::12][score]-target)**2)))
 a=grid[np.argmin(loss)];an=grid[np.argmin(nl)]
 full=PreparedSlab(t-t[0],T[:,0],np.zeros(len(t)),zz,initz,initT,modes=modes,domain_provenance=PROVENANCE)
 g=full.bottom_downward_gradient(a).gradient_K_m
 n=node(t-t[0],T[:,0],d-z[0],an,np.interp((z[0]+d)/2,initz,initT))
 return dict(t=t,grad={'pde':g,'quasi_steady':T[:,0]/(d-z[0]),'one_node':2*n/(d-z[0])},alpha={'pde':float(a),'one_node':float(an)},alpha_boundary={'pde':bool(np.argmin(loss) in [0,60]),'one_node':bool(np.argmin(nl) in [0,60])},excluded_depths_m=z[(z>=d)].tolist(),thickness=d)

def integral(t,g,start,end):
 lo=max(start,t[0]);hi=min(end,t[-1]);sel=(t>lo)&(t<hi);xx=np.r_[lo,t[sel],hi];yy=np.interp(xx,t,g)
 cross=np.where(yy[:-1]*yy[1:]<0)[0]
 if len(cross):
  zeros=xx[cross]+(xx[cross+1]-xx[cross])*(-yy[cross])/(yy[cross+1]-yy[cross]);oldx=xx;oldy=yy;xx=np.sort(np.r_[xx,zeros]);yy=np.interp(xx,oldx,oldy)
 return float(np.trapz(np.maximum(yy,0),xx)),float(np.trapz(np.minimum(yy,0),xx)),lo,hi

def raw2020():
 cols=[]
 for id in ['20457129','20457130']:
  with (ROOT/'original_2022'/(id+'.csv')).open(encoding='utf-8-sig') as f:r=csv.reader(f);next(r);next(r);a=list(r)
  cols.append((np.array([(datetime.strptime(x[1],'%m/%d/%y %H:%M:%S')+timedelta(hours=8)).replace(tzinfo=timezone.utc).timestamp() for x in a]),np.array([float(x[2]) for x in a])))
 common,a,b=np.intersect1d(cols[0][0],cols[1][0],return_indices=True);assert np.all(np.diff(common)==300);return common,cols[0][1][a],cols[1][1][b]

def run():
 output={'scope':'Conditional retrospective ablation transfer; not cliff/hydraulic field validation','LF_J_kg':LF,'modes':MODES,'source_target_comparisons':[], 'columns':[]}
 rt,r0,r1=raw2020()
 for geometry,ds in [('profile',PROFILE),('stake',STAKE)]:
  cases=[make(i,d) for i,d in zip(IDS,ds)]
  output['columns'] += [dict(profile=i,geometry=geometry,alpha=x['alpha'],alpha_boundary=x['alpha_boundary'],excluded_depths_m=x['excluded_depths_m']) for i,x in zip(IDS,cases)]
  for early in ['current','earlier']:
   for sh in [-12,0,12]:
    for eh in [-12,0,12]:
     integrals=[]
     for j,c in enumerate(cases):
      year=2020 if j==0 else 2023; s=datetime(year,7,25,8) if j==0 else datetime(year,6,8,8);e=datetime(year,9,6,8) if j==0 else datetime(year,9,2,8)
      start=(s+timedelta(hours=sh)).replace(tzinfo=timezone.utc).timestamp();end=(e+timedelta(hours=eh)).replace(tzinfo=timezone.utc).timestamp();rec={}
      for fam,g0 in c['grad'].items():
       t=c['t'];g=g0
       if j==0:
        # Original forcing covers initial mass period; early domain is separately initialized from actual two probes, then restarted at confirmed processed geometry. Explicit installation/initial-state scenario.
        sel=(rt>=start)&(rt<t[0]);topdepth=.05 if early=='current' else .005
        lt=rt[sel];ltop=r0[sel];lowerdepth=.15 if early=='current' else .18
        if fam=='quasi_steady':leading=ltop/(ds[j]-topdepth)
        elif fam=='one_node':
         initialnode=float(np.interp((topdepth+ds[j])/2,[topdepth,lowerdepth,ds[j]],[r0[sel][0],r1[sel][0],0]))
         leading=2*node(lt-lt[0],ltop,ds[j]-topdepth,c['alpha']['one_node'],initialnode)/(ds[j]-topdepth)
        else:
         leadingSlab=PreparedSlab(lt-lt[0],ltop,np.zeros(len(lt)),[lowerdepth],[topdepth,lowerdepth,ds[j]],[r0[sel][0],r1[sel][0],0.],modes=MODES,domain_provenance=PROVENANCE)
         leading=leadingSlab.bottom_downward_gradient(c['alpha']['pde']).gradient_K_m
        t=np.r_[rt[sel],t];g=np.r_[leading,g]
       pos,neg,lo,hi=integral(t,g,start,end);missing=max(0,lo-start)
       # First complete observed24h pattern repeated by source time of day; integral is phase-matched.
       first=c['t'][0];native=c['t'];gg=g0;mask=native<=first+86400
       if missing:
        missingtimes=np.arange(start,lo,300.);xx=np.r_[missingtimes,lo];phase=(xx-first)%86400;pat=np.interp(first+phase,native[mask],np.maximum(gg[mask],0));gap=float(np.trapz(pat,xx))
       else:gap=0.
       days=[integral(native,gg,first+k*86400,first+(k+1)*86400)[0] for k in range(14)];upper=max(days)*missing/86400
       rec[fam]=dict(covered=pos,negative=neg,point=pos+gap,upper=upper,missing_s=missing)
      rec['duration']=end-start;integrals.append(rec)
     for source in range(4):
      for target in range(4):
       if source==target:continue
       for fam in ['pde','quasi_steady','one_node','constant_daily']:
        if fam=='constant_daily':pred=MASS[source]*integrals[target]['duration']/integrals[source]['duration'];r=dict(predicted_full_kg_m2=pred,relative_error=pred/MASS[target]-1)
        else:
         src=integrals[source][fam];dst=integrals[target][fam];k=LF*MASS[source]/src['point'];pred=k*dst['point']/LF
         r=dict(k_source_W_m_K=k,source_effective_Cv_J_m3_K=(k/cases[source]['alpha'][fam] if fam in ['pde','one_node'] else None),predicted_full_kg_m2=pred,relative_error=pred/MASS[target]-1,predicted_covered_kg_m2=k*dst['covered']/LF,covered_exceeds_full=k*dst['covered']/LF>MASS[target],k_conditional_scenario_endpoints=[LF*MASS[source]/(src['covered']+src['upper']),LF*MASS[source]/src['covered']],target_missing_s=dst['missing_s'],target_negative_unit_gradient_integral=dst['negative'])
        output['source_target_comparisons'].append(dict(source=IDS[source],target=IDS[target],geometry=geometry,early_geometry=early,start_shift_h=sh,end_shift_h=eh,family=fam,observed_kg_m2=MASS[target],**r))
  Path(__file__).with_name('results.json' if MODES==64 else 'results_modes'+str(MODES)+'.json').write_text(json.dumps(output,indent=2));print('completed',geometry,flush=True)
if __name__=='__main__':run()
