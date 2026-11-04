"""Bounded thermal response inference. Never reads ablation or outcome loaders."""
import cmath,math
OMEGA=2*math.pi/86400

def transfer(alpha,b,omega,L=.17,x=.10):
 if not all(math.isfinite(v) for v in [alpha,b,omega,L,x]) or alpha<=0 or b<0 or omega<0 or not 0<=x<=L:raise ValueError('Invalid thermal response parameters')
 if omega==0:return complex((L-x+b)/(L+b),0)
 g=cmath.sqrt(1j*omega/alpha)
 return (cmath.sinh(g*(L-x))+b*g*cmath.cosh(g*(L-x)))/(cmath.sinh(g*L)+b*g*cmath.cosh(g*L))

def invert(a):
 n=len(a);v=[list(map(float,row))+[float(i==j) for j in range(n)] for i,row in enumerate(a)]
 for j in range(n):
  i=max(range(j,n),key=lambda i:abs(v[i][j]))
  if abs(v[i][j])<1e-14:raise ValueError('Rank-deficient OLS design')
  v[i],v[j]=v[j],v[i];p=v[j][j];v[j]=[x/p for x in v[j]]
  for i in range(n):
   if i==j:continue
   f=v[i][j];v[i]=[x-f*y for x,y in zip(v[i],v[j])]
 return [row[n:] for row in v]

def harmonic_extract(times,upper,lower,accuracy_C=.2):
 if len(times)!=len(upper) or len(times)!=len(lower) or len(times)<12:raise ValueError('Paired thermal series required')
 if any(not math.isfinite(x) for a in [times,upper,lower] for x in a) or any(b<=a for a,b in zip(times,times[1:])):raise ValueError('Finite strictly increasing times required')
 if not math.isfinite(accuracy_C) or accuracy_C<=0:raise ValueError('Positive source error bound required')
 t0=times[0];mid=(times[0]+times[-1])/2
 X=[[1,(t-mid)/86400,math.cos(OMEGA*(t-t0)),math.sin(OMEGA*(t-t0)),math.cos(2*OMEGA*(t-t0)),math.sin(2*OMEGA*(t-t0))] for t in times]
 normal=[[sum(r[i]*r[j] for r in X) for j in range(6)] for i in range(6)];inv=invert(normal);W=[[sum(inv[i][j]*r[j] for j in range(6)) for r in X] for i in range(6)];cu=[sum(w*v for w,v in zip(row,upper)) for row in W];cy=[sum(w*v for w,v in zip(row,lower)) for row in W]
 def pairfreq(i):return complex(cu[i],-cu[i+1]),complex(cy[i],-cy[i+1]),accuracy_C*math.hypot(sum(abs(w) for w in W[i]),sum(abs(w) for w in W[i+1]))
 obs=[]
 for omega,U,V,e in [(0,complex(cu[0]),complex(cy[0]),accuracy_C*sum(abs(w) for w in W[0])),(OMEGA,*pairfreq(2)),(2*OMEGA,*pairfreq(4))]:
  if abs(U)<=e:raise ValueError('Input harmonic/mean not resolved by declared source-error bound')
  ratio=V/U;scale=(e+abs(ratio)*e)/(abs(U)-e)
  obs.append({'omega_rad_s':omega,'ratio_real':ratio.real,'ratio_imag':ratio.imag,'upper_amplitude_C':abs(U),'lower_amplitude_C':abs(V),'coefficient_complex_error_bound_C':e,'ratio_error_scale':scale})
 residual={}
 for key,y,c in [('upper',upper,cu),('lower',lower,cy)]:
  v=[sum(a*b for a,b in zip(row,c))-val for row,val in zip(X,y)];residual[key]={'RMSE_C':math.sqrt(sum(x*x for x in v)/len(v)),'max_absolute_C':max(abs(x) for x in v)}
 return {'observations':obs,'upper_coefficients':cu,'lower_coefficients':cy,'basis':'constant,centered trendC/day,cos1,sin1,cos2,sin2','fit_reconstruction':residual,'n':len(times),'source_accuracy_C':accuracy_C,'scale_interpretation':'deterministic ratio perturbation scale, not a standard deviation'}

def loss(logalpha,b,obs):
 alpha=math.exp(logalpha);total=0.
 for o in obs:
  h=transfer(alpha,b,o['omega_rad_s']);r=complex(o['ratio_real'],o['ratio_imag']);total+=abs(h-r)**2/o['ratio_error_scale']**2
 return total

def bounded_fit(observations,*,null=False):
 amin,amax=math.log(1e-8),math.log(5e-6);blo,bhi=0.,.3;lo,hi=amin,amax;best=None
 for level in range(4):
  na=61 if level==0 else 21;nb=1 if null else (61 if level==0 else 21);da=(hi-lo)/(na-1);db=0 if null else (bhi-blo)/(nb-1)
  for i in range(na):
   a=lo+i*da
   for j in range(nb):
    b=0. if null else blo+j*db;score=loss(a,b,observations)
    if best is None or (score,a,b)<(best[0],best[1],best[2]):best=(score,a,b)
  lo=max(amin,best[1]-da);hi=min(amax,best[1]+da)
  if not null:blo=max(0.,best[2]-db);bhi=min(.3,best[2]+db)
 alpha=math.exp(best[1]);b=best[2];a_edge=min(abs(best[1]-amin),abs(best[1]-amax))<1e-6;b_edge=null or min(b,.3-b)<1e-7
 return {'alpha_m2_s':alpha,'b_kR_m':b,'weighted_objective':best[0],'alpha_at_bound':a_edge,'b_at_bound':b_edge,'bounds_engineering_not_measurement':True}

def jacobian_rank(alpha,b,obs,h=1e-4,L=.17):
 def vec(a,bb):
  v=[]
  for o in obs:
   z=transfer(a,bb,o['omega_rad_s'])/o['ratio_error_scale'];v.append(z.real)
   if o['omega_rad_s']:v.append(z.imag)
  return v
 va,vb=vec(alpha*math.exp(h),b),vec(alpha*math.exp(-h),b);j1=[(a-z)/(2*h) for a,z in zip(va,vb)]
 lower=max(0,b-h*L);upper=b+h*L;vp,vm=vec(alpha,upper),vec(alpha,lower);j2=[(a-z)/((upper-lower)/L) for a,z in zip(vp,vm)]
 a=sum(x*x for x in j1);d=sum(x*x for x in j2);c=sum(x*y for x,y in zip(j1,j2));disc=math.sqrt(max(0,(a-d)**2+4*c*c));l1=(a+d+disc)/2;l2=max(0,(a+d-disc)/2)
 return {'coordinates':['logalpha','b/L'],'singular_values':[math.sqrt(l1),math.sqrt(l2)],'condition_number':math.sqrt(l1/l2) if l2 else None,'rank_scope':'local response sensitivity conditional on fixedgeometry/model/error scales; not global physical identifiability or confidence'}
