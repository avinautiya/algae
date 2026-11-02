"""Conservative two-layer conduction; classical FV, no phase/flux calibration."""
import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.signal import lfilter

class PreparedLayered:
 def __init__(self,times_s,top,bottom,depths,initial_depths,initial_temperatures,cells=40,max_gap_s=300.):
  arrays=[np.array(x,dtype=float,copy=True) for x in [times_s,top,bottom,depths,initial_depths,initial_temperatures]]
  if any(x.ndim!=1 or not len(x) or not np.isfinite(x).all() for x in arrays):raise ValueError('finite nonempty vectors required')
  self.t,self.top,self.bottom,self.depths,z,T=arrays
  if len(self.t)!=len(self.top) or len(self.t)!=len(self.bottom) or len(z)!=len(T) or len(z)<2:raise ValueError('length mismatch')
  if np.any(np.diff(z)<=0) or np.any(np.diff(self.t)<=0):raise ValueError('ordered depths/times required')
  if not np.isfinite(max_gap_s) or max_gap_s<=0 or np.any(np.diff(self.t)>max_gap_s+1e-7):raise ValueError('forcing gap unsupported')
  if isinstance(cells,bool) or not isinstance(cells,int) or cells<4 or cells%2 or cells>512:raise ValueError('even cells4..512 required')
  self.z0=z[0];self.L=z[-1]-z[0];self.n=cells;self.dx=self.L/cells
  if np.any(self.depths<=z[0]) or np.any(self.depths>=z[-1]):raise ValueError('targets strictly interior')
  if abs(T[0]-self.top[0])>1e-12 or abs(T[-1]-self.bottom[0])>1e-12:raise ValueError('initial endpoint mismatch')
  self.x=z[0]+(np.arange(cells)+.5)*self.dx
  self.initial=np.interp(self.x,z,T);self.initial_targets=np.interp(self.depths,z,T)
  self.dt=np.diff(self.t);self.uniform=not len(self.dt) or np.allclose(self.dt,self.dt[0],rtol=1e-10,atol=1e-7)
 def operator(self,a,b):
  if not np.isfinite([a,b]).all() or min(a,b)<=0:raise ValueError('positive finite diffusivities required')
  alpha=np.r_[np.full(self.n//2,a),np.full(self.n//2,b)]
  face=np.empty(self.n+1);face[0]=2*a;face[-1]=2*b;face[1:-1]=2*alpha[:-1]*alpha[1:]/(alpha[:-1]+alpha[1:])
  diag=-(face[:-1]+face[1:])/self.dx**2;off=face[1:-1]/self.dx**2
  return diag,off,face
 def history(self,alpha_top,alpha_bottom):
  d,e,f=self.operator(alpha_top,alpha_bottom);w,V=eigh_tridiagonal(d,e)
  rate=-w
  # Algebraic steady profile for each boundary, harmonic resistance exactly.
  resistance=1/f;cum=np.cumsum(resistance)
  fraction=(cum[:-1])/cum[-1]
  # cum index includes top half-face resistance; interior cell j is after j+1 faces.
  P=1-fraction;B=fraction
  steady=self.top[:,None]*P+self.bottom[:,None]*B
  initial=self.initial-steady[0];c0=V.T@initial
  delta=np.diff(self.top)[:,None]*P+np.diff(self.bottom)[:,None]*B
  forcing=delta@V
  modes=np.empty((len(self.t),self.n));modes[0]=c0
  if len(self.dt):
   if self.uniform:
    exponent=rate*self.dt[0];r=np.exp(-exponent);phi=-np.expm1(-exponent)/exponent
    for j in range(self.n):modes[1:,j]=lfilter([1.],[1.,-r[j]],-phi[j]*forcing[:,j],zi=[r[j]*c0[j]])[0]
   else:
    for i,dt in enumerate(self.dt,1):
     exponent=rate*dt;modes[i]=np.exp(-exponent)*modes[i-1]+np.expm1(-exponent)/exponent*forcing[i-1]
  return steady+modes@V.T
 def predict(self,alpha_top,alpha_bottom):
  history=self.history(alpha_top,alpha_bottom)
  z=np.r_[self.z0,self.x,self.z0+self.L]
  full=np.column_stack([self.top,history,self.bottom])
  resistance=lambda zz: np.minimum(zz-self.z0,self.L/2)/alpha_top+np.maximum(zz-self.z0-self.L/2,0)/alpha_bottom
  rz=resistance(z);rd=resistance(self.depths)
  right=np.searchsorted(z,self.depths);left=right-1;weight=(rd-rz[left])/(rz[right]-rz[left])
  out=(1-weight)*full[:,left]+weight*full[:,right];out[0]=self.initial_targets
  if not np.isfinite(out).all():raise ValueError('nonfinite prediction')
  return out
