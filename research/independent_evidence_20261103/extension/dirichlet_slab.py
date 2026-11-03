"""Boundary-forced one-dimensional rock/debris heat diagnostic.

Both boundaries are measured probes, not assumed ice/surface temperatures.
Interior temperatures are predicted by constant-diffusivity heat conduction.
Negative temperatures are permitted. This is not a melting/enthalpy model.
"""
from __future__ import annotations
import math
from dataclasses import dataclass
import numpy as np
from scipy.signal import lfilter


def _vector(values,name):
    raw=np.asarray(values)
    if raw.ndim!=1 or raw.size==0 or raw.dtype.kind not in 'iuf':
        raise ValueError(name+' must be a nonempty one-dimensional numeric array')
    out=np.array(raw,dtype=float,copy=True)
    if not np.isfinite(out).all():raise ValueError(name+' has missing/nonfinite support')
    return out


def _positive(value,name):
    if isinstance(value,(bool,np.bool_)) or not isinstance(value,(int,float,np.integer,np.floating)) or not math.isfinite(value) or value<=0:
        raise ValueError('positive finite '+name+' required')
    return float(value)


def _phi(values):
    # expm1 retains the exact small-a limit without subtracting near-equals.
    values=np.asarray(values,float)
    out=np.ones_like(values)
    np.divide(-np.expm1(-values),values,out=out,where=values!=0)
    return out


@dataclass(frozen=True)
class Prediction:
    temperatures_C: np.ndarray
    diagnostics: dict


@dataclass(frozen=True)
class GradientPrediction:
    gradient_K_m: np.ndarray
    diagnostics: dict


@dataclass(frozen=True)
class IntervalGradientPrediction:
    integral_K_s_m: np.ndarray
    interval_mean_K_m: np.ndarray
    temperature_moment_K_m: np.ndarray
    diagnostics: dict


class PreparedSlab:
    """Prepare geometry/boundaries once, then predict without fitting alpha.

    Depth coordinates have units metres and must increase downwards. The domain
    is exactly [initial_depths[0], initial_depths[-1]], i.e. the shallowest and
    deepest actual supplied probes. Initial profile is their piecewise-linear
    interpolation at the FIRST supplied time only. No later interior values
    enter this object. Endpoints at time0 must agree with the supplied profile.
    """
    def __init__(self,times_s,top,bottom,depths,initial_depths,initial_temperatures,
                 modes=64,max_gap_s=300.,domain_provenance=None):
        self.times=_vector(times_s,'times_s');self.top=_vector(top,'top');self.bottom=_vector(bottom,'bottom')
        self.depths=_vector(depths,'target depths');self.initial_depths=_vector(initial_depths,'initial depths')
        self.initial_temperatures=_vector(initial_temperatures,'initial temperatures')
        if len(self.times)!=len(self.top) or len(self.times)!=len(self.bottom):raise ValueError('boundary/time lengths differ')
        if len(self.initial_depths)!=len(self.initial_temperatures) or len(self.initial_depths)<2:
            raise ValueError('initial profile needs matching depths/temperatures including both real endpoints')
        if np.any(np.diff(self.initial_depths)<=0) or np.any(np.diff(self.depths)<=0):raise ValueError('strictly ordered unique depths required')
        self.dt=np.diff(self.times)
        if np.any(self.dt<=0):raise ValueError('strictly increasing timestamps required')
        gap=_positive(max_gap_s,'maximum supported gap')
        tolerance=1e-7+gap*1e-10
        if np.any(self.dt>gap+tolerance):raise ValueError('forcing support gap exceeds declared maximum; no gap fill or extrapolation')
        if isinstance(modes,(bool,np.bool_)) or not isinstance(modes,(int,np.integer)) or not 1<=modes<=2048:
            raise ValueError('modes must be an integer from1 to2048')
        if domain_provenance is None:
            domain_provenance={'top_boundary':'measured_probe','bottom_boundary':'measured_probe','source':'caller-supplied real endpoint probe records'}
        if not isinstance(domain_provenance,dict) or set(domain_provenance)!={'top_boundary','bottom_boundary','source'} or any(not isinstance(v,str) or not v.strip() for v in domain_provenance.values()):
            raise ValueError('explicit top/bottom boundary and source provenance required')
        if domain_provenance['bottom_boundary']=='assumed_ice_Tm0' and np.any(self.bottom!=0):
            raise ValueError('assumed melting-point bottom must be explicitly zero')
        self.domain_provenance=dict(domain_provenance)
        self.modes=int(modes)
        if self.modes*len(self.times)>50_000_000:raise ValueError('requested spectral history exceeds bounded memory policy')
        self.z0=float(self.initial_depths[0]);self.z1=float(self.initial_depths[-1]);self.L=self.z1-self.z0
        if not math.isfinite(self.L) or self.L<=0:raise ValueError('invalid measured endpoint separation')
        if np.any(self.depths<=self.z0) or np.any(self.depths>=self.z1):raise ValueError('prediction depths must lie strictly between measured endpoint probes')
        endpoint_error=max(abs(self.initial_temperatures[0]-self.top[0]),abs(self.initial_temperatures[-1]-self.bottom[0]))
        temp_scale=max(1.,abs(self.top[0]),abs(self.bottom[0]))
        if endpoint_error>16*np.finfo(float).eps*temp_scale:raise ValueError('initial profile must belong to the FIRST boundary observation; no test-boundary reset')
        self.x=(self.depths-self.z0)/self.L
        knots=(self.initial_depths-self.z0)/self.L
        self.k=np.arange(1,self.modes+1,dtype=float)*math.pi
        self.sign=np.where(np.arange(1,self.modes+1)%2==1,1.,-1.)
        self.sines=np.sin(self.k[:,None]*self.x[None,:])
        self.linear=(1-self.x[None,:])*self.top[:,None]+self.x[None,:]*self.bottom[:,None]
        residual=self.initial_temperatures-((1-knots)*self.top[0]+knots*self.bottom[0])
        # Endpoint discrepancies within representation precision are normalized
        # to the prescribed Dirichlet values; magnitude is retained below.
        residual[0]=residual[-1]=0.
        slopes=np.diff(residual)/np.diff(knots)
        knot_sines=np.sin(self.k[:,None]*knots[None,:]);knot_sines[:,0]=0.;knot_sines[:,-1]=0.
        self.initial_coefficients=2/(self.k*self.k)*(np.diff(knot_sines,axis=1)@slopes)
        self.boundary_increments=(2/self.k[:,None])*(np.diff(self.top)[None,:]+self.sign[:,None]*np.diff(self.bottom)[None,:])
        self.initial_values=np.interp(self.depths,self.initial_depths,self.initial_temperatures)
        modal_initial=self.linear[0]+self.initial_coefficients@self.sines
        self.initial_modal_error=float(np.max(np.abs(modal_initial-self.initial_values)))
        self.endpoint_representation_error=float(endpoint_error)
        self.uniform=bool(len(self.dt)==0 or np.allclose(self.dt,self.dt[0],rtol=1e-10,atol=1e-7))
        arrays=[self.times,self.top,self.bottom,self.depths,self.initial_depths,self.initial_temperatures,self.dt,self.k,self.sign,self.sines,self.linear,self.initial_coefficients,self.boundary_increments,self.initial_values]
        if not all(np.isfinite(a).all() for a in arrays):raise ValueError('derived initialization or geometry overflow')
        for a in arrays:a.flags.writeable=False

    def _modal_history(self,alpha):
        alpha=_positive(alpha,'diffusivity in m2/s')
        rates=alpha*(self.k/self.L)**2
        if not np.isfinite(rates).all():raise ValueError('diffusion eigenvalues overflow')
        history=np.empty((self.modes,len(self.times)),float);history[:,0]=self.initial_coefficients
        if len(self.dt):
            if self.uniform:
                exponent=rates*self.dt[0]
                if not np.isfinite(exponent).all():raise ValueError('diffusion interval exponent overflow')
                decay=np.exp(-exponent);factor=_phi(exponent)
                for n in range(self.modes):
                    forcing=-factor[n]*self.boundary_increments[n]
                    history[n,1:]=lfilter([1.],[1.,-decay[n]],forcing,zi=[decay[n]*self.initial_coefficients[n]])[0]
            else:
                for j,dt in enumerate(self.dt,1):
                    exponent=rates*dt
                    if not np.isfinite(exponent).all():raise ValueError('diffusion interval exponent overflow')
                    history[:,j]=np.exp(-exponent)*history[:,j-1]-_phi(exponent)*self.boundary_increments[:,j-1]
        return alpha,history

    def predict_with_diagnostics(self,alpha):
        alpha,history=self._modal_history(alpha)
        values=self.linear+history.T@self.sines
        # t0 is an explicitly supplied initial condition, not a fitted/scored
        # assimilation. Modal truncation of that profile is separately exposed.
        values[0]=self.initial_values
        if not np.isfinite(values).all():raise ValueError('nonfinite modal prediction')
        return Prediction(values,dict(diffusivity_m2_s=alpha,modes=self.modes,
            uniform_time_grid=self.uniform,algorithm='exact per-mode recurrence for linearly interpolated endpoint forcing',
            domain_depths_m=[self.z0,self.z1],domain_length_m=self.L,domain_provenance=dict(self.domain_provenance),
            measured_domain_depths_m=[self.z0,self.z1] if self.domain_provenance['bottom_boundary']=='measured_probe' else None,
            initial_modal_reconstruction_max_error_C=self.initial_modal_error,
            endpoint_representation_normalization_C=self.endpoint_representation_error,
            initialization='all supplied probe values at first input timestamp only; no subsequent interior assimilation',
            boundary_support='piecewise-linear within validated adjacent timestamps only; no gap fill/extrapolation',
            scope='boundary-forced temperature diagnostic; no heat conductivity, flux, phase change, melt or glacier forecasting inferred'))

    @classmethod
    def with_assumed_ice_bottom(cls,times_s,top,depths,ice_depth_m,measured_initial_depths,
                               measured_initial_temperatures,provenance,modes=64,max_gap_s=300.):
        """Append an ASSUMED0 C ice boundary, never relabel it a measurement."""
        z=_vector(measured_initial_depths,'measured initial depths')
        initial=_vector(measured_initial_temperatures,'measured initial temperatures')
        ice=_positive(ice_depth_m,'assumed ice depth in metres')
        if len(z)!=len(initial) or np.any(np.diff(z)<=0) or ice<=z[-1]:
            raise ValueError('independently supplied ice depth must lie below all measured initial probes')
        if not isinstance(provenance,str) or not provenance.strip():raise ValueError('assumed ice-boundary provenance required')
        times=_vector(times_s,'times_s')
        return cls(times,top,np.zeros(len(times)),depths,np.r_[z,ice],np.r_[initial,0.],
            modes=modes,max_gap_s=max_gap_s,domain_provenance={'top_boundary':'measured_probe',
                'bottom_boundary':'assumed_ice_Tm0','source':provenance})

    def bottom_downward_gradient(self,alpha):
        """Signed q_bottom/k in K/m; no conductivity or latent heat inferred.

        This is a derivative of the finite modal solution, not a validated melt
        rate. Its convergence can be slower than interior temperatures.
        """
        alpha,history=self._modal_history(alpha)
        cosine=-self.sign
        gradient=(self.top-self.bottom)/self.L-history.T@(self.k/self.L*cosine)
        if not np.isfinite(gradient).all():raise ValueError('nonfinite bottom gradient')
        return GradientPrediction(gradient,dict(diffusivity_m2_s=alpha,modes=self.modes,
            units='K/m = q_downward / conductivity',sign='positive downwards into bottom; negative values retained',
            negative_samples=int(np.count_nonzero(gradient<0)),domain_provenance=dict(self.domain_provenance),
            measured_domain_depths_m=[self.z0,self.z1] if self.domain_provenance['bottom_boundary']=='measured_probe' else None,
            domain_depths_m=[self.z0,self.z1],
            caveat='instantaneous derivative Fourier tails converge slower than temperatures; no clipping, conductivity, phase or field-melt validation'))

    def bottom_interval_gradient(self,alpha):
        """Signed time-integrated q/k via an independent temperature moment.

        For coordinate z relative to the top probe, M=integral(z/L*T dz).
        Integral(q/k dt)=Integral((Ttop-Tbot)/L dt)-Delta(M)/alpha.
        The initial modal projection is retained, matching the propagated state;
        no temperature observation is assimilated after the first timestamp.
        """
        alpha,history=self._modal_history(alpha)
        moment=self.L*(self.top/6+self.bottom/3+history.T@(self.sign/self.k))
        quasi=(self.top-self.bottom)/self.L
        integral=.5*(quasi[1:]+quasi[:-1])*self.dt-np.diff(moment)/alpha
        mean=integral/self.dt
        if not np.isfinite(integral).all():raise ValueError('nonfinite moment-integrated bottom gradient')
        return IntervalGradientPrediction(integral,mean,moment,dict(diffusivity_m2_s=alpha,modes=self.modes,
            units='K*s/m = time-integrated q_downward / conductivity',sign='signed net heat; negative intervals retained',
            negative_intervals=int(np.count_nonzero(integral<0)),net_integral_K_s_m=float(np.sum(integral)),
            domain_provenance=dict(self.domain_provenance),domain_depths_m=[self.z0,self.z1],
            operator='linear-boundary quasi-steady integral minus exact modal temperature-moment increment / alpha',
            caveat='not a positive-only integral or latent-melt estimate; finite-temperature-mode convergence must still be checked'))

    def predict(self,alpha):
        return self.predict_with_diagnostics(alpha).temperatures_C


def predict(times_s,top,bottom,depths,alpha,initial_depths,initial_temperatures,
            modes=64,max_gap_s=300.):
    """Return [time,interior_depth] Celsius predictions without fitting alpha."""
    return PreparedSlab(times_s,top,bottom,depths,initial_depths,initial_temperatures,
                        modes=modes,max_gap_s=max_gap_s).predict(alpha)
