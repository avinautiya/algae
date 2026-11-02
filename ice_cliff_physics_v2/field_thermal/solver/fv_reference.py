"""Independent Crank-Nicolson finite-volume reference for solver checks only."""
import numpy as np
from scipy.linalg import solve_banded


def predict_fv(times_s,top,bottom,depths,alpha,initial_depths,initial_temperatures,
               cells=200,max_step_s=.005):
    times=np.asarray(times_s,float);top=np.asarray(top,float);bottom=np.asarray(bottom,float)
    z=np.asarray(initial_depths,float);initial=np.asarray(initial_temperatures,float)
    z0=z[0];L=z[-1]-z0;dx=L/cells
    centers=z0+(np.arange(cells)+.5)*dx
    state=np.interp(centers,z,initial);outputs=[np.interp(depths,z,initial)]
    for j in range(1,len(times)):
        span=times[j]-times[j-1];steps=max(1,int(np.ceil(span/max_step_s)));dt=span/steps
        r=alpha*dt/(2*dx*dx)
        band=np.zeros((3,cells));band[0,1:]=-r;band[2,:-1]=-r;band[1,:]=1+2*r;band[1,[0,-1]]=1+3*r
        for sub in range(steps):
            f0=sub/steps;f1=(sub+1)/steps
            ta=top[j-1]+f0*(top[j]-top[j-1]);tb=top[j-1]+f1*(top[j]-top[j-1])
            ba=bottom[j-1]+f0*(bottom[j]-bottom[j-1]);bb=bottom[j-1]+f1*(bottom[j]-bottom[j-1])
            rhs=state*(1-2*r);rhs[1:]+=r*state[:-1];rhs[:-1]+=r*state[1:]
            rhs[0]-=r*state[0];rhs[-1]-=r*state[-1]
            rhs[0]+=2*r*(ta+tb);rhs[-1]+=2*r*(ba+bb)
            state=solve_banded((1,1),band,rhs,check_finite=True)
        outputs.append(np.interp(depths,centers,state))
    return np.asarray(outputs)
