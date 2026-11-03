"""Standard Thomas tridiagonal elimination; no third-party dependencies."""
import math


def tridiagonal(lower,diagonal,upper,rhs):
    n=len(diagonal);b=list(diagonal);r=list(rhs)
    for j in range(1,n):
        factor=lower[j-1]/b[j-1];b[j]-=factor*upper[j-1];r[j]-=factor*r[j-1]
    out=[0.]*n;out[-1]=r[-1]/b[-1]
    for j in range(n-2,-1,-1):out[j]=(r[j]-upper[j]*out[j+1])/b[j]
    if not all(math.isfinite(x) for x in out):raise ValueError('nonfinite linear solve')
    return out
