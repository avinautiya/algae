"""Small dense linear-response diagnostics, standard library only.

Computes ||C exp(A*t) B||_2^2 for an EXTERNALLY derived, dimensionless-state
linearization. This does not derive a glacier operator or establish novelty.
State/input/output normalization must be documented by the operator exporter.
Limited to matrices of dimension <=32; this is not a large CFD eigensolver.
"""
from math import ceil, isfinite, log2, sqrt


def matrix(a, name):
    if not a or not a[0] or len(a)>32 or len(a[0])>32:
        raise ValueError(f"{name}: nonempty dimensions <=32 required")
    cols=len(a[0])
    if any(len(row)!=cols for row in a) or any(not isfinite(x) for row in a for x in row):
        raise ValueError(f"{name}: finite rectangular matrix required")
    return [[float(x) for x in row] for row in a]


def identity(n):
    return [[float(i==j) for j in range(n)] for i in range(n)]


def multiply(a,b):
    if len(a[0])!=len(b):
        raise ValueError("Incompatible matrix dimensions")
    return [[sum(a[i][k]*b[k][j] for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def max_row_sum(a):
    return max(sum(abs(x) for x in row) for row in a)


def exponential(a,time):
    a=matrix(a,"A")
    if len(a)!=len(a[0]) or not isfinite(time) or time<0:
        raise ValueError("A must be square and time finite/nonnegative")
    norm=max_row_sum(a)*time
    if not isfinite(norm) or norm>1e6:
        raise ValueError("A*t is outside this small diagnostic's supported norm")
    scaling=max(0,int(ceil(log2(norm/0.5)))) if norm>0.5 else 0
    scaled=[[x*time/(2**scaling) for x in row] for row in a]
    result=identity(len(a)); term=identity(len(a))
    for k in range(1,101):
        term=[[x/k for x in row] for row in multiply(term,scaled)]
        result=[[x+y for x,y in zip(row,termrow)] for row,termrow in zip(result,term)]
        if max_row_sum(term)<=1e-16*max(1,max_row_sum(result)):
            break
    else:
        raise ValueError("Matrix exponential series did not converge")
    for _ in range(scaling):
        result=multiply(result,result)
        if any(not isfinite(x) for row in result for x in row):
            raise ValueError("Overflow in exponential; rescale/reduce time")
    return result


def largest_symmetric_eigenvalue(a):
    """Jacobi diagonalization, with explicit symmetry and residual controls."""
    a=matrix(a,"symmetric matrix"); n=len(a)
    if len(a[0])!=n:
        raise ValueError("Square matrix required")
    scale=max_row_sum(a)
    if any(abs(a[i][j]-a[j][i])>1e-12*max(scale,1e-300) for i in range(n) for j in range(n)):
        raise ValueError("Matrix is not symmetric")
    if n==1:
        return a[0][0]
    for _ in range(100*n*n):
        p,q=max(((i,j) for i in range(n) for j in range(i+1,n)),key=lambda ij:abs(a[ij[0]][ij[1]]))
        off=a[p][q]
        if abs(off)<=1e-14*max(scale,1e-300):
            return max(a[i][i] for i in range(n))
        tau=(a[q][q]-a[p][p])/(2*off)
        tangent=(1 if tau>=0 else -1)/(abs(tau)+sqrt(1+tau*tau))
        cosine=1/sqrt(1+tangent*tangent); sine=tangent*cosine
        app,aqq=a[p][p],a[q][q]
        for k in range(n):
            if k not in (p,q):
                akp,akq=a[k][p],a[k][q]
                a[k][p]=a[p][k]=cosine*akp-sine*akq
                a[k][q]=a[q][k]=sine*akp+cosine*akq
        a[p][p]=app-tangent*off; a[q][q]=aqq+tangent*off
        a[p][q]=a[q][p]=0.0
    raise ValueError("Symmetric eigensolve failed to reach tolerance")


def response_gain(a,time,b=None,c=None):
    """Squared singular gain from normalized input perturbations to outputs."""
    a=matrix(a,"A"); n=len(a)
    b=identity(n) if b is None else matrix(b,"B")
    c=identity(n) if c is None else matrix(c,"C")
    if len(b)!=n or len(c[0])!=n:
        raise ValueError("B rows and C columns must match the state dimension")
    m=multiply(multiply(c,exponential(a,time)),b)
    gram=[[sum(row[i]*row[j] for row in m) for j in range(len(m[0]))]
          for i in range(len(m[0]))]
    if any(not isfinite(x) for row in gram for x in row):
        raise ValueError("Overflow in singular-gain calculation")
    return max(0.0,largest_symmetric_eigenvalue(gram))


def analyze(data):
    required=("A_per_s","times_s","state_normalization","input_normalization","output_normalization","operator_provenance")
    if any(key not in data for key in required):
        raise ValueError("Missing operator normalization, provenance, matrix or time fields")
    if any(not isinstance(data[key],str) or not data[key].strip() for key in required[2:]):
        raise ValueError("Document nonempty normalization and provenance strings")
    times=data["times_s"]
    if not times or any(not isfinite(t) or t<0 for t in times):
        raise ValueError("Supply finite nonnegative times")
    gains=[dict(time_s=t,squared_gain=response_gain(data["A_per_s"],t,data.get("B"),data.get("C"))) for t in times]
    return dict(gains=gains,max_sampled_squared_gain=max(r["squared_gain"] for r in gains),
        operator_provenance=data["operator_provenance"],
        state_normalization=data["state_normalization"],input_normalization=data["input_normalization"],
        output_normalization=data["output_normalization"],
        limits=["Sampling times need refinement; sampled maximum is not a global maximum",
                "Gain depends on the declared norm and input/output maps",
                "No eigenvalue stability, CFD linearization, finite-amplitude validity or novelty is certified",
                "Neutral modes, conserved quantities and jump displacement require physical treatment before export"])
