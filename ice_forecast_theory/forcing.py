"""Finite-horizon response to piecewise-constant, budget-constrained forcing.

Known linear-systems mathematics, implemented for future physically exported
operators. The current examples are manufactured; no glacier operator is fitted.
Depends only on the project's previously verified ice_cfd_checks package.
"""
from math import isfinite, sqrt

from ice_cfd_checks.linear_response import exponential, identity, matrix, max_row_sum, multiply


def transpose(a):
    return [list(x) for x in zip(*a)]


def add(a,b):
    return [[x+y for x,y in zip(ar,br)] for ar,br in zip(a,b)]


def zero_order_hold(a,b,dt):
    """Return Phi, Gamma without inverting A, including singular/neutral modes."""
    a=matrix(a,"A"); b=matrix(b,"B")
    n,m=len(a),len(b[0])
    if len(a[0])!=n or len(b)!=n or n+m>32:
        raise ValueError("A/B dimensions disagree or augmented dimension exceeds 32")
    if not isfinite(dt) or dt<=0:
        raise ValueError("dt must be finite and positive")
    augmented=[a[i]+b[i] for i in range(n)]+[[0.0]*(n+m) for _ in range(m)]
    propagator=exponential(augmented,dt)
    return [r[:n] for r in propagator[:n]],[r[n:] for r in propagator[:n]]


def constraint_projector(rows,dimension):
    """Orthogonal projector onto R*z=0 with redundant rows handled by QR.

Two-pass modified Gram-Schmidt, scaled by each row's maximum entry. Constraints
must already be expressed in the normalized forcing coordinates. Near-dependent
rows require a sensitivity check; this is a numerical projector, not a proof.
"""
    basis=[]
    for row in rows:
        if len(row)!=dimension or any(not isfinite(x) for x in row):
            raise ValueError("Constraint dimensions or values invalid")
        scale=max(map(abs,row),default=0.0)
        if scale==0:
            continue
        v=[x/scale for x in row]
        for _ in range(2):
            for q in basis:
                projection=sum(a*b for a,b in zip(v,q))
                v=[a-projection*b for a,b in zip(v,q)]
        norm=sqrt(sum(x*x for x in v))
        if norm>1e-12:
            basis.append([x/norm for x in v])
    p=identity(dimension)
    for q in basis:
        p=[[p[i][j]-q[i]*q[j] for j in range(dimension)] for i in range(dimension)]
    return p,len(basis)


def dominant_symmetric_pair(a):
    """Jacobi eigenvectors for a small symmetric positive-semidefinite Gram matrix."""
    a=matrix(a,"Gram matrix"); n=len(a)
    if len(a[0])!=n:
        raise ValueError("Gram matrix must be square")
    vectors=identity(n); scale=max_row_sum(a)
    if any(abs(a[i][j]-a[j][i])>1e-12*max(scale,1e-300) for i in range(n) for j in range(n)):
        raise ValueError("Gram matrix is not symmetric")
    for _ in range(100*n*n):
        if n==1:
            break
        p,q=max(((i,j) for i in range(n) for j in range(i+1,n)),key=lambda ij:abs(a[ij[0]][ij[1]]))
        off=a[p][q]
        if abs(off)<=1e-14*max(scale,1e-300):
            break
        tau=(a[q][q]-a[p][p])/(2*off)
        t=(1 if tau>=0 else -1)/(abs(tau)+sqrt(1+tau*tau))
        c=1/sqrt(1+t*t); s=t*c; app,aqq=a[p][p],a[q][q]
        for k in range(n):
            if k not in (p,q):
                kp,kq=a[k][p],a[k][q]
                a[k][p]=a[p][k]=c*kp-s*kq
                a[k][q]=a[q][k]=s*kp+c*kq
            vp,vq=vectors[k][p],vectors[k][q]
            vectors[k][p]=c*vp-s*vq; vectors[k][q]=s*vp+c*vq
        a[p][p]=app-t*off; a[q][q]=aqq+t*off; a[p][q]=a[q][p]=0.0
    else:
        raise ValueError("Gram eigensolve did not converge")
    largest=max(range(n),key=lambda i:a[i][i])
    return max(a[largest][largest],0.0),[vectors[k][largest] for k in range(n)]


def optimal_forcing(data):
    """Optimize zero-initial-state, time-average-L2 forcing under linear budgets.

A_j/B_j act on consistently normalized states and input amplitudes throughout.
All intervals share one state basis and the same number/order of input channels.
No amplitude box constraints or nonlinear physical validity are certified.
"""
    for key in ("intervals","C","state_normalization","input_normalization","output_normalization","operator_provenance"):
        if key not in data:
            raise ValueError(f"Missing {key}")
    for key in ("state_normalization","input_normalization","output_normalization","operator_provenance"):
        if not isinstance(data[key],str) or not data[key].strip():
            raise ValueError("Nonempty normalization and provenance strings required")
    intervals=data["intervals"]
    if not intervals:
        raise ValueError("No forcing intervals")
    total=sum(r["duration_s"] for r in intervals)
    if not isfinite(total) or total<=0:
        raise ValueError("Invalid total horizon")
    n=len(intervals[0]["A_per_s"]); m=len(intervals[0]["B_per_s"][0])
    dimension=m*len(intervals)
    if dimension>32:
        raise ValueError("This dense utility supports at most 32 forcing coordinates")
    response=None; blocks=[]; scales=[]
    for interval in intervals:
        a=matrix(interval["A_per_s"],"A"); b=matrix(interval["B_per_s"],"B")
        if len(a)!=n or len(a[0])!=n or len(b)!=n or len(b[0])!=m:
            raise ValueError("State/input dimensions or basis ordering changed")
        dt=interval["duration_s"]
        phi,gamma=zero_order_hold(a,b,dt)
        scale=sqrt(total/dt)
        scaled=[[scale*x for x in row] for row in gamma]
        if response is None:
            response=scaled
        else:
            response=[r+s for r,s in zip(multiply(phi,response),scaled)]
        scales.extend([scale]*m); blocks.append((phi,gamma))
    c=matrix(data["C"],"C")
    if len(c[0])!=n:
        raise ValueError("Output map dimensions disagree")
    h=multiply(c,response)
    raw_constraints=data.get("budget_rows_on_input_amplitudes",[])
    for row in raw_constraints:
        if len(row)!=dimension or any(not isfinite(x) for x in row):
            raise ValueError("Budget rows must match all interval/channel amplitudes")
    normalized=[[r*s for r,s in zip(row,scales)] for row in raw_constraints]
    p,rank=constraint_projector(normalized,dimension)
    constrained=multiply(h,p)
    gram=multiply(transpose(constrained),constrained)
    value,vector=dominant_symmetric_pair(gram)
    z=[row[0] for row in multiply(p,[[v] for v in vector])]
    norm=sqrt(sum(v*v for v in z))
    if value==0 or norm<1e-12:
        z=[0.0]*dimension
    else:
        z=[v/norm for v in z]
    amplitudes=[v*s for v,s in zip(z,scales)]
    # Independent state recurrence replays the returned input sequence.
    state=[[0.0] for _ in range(n)]
    for j,(phi,gamma) in enumerate(blocks):
        inputs=[[v] for v in amplitudes[j*m:(j+1)*m]]
        state=add(multiply(phi,state),multiply(gamma,inputs))
    output=[row[0] for row in multiply(c,state)]
    realized=sum(v*v for v in output)
    constraints=[sum(r*a for r,a in zip(row,amplitudes)) for row in raw_constraints]
    return dict(horizon_s=total,forcing_norm="(1/T) integral |f(t)|^2 dt",
        squared_gain=value,linear_norm_bound_per_unit_input=sqrt(value),
        optimal_interval_amplitudes=[amplitudes[j*m:(j+1)*m] for j in range(len(intervals))],
        normalized_optimal_coordinates=z,realized_output=output,realized_squared_output=realized,
        gain_replay_difference=realized-value,budget_residuals=constraints,budget_rank=rank,
        admissible_coordinate_dimension=dimension-rank,
        normalized_response_matrix=h,operator_provenance=data["operator_provenance"],
        normalization={k:data[k] for k in ("state_normalization","input_normalization","output_normalization")},
        limitations=["Reduced linear model only; not a physical forecast until operator and nonlinear response are validated",
            "Budget constraints are linearized; nonlinear heat-flow products require explicit paired-case checks",
            "Optimal signed perturbations need amplitude and positivity bounds before physical use",
            "All intervals require a common state basis; no hidden projection/remap is allowed",
            "No claim of a new optimization theorem or an actual glacier discovery"])
