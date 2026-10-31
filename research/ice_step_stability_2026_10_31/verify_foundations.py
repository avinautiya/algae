"""Small verification calculations; no CFD, eigenvalue campaign or field prediction.

Python 3.9+ standard library only. Run from any directory.
"""
import cmath
import csv
import json
import math
import random
from pathlib import Path

OUT = Path(__file__).resolve().parent
G = 9.81


def thermal_response(s, k, depth, diffusivity):
    """Dirichlet-to-Neumann multiplier p*coth(p*D), units 1/m.

    Constant-property flat slab, fixed lower temperature; Re(s)>=0 here.
    Not a substitute for the operator on a curved traveling ice boundary.
    """
    p = cmath.sqrt(k*k + s/diffusivity)
    x = p*depth
    if abs(x) < 1e-4:
        return (1 + x*x/3 - x**4/45 + 2*x**6/945)/depth
    if x.real > 30:
        return p
    return p/cmath.tanh(x)


def fd_response(s, k, depth, diffusivity, n):
    """Independent second-order finite-difference boundary-value solve.

    psi(-D)=0, psi(0)=1, psi''=(k^2+s/kappa)*psi.
    Tridiagonal elimination; second-order one-sided boundary derivative.
    """
    dz = depth/n
    a = k*k+s/diffusivity
    diag = [complex(2+a*dz*dz)]*(n-1)
    rhs = [0j]*(n-1)
    rhs[-1] = 1
    for j in range(1, n-1):
        factor = -1/diag[j-1]
        diag[j] += factor
        rhs[j] -= factor*rhs[j-1]
    y = [0j]*(n-1)
    y[-1] = rhs[-1]/diag[-1]
    for j in range(n-3, -1, -1):
        y[j] = (rhs[j]+y[j+1])/diag[j]
    return (3-4*y[-1]+y[-2])/(2*dz)


def flux(v):
    h, q = v
    return [q, q*q/h+0.5*G*h*h]


def flux_jacobian_action(v, dv):
    h, q = v
    dh, dq = dv
    return [dq, (G*h-q*q/(h*h))*dh+2*q/h*dq]


def shock_checks():
    rng = random.Random(619)
    errors, ratios = [], []
    energy_errors = []
    for _ in range(200):
        hm = rng.uniform(0.02, 0.3)
        fr = rng.uniform(1.2, 5)
        c = rng.uniform(-0.2, 0.2)
        um_rel = fr*math.sqrt(G*hm)
        hp = hm*0.5*(math.sqrt(1+8*fr*fr)-1)
        j = hm*um_rel
        states = [[hm, j+c*hm], [hp, j+c*hp]]
        speeds = [v[1]/v[0]-c for v in states]
        energy_difference = (hm+speeds[0]**2/(2*G))-(hp+speeds[1]**2/(2*G))
        exact_loss = (hp-hm)**3/(4*hp*hm)
        energy_errors.append(abs(energy_difference/exact_loss-1))
        perturb = [[rng.uniform(-0.01, 0.01), rng.uniform(-0.01, 0.01)] for _ in range(2)]
        gradients = [[rng.uniform(-0.1, 0.1), rng.uniform(-0.1, 0.1)] for _ in range(2)]
        zeta, sigma = rng.uniform(-0.2, 0.2), rng.uniform(-1, 1)
        dv = [[p+zeta*g for p,g in zip(v,w)] for v,w in zip(perturb,gradients)]
        act = [flux_jacobian_action(v,d) for v,d in zip(states,dv)]
        derivative = [(act[1][i]-c*dv[1][i])-(act[0][i]-c*dv[0][i])
                      -sigma*zeta*(states[1][i]-states[0][i]) for i in range(2)]
        def residual(eps):
            vals = [[x+eps*y for x,y in zip(v,d)] for v,d in zip(states,dv)]
            s = c+eps*sigma*zeta
            f = [flux(v) for v in vals]
            return [f[1][i]-s*vals[1][i]-f[0][i]+s*vals[0][i] for i in range(2)]
        eps = 1e-5
        rp, rm = residual(eps), residual(-eps)
        fd = [(a-b)/(2*eps) for a,b in zip(rp,rm)]
        errors.append(max(abs(x-y) for x,y in zip(fd,derivative))/max(1,*(abs(x) for x in derivative)))
        r0 = residual(0)
        ratios.append(max(abs(x) for x in r0))
    assert max(errors) < 1e-8
    assert max(energy_errors) < 1e-10
    assert max(ratios) < 1e-12
    return dict(cases=200, max_scaled_linearization_error=max(errors),
                max_base_jump_residual=max(ratios), max_relative_jump_head_loss_error=max(energy_errors))


def thermal_checks():
    # Nondimensional examples: kappa=D=1. Physical units restored in screen.
    cases = [(0j, 0), (0j, 2), (1j, 0), (10j, 1), (3+10j, 2)]
    rows = []
    for s,k in cases:
        target = thermal_response(s,k,1,1)
        for n in (32,64,128,256,512):
            got = fd_response(s,k,1,1,n)
            rows.append(dict(s_real=s.real,s_imag=s.imag,k=k,n=n,
                             relative_error=abs(got-target)/abs(target)))
    finest = [r['relative_error'] for r in rows if r['n']==512]
    assert max(finest)<5e-5
    nontrivial = [r for r in rows if (r['s_real'],r['s_imag'],r['k']) != (0,0,0)]
    orders=[]
    for s,k in cases[1:]:
        seq=[r for r in nontrivial if r['s_real']==s.real and r['s_imag']==s.imag and r['k']==k]
        orders.append(math.log(seq[-2]['relative_error']/seq[-1]['relative_error'],2))
    assert min(orders)>1.95
    rng=random.Random(620)
    positivity=[]
    for _ in range(1000):
        s=complex(10**rng.uniform(-5,3),rng.uniform(-100,100))
        positivity.append(thermal_response(s,10**rng.uniform(-4,1),10**rng.uniform(-2,1),1).real)
    assert min(positivity)>0
    assert abs(thermal_response(0,0,0.2,1)-5)<1e-12
    return rows, dict(finest_relative_error=max(finest),observed_orders=orders,
                      passivity_samples=1000,min_sampled_real_response=min(positivity),
                      note='The positivity proof is in MATHEMATICS.md; sampling is only a regression check.')


def dimensional_screen():
    # Deliberately illustrative, NOT measured Himalayan boundary conditions.
    kappa=1.2e-6
    rows=[]
    for d in (0.02,0.08,0.2,0.5):
        for period_h in (1,6,24):
            omega=2*math.pi/(period_h*3600)
            r=thermal_response(1j*omega,0,d,kappa)
            rows.append(dict(depth_m=d,period_hours=period_h,
                             diffusion_time_hours=d*d/kappa/3600,
                             omega_D2_over_kappa=omega*d*d/kappa,
                             penetration_depth_m=math.sqrt(2*kappa/omega),
                             gain_relative_to_static=abs(r)*d,
                             phase_deg=cmath.phase(r)*180/math.pi))
    return rows


def write_csv(name, rows):
    with (OUT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


if __name__=='__main__':
    rows,thermal=thermal_checks()
    write_csv('thermal_refinement.csv',rows)
    write_csv('illustrative_timescales.csv',dimensional_screen())
    result=dict(status='PASS',scope='Verification of analytical subproblems only; no ice-step stability result.',
                thermal_response=thermal,moving_hydraulic_jump=shock_checks())
    (OUT/'verification_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
