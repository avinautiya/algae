"""Analytical prediction-design checks, not a CFD or glacier prediction solver.

Uses only Python 3.9+ standard library. Outputs are reproducible synthetic checks.
"""
from fractions import Fraction as F
from pathlib import Path
import csv
import json
import math

ROOT = Path(__file__).resolve().parent
RHO_I, RHO_W, LF, CP = 917.0, 1000.0, 334000.0, 4184.0
G, NU, K_W, SIGMA = 9.81, 1.8e-6, 0.56, 0.075


def integral(p):
    return sum((c/F(i+1) for i,c in enumerate(p)), F(0))


def multiply(p,q):
    r=[F(0)]*(len(p)+len(q)-1)
    for i,a in enumerate(p):
        for j,b in enumerate(q):r[i+j]+=a*b
    return r


def evaluate(p,y):
    return sum(float(c)*y**i for i,c in enumerate(p))


def closure_counterexample():
    # y in [0,1], no-slip/parabolic stress-free velocity shape w=2y-y².
    # Both temperature profiles have equal depth mean, flow-weighted mean,
    # ice-interface temperature, top temperature and top conductive heat flux.
    base=[F(0),F(2),F(-1)]
    f=[F(0),F(1),-F(17,2),F(22),-F(45,2),F(8)]
    assert integral(f)==0
    assert integral(multiply(base,f))==0
    assert sum(f)==0
    assert sum(i*c for i,c in enumerate(f))==0
    eps=F(1,2)
    profiles=[]
    for sign in (-1,1):
        p=[(base[i] if i<len(base) else F(0))+sign*eps*f[i] for i in range(len(f))]
        # Positivity also follows from bound in the derivation; sampling catches coding errors.
        assert min(evaluate(p,j/10000) for j in range(10001))>=-1e-12
        profiles.append(p)
    assert integral(profiles[0])==integral(profiles[1])==F(2,3)
    assert integral(multiply(base,profiles[0]))/integral(base)==F(4,5)
    assert integral(multiply(base,profiles[1]))/integral(base)==F(4,5)
    result=dict(depth_mean_over_scale=str(F(2,3)),flux_mean_over_scale=str(F(4,5)),
                bottom_gradient_over_scale=[str(p[1]) for p in profiles],
                flux_ratio=float(profiles[1][1]/profiles[0][1]),
                warning='Constructed admissible instantaneous fields, not evidence of a persistent environmental effect.')
    rows=[dict(y_over_depth=j/100,T_minus=evaluate(profiles[0],j/100),T_plus=evaluate(profiles[1],j/100)) for j in range(101)]
    return result,rows


def error_bound(e0, flux_error, lipschitz_flux, duration_s):
    """Conditional bound: K is W/m³, epsilon_q W/m²; see derivation assumptions."""
    rate=lipschitz_flux/(RHO_I*LF)
    if lipschitz_flux==0:return e0+flux_error*duration_s/(RHO_I*LF)
    return e0*math.exp(rate*duration_s)+flux_error/lipschitz_flux*math.expm1(rate*duration_s)


def projection_check():
    # Exchange rates in watts: exactly opposite terms must use same computed value.
    incoming,outgoing,air,mechanical,storage,conduction=12000.,4000.,500.,100.,600.,2000.
    water_to_ice=incoming-outgoing+air+mechanical-storage
    latent=water_to_ice-conduction
    residual=incoming-outgoing+air+mechanical-storage-conduction-latent
    assert residual==0
    m3_per_s=latent/(RHO_I*LF)
    assert abs(RHO_W*(RHO_I/RHO_W)*m3_per_s-latent/LF)<1e-14
    return dict(water_to_interface_W=water_to_ice,latent_W=latent,
                ice_volume_m3_s=m3_per_s,ledger_residual_W=residual,
                warning='Manufactured ledger identity; not a discrete moving-mesh conservation test.')


def regime(h,u,width,theta):
    c=math.sqrt(G*h)
    return dict(depth_m=h,velocity_m_s=u,width_m=width,temperature_scale_K=theta,
                Fr=u/c,Re_Uh=u*h/NU,We_h=RHO_W*u*u*h/SIGMA,
                Bo_h=RHO_W*G*h*h/SIGMA,Ste_sensible_over_latent=CP*theta/LF,
                inverse_Ste=LF/(CP*theta),
                upstream_characteristic_m_s=u-c,
                upstream_small_disturbance_possible=(u-c<0),
                latent_bias_for_2cm_per_day_W_m2=RHO_I*LF*0.02/86400)


def main():
    counterexample,profiles=closure_counterexample()
    with (ROOT/'temperature_counterexample.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(profiles[0]));w.writeheader();w.writerows(profiles)
    # For the linear worst-case scalar problem the bound is attained exactly.
    assert abs(error_bound(0,100,0,86400)-100*86400/(RHO_I*LF))<1e-15
    assert abs(error_bound(0,100,1e-12,86400)/error_bound(0,100,0,86400)-1)<1e-12
    target_q=1.329e-5
    examples={
        'published_flume_hydraulics_only':regime(.00165,target_q/(.02*.00165),.02,.7),
        'synthetic_supercritical_reach':regime(.05,1,.5,.2),
        'synthetic_subcritical_reach':regime(.1,.5,.5,.2)}
    assert examples['synthetic_supercritical_reach']['upstream_small_disturbance_possible'] is False
    assert examples['synthetic_subcritical_reach']['upstream_small_disturbance_possible'] is True
    result=dict(status='PASS',scope='Analytical and dimensional design checks only; not a predictive CFD result.',
                constructed_closure_counterexample=counterexample,energy_projection=projection_check(),
                regime_examples=examples,
                conditional_bed_error_m_at_100W_m2_24h=error_bound(0,100,0,86400),
                warning='Example material constants are approximate design values, not site calibrations.')
    (ROOT/'design_check_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
