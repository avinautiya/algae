"""Independent ideal-control counterexample, not field validation or a new theorem.
No project files are imported or changed. Python standard library only.
"""
import json, math
from pathlib import Path

def partition(Q,g,K):
    """Bypass crest0; toe crestg. Returns actual mass-balanced ideal split."""
    lo=min(0.,g);hi=max(0.,g)+(Q/K)**(2/3)
    for _ in range(100):
        eta=(lo+hi)/2;total=K*max(eta-g,0.)**1.5+K*max(eta,0.)**1.5
        if total<Q:lo=eta
        else:hi=eta
    eta=(lo+hi)/2
    return dict(toe_m3_s=K*max(eta-g,0.)**1.5,bypass_m3_s=K*max(eta,0.)**1.5,stage_m=eta)

def main():
    Q0=.07;K=(2/3)**1.5*math.sqrt(9.81)*.5;g=.25;H0=(Q0/K)**(2/3)
    before=partition(Q0,g,K);pulse=partition(2*Q0,g,K);after=partition(Q0,g,K)
    rhoL=917*334000.;h=200.;Tin=.5;maxE=h*Tin/rhoL
    rows=[]
    for ratio in [.25,.5,.9,1,1.1,1.5,2,5]:
        gap=ratio*H0;r=partition(Q0,gap,K)
        rows.append(dict(gap_over_critical_head=ratio,toe_flow_fraction=r['toe_m3_s']/Q0,
           algebraic_contact=gap<H0,root_contact=r['toe_m3_s']>1e-15))
    assert all(r['algebraic_contact']==r['root_contact'] for r in rows)
    assert before['toe_m3_s']==0 and pulse['toe_m3_s']>0 and after['toe_m3_s']==0
    deficit=g-H0;duration=3600.
    assert maxE*duration<deficit
    out=dict(scope='exact ideal-controls check and necessary bound only',K_mhalf_s=K,
       baseline_Q_m3_s=Q0,critical_head_m=H0,remaining_gap_m=g,
       hydraulic_reopening_Q_threshold_m3_s=K*g**1.5,
       initial=before,doubled_discharge_pulse_fixed_geometry=pulse,return_to_baseline_fixed_geometry=after,
       grid_boundary_checks=rows,
       one_hour_maximum_possible_gap_repair_m_assuming_zero_bypass_melt=maxE*duration,
       required_gap_repair_m=deficit,
       minimum_repair_time_days_under_optimistic_flux_bound=deficit/maxE/86400,
       max_melt_rate_m_day=maxE*86400,
       elementary_constant_rate_removal_threshold_days=(H0+.08)/maxE/86400,
       conclusions=['A discharge pulse can reopen a fixed geometry hydraulically.',
       'Reopening alone does not restore baseline-flow contact.',
       'Under declared h/T upper bound, even zero bypass melt cannot repair this example in one hour.',
       'These follow from established ideal-control and Stefan equations, not new physics.'])
    p=Path(__file__).with_name('FALSIFIER_RESULT.json');p.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
