"""Independent analytic comparison of the already saved slab mini-campaign."""
from pathlib import Path
import hashlib,json,math
from robin_reference import RobinSlab
ROOT=Path(__file__).resolve().parent

def main():
    saved=json.loads((ROOT/'SLAB_COMPARISON.json').read_text());rows=[];worst=0.
    for case in saved['main_cases']:
        G=case['Bi_nominal']*20
        a=RobinSlab.from_contact(90000,20,400,G,.5,0,256)
        b=RobinSlab.from_contact(90000,20,400,G,.5,0,512)
        t=case['time_s'];ra=a.response(t);rb=b.response(t)
        conditioning=rb['roundoff_cancellation_scale_J'];heat=rb['ice_heat_J']
        meaningful=abs(heat)>100*conditioning
        scale=max(1.,abs(rb['boundary_heat_J']),abs(rb['stored_energy_J']))
        error=abs(rb['energy_residual_J'])/scale;worst=max(worst,error)
        rows.append(dict(Bi_nominal=case['Bi_nominal'],Robin_B=b.B,time_over_diffusion_time=t/b.tau,
            time_s=t,analytic_ice_heat_J=heat,analytic_256_512_difference_J=ra['ice_heat_J']-heat,
            analytic_energy_residual_J=rb['energy_residual_J'],analytic_roundoff_cancellation_scale_J=conditioning,
            reference_heat_large_enough_for_relative_comparison=meaningful,
            finite_volume_ice_heat_J=case['slab_covered_heat_J'],lumped_ice_heat_J=case['lumped_covered_heat_J'],
            finite_volume_relative_error=None if not meaningful else (case['slab_covered_heat_J']-heat)/heat,
            lumped_relative_error=None if not meaningful else (case['lumped_covered_heat_J']-heat)/heat,
            finite_volume_absolute_error_J=case['slab_covered_heat_J']-heat))
    result=dict(classification='classical analytic component verification; no new theorem or field accuracy claim',
        zero_initial_temperature=True,constant_coefficients=True,cases=rows,
        max_normalized_analytic_energy_residual=worst,
        cancelled_early_relative_errors='not reported when reference heat <=100 times arithmetic conditioning scale; this scale is not a certified interval bound',
        saved_campaign_sha256=hashlib.sha256((ROOT/'SLAB_COMPARISON.json').read_bytes()).hexdigest(),
        code_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'robin_reference.py',ROOT/'test_robin_reference.py',Path(__file__)]})
    (ROOT/'ROBIN_REFERENCE_VERIFICATION.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
