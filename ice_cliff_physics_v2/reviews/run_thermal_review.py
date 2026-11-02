"""Independent bounded coupling review. No edits to integrated source."""
from pathlib import Path
import copy,hashlib,json,random,sys,traceback
ROOT=Path('/private/tmp/cliff-upgrade-20261101/package')
OUT=Path('/private/tmp/cliff-upgrade-20261101/reviews')
sys.path.insert(0,str(ROOT))
from model import forecast

def main():
    base=json.loads((ROOT/'example.json').read_text());rng=random.Random(20261101);cases=[]
    for i in range(20):
        c=copy.deepcopy(base);p=c['physics']
        p.update(initial_surface_cold_J_m2=0.,new_surface_cold_J_m2=0.,ice_subcooling_K=0.,
                 face_cover_fraction=rng.uniform(.01,.1),debris_capacity_J_m2_K=rng.uniform(20000,100000),
                 debris_conductance_W_m2_K=rng.uniform(5,100),initial_debris_temperature_C=rng.uniform(0,8),
                 new_debris_temperature_C=rng.uniform(0,4))
        q=rng.uniform(.02,.1);tw=rng.uniform(0,1);flux=rng.uniform(0,250)
        c['periods']=[dict(duration_s=3600.,inflow_m3_s=q,water_superheat_K=tw,
            face_atmospheric_net_W_m2=flux,floor_atmospheric_net_W_m2=flux,bypass_atmospheric_net_W_m2=flux),
            dict(duration_s=3600.,inflow_m3_s=q*.9,water_superheat_K=tw*.7,
            face_atmospheric_net_W_m2=flux*.6,floor_atmospheric_net_W_m2=flux*.6,bypass_atmospheric_net_W_m2=flux*.6)]
        c['step_s']=120.;c['removal_time_s']=3600. if i%2 else None
        c['upgrades']=dict(wet_debris='lumped_temperate' if i%4==0 else 'resolved_temperate',
            face_cover_location='lower' if i%2 else 'upper',covered_heat_policy='prescribed_all_cover' if i%2 else 'dry_only',
            debris_cells=8,debris_substep_fraction=.01)
        if i%3==0:c['upgrades'].update(face_profile_bands=8,ledge_heat_policy='adiabatic')
        if i%2==0:
            floor=c['geometry']['floor_z_m'];height=c['geometry']['crest_z_m']-floor
            lo=floor+.01;hi=lo+p['face_cover_fraction']*height
            c['upgrades'].update(fixed_cover_interval_m=[lo,hi],cover_interval_provenance='SYNTHETIC fixed elevation mask; independently prescribed from geometry before run')
        c['provenance']='SYNTHETIC independent thermal integration review; seeded warm inputs, not field data'
        cases.append(dict(id=f'T{i+1:02}',configuration=c))
    sources=[ROOT/'model.py',ROOT/'face_thermal.py',ROOT/'thermal/wet_debris.py',ROOT/'thermal/resolved_slab.py',ROOT/'thermal/column_partition.py']
    source_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    frozen=dict(seed=20261101,n_cases=20,cases=cases,source_sha256=source_hashes)
    (OUT/'thermal_review_cases_frozen.json').write_text(json.dumps(frozen,indent=2)+'\n')
    results=[];worst=0.;fail=[]
    for case in cases:
        c=case['configuration']
        try:
            r=forecast(c);energy_scale=max(1.,abs(r['totals']['atmospheric_J'])+r['totals']['water_in_J'])
            normalized=abs(r['checks']['global_energy_residual_J'])/energy_scale;worst=max(worst,normalized)
            reservoir_error=[];wetsteps=0
            for row in r['series']:
                therm=row['wet_debris_thermal']
                if row['wet_covered_face_area_m2']>0:wetsteps+=1
                if therm is not None and 'partition' in therm:
                    pp=therm['partition'];n=len(pp['wet_profile_C']) if pp['wet_profile_C'] is not None else len(pp['dry_profile_C'])
                    e=pp['capacity_J_m2_K']/n*((pp['wet_area_m2']*sum(pp['wet_profile_C']) if pp['wet_profile_C'] is not None else 0.)+(pp['dry_area_m2']*sum(pp['dry_profile_C']) if pp['dry_profile_C'] is not None else 0.))
                    state=row['port_states']['face_debris'];expected_area=state['area_m2']
                    # Row thermal is BEFORE the geometry's area expansion; its
                    # post-remap scalar reservoir may include declared imports.
                    reservoir_error.append(dict(thermal_area_m2=pp['wet_area_m2']+pp['dry_area_m2'],post_geometry_area_m2=expected_area,
                        thermal_stored_J=e,post_geometry_stored_J=state['debris_energy_j']))
            results.append(dict(id=case['id'],status='completed',mode=c['upgrades']['wet_debris'],
                mask='fixed' if 'fixed_cover_interval_m' in c['upgrades'] else 'relative_fraction',
                heat_policy=c['upgrades']['covered_heat_policy'],checks=r['checks'],normalized_energy_residual=normalized,
                steps=len(r['series']),wet_debris_steps=wetsteps,
                face_melt_m3=r['totals']['face_melt_m3'],max_phase_water_ratio=max(x['meltwater_to_inflow_ratio'] or 0 for x in r['series']),
                reservoir_samples=reservoir_error[:1]+reservoir_error[-1:]))
        except ValueError as exc:
            msg=str(exc)
            expected=any(x in msg for x in ['more head','cold','refreez','overtop','domain exhausted','scope','hydraulic feedback','local route meltwater feedback','drainage capacity'])
            results.append(dict(id=case['id'],status='expected_scope_refusal' if expected else 'unexpected_refusal',reason=msg))
            if not expected:fail.append(case['id'])
        except Exception as exc:
            results.append(dict(id=case['id'],status='failure',reason=repr(exc),traceback=traceback.format_exc()));fail.append(case['id'])
    # Separate explicit refusal/zero-inflow contract checks, not counted as the
    # twenty warm development configurations.
    guard=[]
    c=copy.deepcopy(cases[1]['configuration']);c['physics']['initial_surface_cold_J_m2']=1
    try:forecast(c);guard.append(dict(test='cold_face',expected='refuse',actual='accepted'))
    except ValueError as e:guard.append(dict(test='cold_face',expected='refuse',actual='refused',reason=str(e)))
    for allow in [False,True]:
        c=copy.deepcopy(cases[2]['configuration']);c['upgrades'].pop('face_profile_bands',None);c['upgrades'].pop('ledge_heat_policy',None)
        for row in c['periods']:row['inflow_m3_s']=0.;row['face_atmospheric_net_W_m2']=100.;row['floor_atmospheric_net_W_m2']=100.;row['bypass_atmospheric_net_W_m2']=100.
        if allow:c['upgrades']['meltwater_policy']=dict(mode='downstream_free_drainage',max_generated_discharge_m3_s=.001,provenance='SYNTHETIC separated downstream meltwater drain, no mixing with imposed warm stream, no backwater, capacity .001 m3/s')
        try:
            r=forecast(c);guard.append(dict(test='zero_inflow_with_contract' if allow else 'zero_inflow_without_contract',expected='accept' if allow else 'refuse',actual='accepted',checks=r['checks']))
        except ValueError as e:guard.append(dict(test='zero_inflow_with_contract' if allow else 'zero_inflow_without_contract',expected='accept' if allow else 'refuse',actual='refused',reason=str(e)))
    now={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    out=dict(classification='independent synthetic integration review, not field validation',seed=20261101,source_sha256_before=source_hashes,source_sha256_after=now,
        source_changed_during_review=source_hashes!=now,results=results,guard_results=guard,
        max_normalized_energy_residual=worst,unexpected_failure_ids=fail)
    (OUT/'thermal_review.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(status_counts={k:sum(x['status']==k for x in results) for k in set(x['status'] for x in results)},worst=worst,failures=fail,guards=guard,source_changed=source_hashes!=now),indent=2),flush=True)

if __name__=='__main__':main()
