"""Independent randomized check of initial-state uncertainty in the contact solver.
No glacier outcomes. Comparison uses the same parameters/forcing in both runs.
"""
import json,random,sys
from pathlib import Path
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R/'kennicott_model'))
from column import simulate
rng=random.Random(833);records=[]
for case in range(100):
 n=4;m=5;C=1.5e6;dz=.2/n;dzi=2./m;Ci=917.*2100.;L=333700.
 D=[rng.uniform(-2,10)for _ in range(n)];I=[rng.uniform(-3,-.5)for _ in range(m)]
 Dh=[x+rng.uniform(0,3)for x in D];Ih=[min(0,x+rng.uniform(0,1))for x in I]
 Ts=[rng.uniform(-6,15)for _ in range(60)]
 kw=dict(surface_temperature_C=Ts,interval_seconds=600.,debris_depth_m=.2,debris_alpha_cells_m2_s=[1e-6]*n,debris_capacity_J_m3K=C,record_intervals=False)
 a=simulate(**kw,initial_debris_temperature_C=D,initial_ice_temperature_C=I)
 b=simulate(**kw,initial_debris_temperature_C=Dh,initial_ice_temperature_C=Ih)
 initial=sum(C*dz*(y-x)for x,y in zip(D,Dh))+sum(Ci*dzi*(y-x)for x,y in zip(I,Ih))
 delta=b['latent_heat_J_m2']-a['latent_heat_J_m2']
 order=min(y-x for x,y in zip(a['final_debris_temperature_C']+a['final_ice_temperature_C'],b['final_debris_temperature_C']+b['final_ice_temperature_C']))
 if delta < -1e-5 or delta>initial+1e-5 or order < -1e-10:raise AssertionError((case,delta,initial,order))
 records.append(dict(case=case,latent_difference_J_m2=delta,initial_energy_difference_J_m2=initial,min_final_temperature_order=order))
out=dict(scope='Synthetic numerical comparison, not new field validation or proof of historical novelty',cases=100,all_passed=True,records=records)
p=Path(__file__).with_name('CONTACT_COMPARISON_AUDIT.json')
with p.open('x')as f:json.dump(out,f,indent=2);f.write('\n')
print({'cases':100,'all_passed':True})
