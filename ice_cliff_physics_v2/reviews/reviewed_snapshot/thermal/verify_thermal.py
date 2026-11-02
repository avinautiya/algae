"""Reproducible implementation/stress verification, not environmental evidence."""
import hashlib,json,random,subprocess,sys
from pathlib import Path
from wet_debris import ThermalState,ThermalParameters,advance
ROOT=Path(__file__).resolve().parent

def main():
    test=subprocess.run([sys.executable,'-m','unittest','discover','-s',str(ROOT),'-p','test_*.py'],capture_output=True,text=True)
    rng=random.Random(20261101);worst=0.;fail=[]
    for i in range(5000):
        C=10**rng.uniform(-12,12);K=10**rng.uniform(-9,5);gb=10**rng.uniform(-9,5);gd=10**rng.uniform(-9,5)
        F=10**rng.uniform(-12,10);dt=10**rng.uniform(-6,8);t0=10**rng.uniform(-6,2);tw=10**rng.uniform(-6,1);air=10**rng.uniform(-9,5)
        try:
            r=advance(ThermalState(t0),ThermalParameters(C,K,gb,gd),dt,tw,F,air)
            scale=max(1.,C*t0,r.water_inlet_heat_J,r.water_outlet_heat_J,air*dt,r.total_ice_heat_J)
            error=abs(r.global_energy_residual_J)/scale;worst=max(worst,error)
            if error>1e-11 or r.water_outlet_heat_J<0:fail.append(dict(case=i,error=error))
        except Exception as e:fail.append(dict(case=i,error=repr(e)))
    result=dict(classification='synthetic implementation and numerical conditioning tests; not field accuracy',
        unit_test_returncode=test.returncode,unit_test_output=test.stdout+test.stderr,
        stress_seed=20261101,stress_cases=5000,stress_failures=fail,stress_max_normalized_energy_residual=worst,
        code_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(ROOT.glob('*.py'))})
    (ROOT/'VERIFICATION.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    if test.returncode or fail:raise SystemExit(1)

if __name__=='__main__':main()
