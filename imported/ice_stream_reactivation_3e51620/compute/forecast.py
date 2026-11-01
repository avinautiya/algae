"""CLI for an explicitly conditional glacier-stream reactivation forecast."""
import argparse,json
from pathlib import Path
from model import Case,evolve
from rescue_certificate import thermal_lockout

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',type=Path,help='JSON object with parameters and optional state/flood objects')
    parser.add_argument('--simulate-days',type=float)
    parser.add_argument('--dt',type=float,default=225.)
    args=parser.parse_args();spec=json.loads(args.input.read_text())
    c=Case(**spec.get('parameters',{}))
    if args.simulate_days is not None:
        result=evolve([c],days=args.simulate_days,dt=args.dt,pulse=spec.get('pulse'))
    else:
        result=thermal_lockout(c,**spec['state'],**spec['flood'])
    def finite(x):
        if isinstance(x,float):return x if __import__('math').isfinite(x) else None
        if isinstance(x,list):return [finite(a) for a in x]
        if isinstance(x,dict):return {k:finite(v) for k,v in x.items()}
        return x
    print(json.dumps(finite(result),indent=2,allow_nan=False))
if __name__=='__main__':main()
