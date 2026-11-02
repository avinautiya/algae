"""Descriptive independent-depth check (no model, no fitting): water-level-cut elevation minus time-interpolated surveyed floor vs author pressure-derived daily stage.
The level file's header says 'P4 water level cuts' while the file is named P5, and its timestamps are a stated hypothesis (only 30 July true); both are reported, not resolved.
python -m ice_integration.diagnosis.clearance_check OUT.json"""
import csv, json, sys, os
from datetime import datetime
R=os.path.join(os.path.dirname(__file__),"..","..","ice_swiss_field_validation","original_mirror")
def tm(s):
    s=s.strip()
    for f in("%Y-%m-%dT%H:%M:%S","%Y-%m-%d %H:%M:%S","%Y-%m-%d"):
        try: return datetime.strptime(s[:19] if "T" in s or " " in s else s,f)
        except ValueError: pass
    if s.count(":")==1 and s.endswith(":0"): return datetime.strptime(s+"0","%Y-%m-%dT%H:%M:%S")
    raise ValueError(s)
def rd(p):
    with open(p,newline="") as f: return list(csv.reader(f,delimiter=";" if ";" in f.readline() else ","))
def main(out):
    floor=[(tm(r[0]),float(r[1]),float(r[2])) for r in rd(os.path.join(R,"data/Channel_elevation/P5_channel_elevation.csv"))[1:]]
    lev=[(tm(r[0]),float(r[1]),float(r[2])) for r in rd(os.path.join(R,"data/Channel_elevation/P5_levels_elevation.csv"))[1:]]
    st={}
    for r in rd(os.path.join(R,"products/stage5_daily.csv")): 
        try: st[tm(r[0]).date()]=float(r[1])
        except Exception: pass
    def zb(t):
        for (t0,z0,s0),(t1,z1,s1) in zip(floor,floor[1:]):
            if t0<=t<=t1: f=(t-t0).total_seconds()/(t1-t0).total_seconds(); return z0+f*(z1-z0),(s0**2+s1**2)**.5
        return None
    rows=[]
    for t,eta,s in lev:
        z=zb(t)
        if z is None: continue
        h=eta-z[0]; sig=(s**2+z[1]**2)**.5; p=st.get(t.date())
        rows.append(dict(t=str(t),eta_m=eta,floor_interp_m=round(z[0],3),depth_independent_m=round(h,3),sigma_m=round(sig,3),pressure_stage_m=p,pressure_minus_independent_m=None if p is None else round(p-h,3)))
    json.dump(dict(rows=rows,caveat="Floor is linearly interpolated between surveys spaced 5-9 days apart during active incision (floor can drop ~0.4 m/day); level cut is at the channel edge and its header names P4. Indicative only; not a clearance measurement."),open(out,"w"),indent=1)
    for r in rows: print(r)
if __name__=="__main__": main(sys.argv[1])
