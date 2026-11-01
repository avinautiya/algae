"""Exploratory Exner–Stefan parallel-route model, independent of Claude.

Critical broad-crested controls; quasi-steady shared head. A confined deposition
patch has mean sediment thickness S/[(1-porosity)WL]. Spatial sorting, seepage,
ice deformation and 2-D avulsion are absent. Grain export is an explicit
single-clast sliding closure, not a validated packed-bed transport law.
"""
from dataclasses import dataclass, asdict
import math
import numpy as np

@dataclass(frozen=True)
class Case:
    discharge: float=.1
    water_temperature: float=.5
    width: float=.5
    length: float=2.
    initial_toe_ice_z: float=-.08
    initial_bypass_z: float=0.
    initial_debris_depth: float=.01
    mantle_thickness: float=.3
    source_width: float=10.
    face_heat_W_m2: float=200.
    grain_sizes: tuple=(.02,.08,.20)
    source_fractions: tuple=(.4,.4,.2)
    friction: float=.05
    drag: float=1.05
    porosity: float=.4
    grain_density: float=2700.
    water_density: float=1000.
    ice_density: float=917.
    latent: float=334000.
    water_cp: float=4217.
    gravity: float=9.81
    h_water: float=200.
    debris_k: float=2.
    export_efficiency: float=1.
    extra_toe_heat_W_m2: float=0.
    allow_export: bool=True
    evolve_routing: bool=True
    allow_shield: bool=True
    footprint_shield: bool=True
    allow_ice_incision: bool=True
    allow_supply: bool=True

    def __post_init__(self):
        vals=asdict(self)
        flags={'allow_export','evolve_routing','allow_shield','allow_ice_incision','allow_supply','footprint_shield'}
        for k in flags:
            if not isinstance(vals[k],bool):raise ValueError('flags must be boolean')
        for k,v in vals.items():
            if k not in flags and isinstance(v,bool):raise ValueError('physical inputs cannot be boolean')
            if isinstance(v,(int,float)) and not isinstance(v,bool) and (not math.isfinite(v)):
                raise ValueError(k)
        positive=('discharge','width','length','water_density','ice_density','latent',
                  'water_cp','gravity','debris_k','grain_density','drag')
        if any(vals[k]<=0 for k in positive):raise ValueError('nonpositive physical input')
        if self.grain_density<=self.water_density:raise ValueError('resting nonbuoyant grains required')
        if not all(math.isfinite(v) and not isinstance(v,bool) for v in (*self.grain_sizes,*self.source_fractions)):raise ValueError('finite grain data required')
        if not 0<=self.porosity<1 or self.water_temperature<0 or self.friction<0 or self.h_water<0:
            raise ValueError('invalid material or forcing input')
        if self.export_efficiency>1:raise ValueError('export efficiency must not exceed one')
        if min(self.mantle_thickness,self.source_width,self.face_heat_W_m2,self.initial_debris_depth,self.export_efficiency,self.extra_toe_heat_W_m2)<0:
            raise ValueError('negative source, stock or export efficiency')
        if len(self.grain_sizes)!=len(self.source_fractions) or not self.grain_sizes or min(self.grain_sizes)<=0 or min(self.source_fractions)<0 or not math.isclose(sum(self.source_fractions),1.):
            raise ValueError('invalid grain distribution')


def partition(Q,z1,z2,K1,K2):
    """Exact prescribed-discharge split through two ideal critical controls."""
    Q=np.asarray(Q);z1=np.asarray(z1);z2=np.asarray(z2)
    low=np.minimum(z1,z2)
    high=np.maximum(z1,z2)+(Q/np.minimum(K1,K2))**(2/3)
    for _ in range(44):
        mid=(low+high)/2
        value=K1*np.maximum(mid-z1,0)**1.5+K2*np.maximum(mid-z2,0)**1.5
        low=np.where(value<Q,mid,low);high=np.where(value<Q,high,mid)
    eta=(low+high)/2
    return K1*np.maximum(eta-z1,0)**1.5,K2*np.maximum(eta-z2,0)**1.5,eta


def evolve(cases, days=30., dt=900., pulse=None, keep_history=True):
    if isinstance(dt,bool) or isinstance(days,bool):raise ValueError('numeric time required')
    if not cases or dt<=0 or days<=0 or not math.isfinite(dt+days):raise ValueError('invalid campaign')
    if len({len(c.grain_sizes) for c in cases})!=1:raise ValueError('matching grain array sizes required')
    n=len(cases)
    def a(name):return np.array([getattr(c,name) for c in cases],dtype=float)
    W,L,phi=a('width'),a('length'),a('porosity');area=W*L;solid_area=(1-phi)*area
    D=np.array([c.grain_sizes for c in cases]);fractions=np.array([c.source_fractions for c in cases])
    S=a('initial_debris_depth')[:,None]*solid_area[:,None]*fractions
    initial=S.copy();input_total=np.zeros_like(S);export_total=np.zeros_like(S)
    b1=a('initial_toe_ice_z');b2=a('initial_bypass_z')
    b10=b1.copy();b20=b2.copy();initial_d=S.sum(axis=1)/solid_area
    latent=a('ice_density')*a('latent')
    K=(2/3)**1.5*np.sqrt(a('gravity'))*W
    source=(1-phi)*a('mantle_thickness')*a('source_width')*a('face_heat_W_m2')/latent
    source*=a('allow_supply')
    Q0=a('discharge');melt1=np.zeros(n);melt2=np.zeros(n);heat_in=np.zeros(n);heat_out=np.zeros(n);external_toe_heat=np.zeros(n)
    wet_seconds=np.zeros(n);abandon=np.full(n,np.nan);history=[]
    steps=math.ceil(days*86400/dt)
    for step in range(steps):
        t=step*dt;h=min(dt,days*86400-t)
        Q=Q0.copy()
        if pulse is not None and pulse['start_day']*86400<=t<(pulse['start_day']+pulse['duration_day'])*86400:
            Q*=pulse['multiplier']
        d=S.sum(axis=1)/solid_area
        z1=b1+d;z2=b2
        z1=np.where(a('evolve_routing')>0,z1,b10+initial_d)
        z2=np.where(a('evolve_routing')>0,z2,b20)
        q1,q2,eta=partition(Q,z1,z2,K,K)
        # Ideal broad-crested control: critical depth is 2/3 of head.
        depth=(2/3)*np.maximum(eta-z1,0)
        velocity=np.divide(q1,W*depth,out=np.zeros(n),where=depth>0)
        Cw=a('water_density')*a('water_cp')
        hwater=a('h_water')
        # Parallel bare/covered patches. Cube footprint = solid volume / D.
        footprints=S/D
        total_footprint=np.sum(footprints,axis=1)
        cover=np.minimum(total_footprint/area,1.)
        fraction=np.divide(footprints,total_footprint[:,None],out=np.zeros_like(S),where=total_footprint[:,None]>0)
        covered_h=1/(1/np.maximum(hwater[:,None],1e-300)+D/a('debris_k')[:,None])
        patch_h=(1-cover)*hwater+cover*np.sum(fraction*covered_h,axis=1)
        uniform_h=1/(1/np.maximum(hwater,1e-300)+d/a('debris_k'))
        heff=np.where(a('footprint_shield')>0,patch_h,uniform_h)
        heff=np.where(a('allow_shield')>0,heff,hwater)
        heff=np.where(a('h_water')>0,heff,0.)
        def power(q,conductance):
            heat_capacity=Cw*q
            ntu=np.divide(conductance*area,heat_capacity,out=np.zeros(n),where=heat_capacity>0)
            return heat_capacity*a('water_temperature')*(-np.expm1(-ntu))
        P1=power(q1,heff)*a('allow_ice_incision');P2=power(q2,a('h_water'))*a('allow_ice_incision')
        external=a('extra_toe_heat_W_m2')*area*a('allow_ice_incision')
        Pwater1=P1.copy();P1+=external
        m1=P1/(latent*area);m2=P2/(latent*area)
        # Insulated-bed control sets P=0; source heat remains advected.
        # A fully/partly submerged cube. This is a sliding criterion only.
        sub=np.minimum(depth[:,None],D)
        normal=a('gravity')[:,None]*D**2*(a('grain_density')[:,None]*D-a('water_density')[:,None]*sub)
        threshold=np.divide(2*a('friction')[:,None]*normal,
                            a('water_density')[:,None]*a('drag')[:,None]*D*sub,
                            out=np.full_like(D,np.inf),where=sub>0)
        Uc=np.sqrt(np.maximum(threshold,0))
        speed=np.maximum(velocity[:,None]-Uc,0)
        active_fraction=np.minimum(1.,D/np.maximum(d[:,None],1e-300))
        rate=a('export_efficiency')[:,None]*a('allow_export')[:,None]*active_fraction*speed/L[:,None]
        incoming=h*source[:,None]*fractions
        newS=(S+incoming)/(1+h*rate)
        removed=S+incoming-newS
        input_total+=incoming;export_total+=removed;S=newS
        b1-=h*m1;b2-=h*m2
        melt1+=h*m1*area;melt2+=h*m2*area
        heat_in+=h*Cw*Q*a('water_temperature')
        heat_out+=h*(Cw*Q*a('water_temperature')-Pwater1-P2)
        external_toe_heat+=h*external
        wet=q1>1e-9*Q
        wet_seconds+=h*wet
        newly=(~wet)&np.isnan(abandon);abandon[newly]=t/86400
        if keep_history and (step%max(1,int(21600/dt))==0 or step==steps-1):
            history.append(dict(day=t/86400,toe_flow=q1.tolist(),debris_depth=d.tolist(),
                                toe_ice_z=(b1+h*m1).tolist(),bypass_z=(b2+h*m2).tolist(),
                                solid_stock_m3=(S-incoming+removed).tolist()))
    residual=S-initial-input_total+export_total
    return dict(cases=[asdict(c) for c in cases],days=days,dt=dt,
                toe_melt_m3=melt1.tolist(),bypass_melt_m3=melt2.tolist(),
                wet_days=(wet_seconds/86400).tolist(),first_dry_day=[float(x) if np.isfinite(x) else None for x in abandon],
                final_stock_m3=S.tolist(),input_solid_m3=input_total.tolist(),export_solid_m3=export_total.tolist(),
                final_toe_z=(b1+S.sum(axis=1)/solid_area).tolist(),final_bypass_z=b2.tolist(),
                max_solid_ledger_residual_m3=float(np.max(np.abs(residual))),
                energy_input_J=(heat_in+external_toe_heat).tolist(),external_toe_heat_J=external_toe_heat.tolist(),energy_advected_out_J=heat_out.tolist(),
                latent_J=(latent*(melt1+melt2)).tolist(),history=history,
                warning='exploratory reduced-model predictions; no independent environmental validation')
