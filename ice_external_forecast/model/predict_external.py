"""Frozen external thermal adapter; no observations, fitting or built-in site geometry.

Finite advective heat under imposed Q and piecewise-uniform rectangular sections.
Point-control outputs are local flux operators, not finite melt volumes.
"""
import argparse
from dataclasses import dataclass
import datetime as dt
import hashlib
import json
import math
from pathlib import Path

RHO_W, CP_W, NU, PR, K = 1000.,4210.,1.8e-6,13.5,.56
RHO_I, LATENT, CP_I, G = 917.,334000.,2100.,9.81
RATE, H_FIXED, DB_SCALE = 3.4944254768079555e-6,6431.4136285363775,4.3956164544491925
PHYSICAL_MODELS=('constant_rate_fitted_heat_limited','constant_h_fitted','db04_unfitted',
                 'db03_unfitted','gnielinski_smooth_unfitted','db04_scale_fitted')
REFERENCE='constant_rate_fitted_uncapped_reference'

class ApplicabilityRefusal(ValueError):
    pass

def number(x,name,positive=False,nonnegative=False):
    if isinstance(x,bool) or not isinstance(x,(int,float)) or not math.isfinite(x):
        raise ApplicabilityRefusal(f'{name} must be finite numeric')
    if positive and x<=0 or nonnegative and x<0:
        raise ApplicabilityRefusal(f'{name} outside supported sign range')
    return float(x)

def text(x,name):
    if not isinstance(x,str) or not x.strip():raise ApplicabilityRefusal(f'{name} source/provenance required')
    return x

@dataclass(frozen=True)
class Section:
    name:str
    width_m:float
    depth_m:float
    length_m:float
    point_control:bool
    source:str
    geometry_basis:str
    potential_drop_m:float=None

    def __post_init__(self):
        text(self.name,'section name');text(self.source,'section source')
        number(self.width_m,'width',positive=True);number(self.depth_m,'depth',positive=True)
        number(self.length_m,'length',nonnegative=True)
        if not isinstance(self.point_control,bool) or (self.point_control and self.length_m!=0) or (not self.point_control and self.length_m<=0):
            raise ApplicabilityRefusal('section must be a positive reach or zero-area point control')
        if self.potential_drop_m is not None:number(self.potential_drop_m,'potential drop',nonnegative=True)

    @classmethod
    def from_input(cls,d,q):
        q=number(q,'Q',positive=True)
        name=text(d.get('name'),'section name');source=text(d.get('source'),'geometry')
        basis=d.get('geometry_basis')
        if basis not in ('measured','source_assumption','synthetic_fixture'):
            raise ApplicabilityRefusal('geometry_basis must be explicitly measured/source_assumption/synthetic_fixture')
        width=number(d.get('width_m'),'width',positive=True)
        kind=d.get('kind','uniform_reach');point=kind=='point_control'
        if kind not in ('point_control','uniform_reach'):raise ApplicabilityRefusal('unsupported section kind')
        length=number(d.get('length_m'),'along-flow section length',nonnegative=True)
        if point and length!=0 or not point and length<=0:raise ApplicabilityRefusal('point controls have zero area; uniform reaches need positive length')
        if d.get('depth_mode')=='critical_control':
            if not point or d.get('critical_control_verified') is not True:
                raise ApplicabilityRefusal('critical depth allowed only at a verified point control')
            text(d.get('critical_control_source'),'critical-control location/width evidence')
            unit_discharge=number(q/width,'critical unit discharge',positive=True)
            depth=unit_discharge**(2/3)/G**(1/3)
        else:
            depth=number(d.get('depth_m'),'section water depth',positive=True)
        drop=d.get('potential_drop_m')
        if drop is not None:
            drop=number(drop,'potential drop scale',nonnegative=True)
            text(d.get('potential_drop_source'),'potential-drop evidence')
        return cls(name,width,depth,length,point,source,basis,drop)

    def geometry(self,q):
        area=number(self.width_m*self.depth_m,'derived area',positive=True)
        perimeter=number(self.width_m+2*self.depth_m,'derived wetted perimeter',positive=True)
        diameter=number(4*(area/perimeter),'derived hydraulic diameter',positive=True)
        velocity=number(q/area,'derived velocity',positive=True)
        re=number(velocity*diameter/NU,'derived Reynolds',positive=True)
        # Common gate keeps identical sections across all physical comparators.
        if not 10000<=re<=5e6:raise ApplicabilityRefusal(f'{self.name}: Re outside common comparator domain')
        return dict(area_m2=area,wetted_perimeter_m=perimeter,hydraulic_diameter_m=diameter,
                    velocity_m_s=velocity,Re=re,Fr=velocity/math.sqrt(G*self.depth_m))

def transfer_coefficient(method,re,diameter):
    if method=='constant_h_fitted':return H_FIXED
    if method in ('db04_unfitted','db03_unfitted','db04_scale_fitted'):
        exponent=.3 if method=='db03_unfitted' else .4
        scale=DB_SCALE if method=='db04_scale_fitted' else 1.
        return scale*K/diameter*.023*re**.8*PR**exponent
    if method=='gnielinski_smooth_unfitted':
        f=(1.82*math.log10(re)-1.64)**-2
        nu=(f/8)*(re-1000)*PR/(1+12.7*math.sqrt(f/8)*(PR**(2/3)-1))
        return K/diameter*nu
    raise ApplicabilityRefusal('unsupported heat-transfer family')

def predict_interval(sections,q_m3_s,inlet_temperature_C,duration_s,ice_temperature_C,method):
    q=number(q_m3_s,'Q',positive=True);tin=number(inlet_temperature_C,'inlet temperature',nonnegative=True)
    duration=number(duration_s,'duration',positive=True);ti=number(ice_temperature_C,'ice temperature')
    if ti>0 or ti<=-273.15:raise ApplicabilityRefusal('ice temperature outside supported physical phase range')
    if method not in PHYSICAL_MODELS:raise ApplicabilityRefusal('physical model name required')
    sensible_per_kg=CP_I*(-ti);cost=LATENT+sensible_per_kg
    capacity=number(RHO_W*CP_W*q,'derived water heat capacity rate',positive=True)
    temperature=tin;total_removed=total_volume=0.;rows=[]
    for section in sections:
        geom=section.geometry(q);p=geom['wetted_perimeter_m'];length=section.length_m
        if method=='constant_rate_fitted_heat_limited':
            demand_flux=RHO_I*cost*RATE
            requested_power=number(demand_flux*p*length,'requested section power',nonnegative=True)
            removed_power=min(capacity*temperature,requested_power)
            out=max(0.,temperature-removed_power/capacity)
            inlet_flux=demand_flux if temperature>0 else 0.
            outlet_flux=demand_flux if out>0 else 0.
            h=None;ntu=None
        else:
            h=transfer_coefficient(method,geom['Re'],geom['hydraulic_diameter_m'])
            ntu=number(h*p*length/capacity,'derived NTU',nonnegative=True)
            removed_power=capacity*temperature*(-math.expm1(-ntu))
            out=temperature*math.exp(-ntu)
            inlet_flux=h*temperature;outlet_flux=h*out;requested_power=None
        removed=removed_power*duration;bed=removed*section.width_m/p;walls=removed-bed
        volume=removed/(RHO_I*cost);bed_volume=bed/(RHO_I*cost);wall_volume=walls/(RHO_I*cost)
        mean_lowering=None if section.point_control else bed_volume/(section.width_m*length)
        # A zero-area point has a local flux/recession but zero total heat/volume.
        local_lowering=inlet_flux*duration/(RHO_I*cost) if section.point_control else None
        potential=None if section.potential_drop_m is None else RHO_W*G*q*section.potential_drop_m
        row=dict(name=section.name,kind='point_control' if section.point_control else 'uniform_reach',
            geometry_source=section.source,geometry_basis=section.geometry_basis,width_m=section.width_m,
            depth_m=section.depth_m,length_m=length,**geom,inlet_temperature_C=temperature,outlet_temperature_C=out,
            transfer_h_W_m2_K=h,NTU=ntu,delivered_heat_W=removed_power,delivered_heat_J=removed,
            requested_heat_W=requested_power,bed_heat_J=bed,sidewall_heat_J=walls,
            inlet_interface_heat_flux_W_m2=inlet_flux,outlet_interface_heat_flux_W_m2=outlet_flux,
            mean_bed_lowering_m=mean_lowering,point_bed_lowering_m=local_lowering,
            removed_ice_volume_m3=volume,bed_ice_volume_m3=bed_volume,sidewall_ice_volume_m3=wall_volume,
            latent_heat_J=volume*RHO_I*LATENT,fresh_ice_sensible_heat_J=volume*RHO_I*sensible_per_kg,
            potential_energy_power_scale_W=potential,
            potential_scale_over_initial_sensible_power=None if potential is None or tin==0 else potential/(capacity*tin))
        rows.append(row);temperature=out;total_removed+=removed;total_volume+=volume
    initial=capacity*tin*duration;final=capacity*temperature*duration
    latent=total_volume*RHO_I*LATENT;sensible=total_volume*RHO_I*sensible_per_kg
    return dict(model=method,sections=rows,inlet_heat_J=initial,outlet_heat_J=final,
        removed_water_heat_J=total_removed,latent_heat_J=latent,fresh_ice_sensible_heat_J=sensible,
        energy_residual_J=initial-final-total_removed,interface_energy_residual_J=total_removed-latent-sensible,
        removed_ice_volume_m3=total_volume,produced_meltwater_kg=total_volume*RHO_I,
        produced_meltwater_to_imposed_input_water_mass_fraction=total_volume*RHO_I/(q*duration*RHO_W),
        meltwater_feedback='reported inventory; not fed back into imposed hydraulic Q or water heat capacity',
        omitted_sources=['atmospheric_exchange','frictional_dissipation'],closed_physical_component_ledger=True)

def uncapped_reference(sections,q,tin,duration,ti):
    cost=LATENT+CP_I*(-ti);rows=[];requested=volume=0.
    for s in sections:
        geom=s.geometry(q);v=RATE*duration*geom['wetted_perimeter_m']*s.length_m
        energy=v*RHO_I*cost;volume+=v;requested+=energy
        rows.append(dict(name=s.name,normal_or_vertical_bed_lowering_m=RATE*duration,
                         requested_ice_volume_m3=v,requested_equivalent_heat_J=energy))
    available=RHO_W*CP_W*q*tin*duration
    return dict(model=REFERENCE,sections=rows,requested_ice_volume_m3=volume,
        requested_equivalent_heat_J=requested,available_inlet_heat_J=available,
        exceeds_inlet_heat_supply=requested>available,closed_physical_component_ledger=False,
        scope='uncapped trained statistical reference; no physical outlet temperature or delivered-heat ledger')

def parse_time(s,clock):
    text(s,'timestamp')
    try:t=dt.datetime.fromisoformat(s.replace('Z','+00:00'))
    except ValueError as exc:raise ApplicabilityRefusal('ISO timestamp required') from exc
    if clock=='UTC':
        if t.tzinfo is None or t.utcoffset()!=dt.timedelta(0):raise ApplicabilityRefusal('UTC requires explicit +00:00 or Z')
    elif clock=='declared_source_naive':
        if t.tzinfo is not None:raise ApplicabilityRefusal('source-naive packet requires consistently naive times')
    else:raise ApplicabilityRefusal('explicit supported clock convention required')
    return t

def predict_packet(packet):
    purpose=packet.get('purpose')
    if purpose not in ('synthetic_verification','conditional_external_hindcast'):
        raise ApplicabilityRefusal('explicit prediction purpose required')
    alignment=packet.get('alignment',{})
    if alignment.get('path_alignment_verified') is not True or alignment.get('geometry_independent_of_bed_target') is not True:
        raise ApplicabilityRefusal('inlet-to-target path and independent geometry matching gate incomplete')
    text(alignment.get('source'),'path alignment evidence')
    if not isinstance(packet.get('records'),list) or not packet['records']:raise ApplicabilityRefusal('supported forcing intervals required')
    outputs=[];previous=None
    for record in packet['records']:
        a=parse_time(record.get('start'),packet.get('clock_convention'));b=parse_time(record.get('end'),packet.get('clock_convention'))
        if b<=a or previous is not None and a!=previous:raise ApplicabilityRefusal('positive contiguous forcing intervals required; no gaps/overlaps/extrapolation')
        previous=b;duration=(b-a).total_seconds();q=number(record.get('Q_m3_s'),'Q',positive=True)
        tin=number(record.get('inlet_temperature_C'),'inlet temperature',nonnegative=True)
        ti=number(record.get('ice_temperature_C'),'explicit ice temperature')
        if ti>0 or ti<=-273.15:raise ApplicabilityRefusal('ice temperature outside supported physical phase range')
        text(record.get('forcing_source'),'forcing')
        if record.get('forcing_basis') not in ('interval_observed','sourced_interval_scenario','synthetic_fixture'):
            raise ApplicabilityRefusal('point readings are not complete interval forcing')
        if purpose=='conditional_external_hindcast' and record.get('forcing_basis')=='synthetic_fixture':
            raise ApplicabilityRefusal('synthetic forcing cannot be relabelled external field input')
        raw_sections=record.get('sections')
        if not isinstance(raw_sections,list) or not raw_sections:raise ApplicabilityRefusal('complete independently sourced thermal route sections required')
        sections=[Section.from_input(d,q) for d in raw_sections]
        if len({s.name for s in sections})!=len(sections):raise ApplicabilityRefusal('section names must be unique for target matching')
        if purpose=='conditional_external_hindcast' and any(s.geometry_basis=='synthetic_fixture' for s in sections):
            raise ApplicabilityRefusal('synthetic geometry cannot be relabelled external field input')
        models={name:predict_interval(sections,q,tin,duration,ti,name) for name in PHYSICAL_MODELS}
        models[REFERENCE]=uncapped_reference(sections,q,tin,duration,ti)
        outputs.append(dict(start=record['start'],end=record['end'],forcing_source=record['forcing_source'],
            forcing_basis=record['forcing_basis'],models=models))
    return dict(purpose=purpose,scope='imposed-flow thermal component; no target observations or scoring',
                clock_convention=packet['clock_convention'],alignment=alignment,intervals=outputs)

def main():
    p=argparse.ArgumentParser();p.add_argument('packet',type=Path);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    try:
        freeze_root=Path(__file__).resolve().parent
        frozen=json.loads((freeze_root/'EXTERNAL_MODEL_FREEZE_v2.json').read_text())
        expected=frozen['parameters']
        if expected!=dict(db04_scale=DB_SCALE,constant_h_W_m2_K=H_FIXED,constant_rate_m_s=RATE):
            raise ApplicabilityRefusal('constants differ from external parameter freeze')
        for name,key in [('EXTERNAL_PROTOCOL_FROZEN.md','protocol_sha256'),('BOUNDARY_POLICY_ADDENDUM.md','boundary_policy_sha256')]:
            if hashlib.sha256((freeze_root/name).read_bytes()).hexdigest()!=frozen[key]:
                raise ApplicabilityRefusal('protocol or boundary policy changed after freeze')
        result=predict_packet(json.loads(a.packet.read_text()))
        result['model_freeze']=frozen
    except (ApplicabilityRefusal,KeyError,TypeError,OSError,json.JSONDecodeError) as exc:
        print(json.dumps(dict(status='REFUSED',reason=str(exc))));raise SystemExit(2)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    result['input_sha256']=hashlib.sha256(a.packet.read_bytes()).hexdigest()
    a.out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(status='PREDICTED_COMPONENT_ONLY',purpose=result['purpose'],intervals=len(result['intervals']))))

if __name__=='__main__':main()
