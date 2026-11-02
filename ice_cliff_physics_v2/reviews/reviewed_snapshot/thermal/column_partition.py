"""Conservative wet/dry classification of normal debris temperature columns.

Profiles are pool averages at each normal cell, not face-height trajectories.
New material is imported into the dry class first by explicit convention. Only
material reclassified by a changing wet area mixes between wet/dry pools.
"""
from __future__ import annotations
from dataclasses import dataclass
import math


def _number(value,name,positive=False):
    if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or value<0 or (positive and value==0):
        raise ValueError('finite nonnegative '+name+' required')
    return float(value)


def _profile(value,name,allow_none=False):
    if value is None:
        if allow_none:return None
        raise ValueError(name+' required')
    if not isinstance(value,(list,tuple)) or not value:raise ValueError('nonempty '+name+' required')
    return tuple(_number(x,name+' temperature') for x in value)


@dataclass(frozen=True)
class ColumnPartition:
    capacity_J_m2_K: float
    wet_area_m2: float
    dry_area_m2: float
    wet_profile_C: tuple | None
    dry_profile_C: tuple | None

    def validate(self):
        _number(self.capacity_J_m2_K,'areal capacity',True)
        for area,profile,name in [(self.wet_area_m2,self.wet_profile_C,'wet'),(self.dry_area_m2,self.dry_profile_C,'dry')]:
            _number(area,name+' area')
            values=_profile(profile,name+' profile',allow_none=area==0)
            if area==0 and values is not None and any(x!=0 for x in values):
                raise ValueError('zero-area '+name+' pool must have None or zero profile')
        n=[len(p) for p in [self.wet_profile_C,self.dry_profile_C] if p is not None]
        if n and len(set(n))!=1:raise ValueError('wet/dry profile grids differ')
        return self

    @property
    def total_area_m2(self):return self.wet_area_m2+self.dry_area_m2

    @property
    def cell_count(self):
        return len(self.wet_profile_C) if self.wet_profile_C is not None else (len(self.dry_profile_C) if self.dry_profile_C is not None else 0)

    def cell_energies_J(self):
        n=self.cell_count
        if n==0:return tuple()
        cap=self.capacity_J_m2_K/n
        return tuple(cap*((self.wet_area_m2*self.wet_profile_C[i] if self.wet_profile_C is not None else 0.)+
                          (self.dry_area_m2*self.dry_profile_C[i] if self.dry_profile_C is not None else 0.)) for i in range(n))

    @property
    def stored_energy_J(self):return math.fsum(self.cell_energies_J())


@dataclass(frozen=True)
class RemapResult:
    partition: ColumnPartition
    imported_energy_J: float
    imported_energy_by_cell_J: tuple
    wetting_transfer_energy_J: float
    drying_transfer_energy_J: float
    energy_residual_J: float
    cell_energy_residuals_J: tuple
    expanded_area_m2: float
    wetting_transfer_area_m2: float
    drying_transfer_area_m2: float


def remap(partition,total_area_new_m2,wet_area_new_m2,incoming_profile_C):
    """Conserve each cell's extensive enthalpy through growth and reclassification.

    Total-area contraction is unsupported and refused. New area first receives
    the explicitly supplied incoming profile in the dry pool. Wetting then
    transfers representative dry-pool material; drying transfers representative
    wet-pool material. Pool-average profiles mix only when such material enters.
    A fixed wet/dry split preserves both profiles exactly.
    """
    if not isinstance(partition,ColumnPartition):raise TypeError('ColumnPartition required')
    p=partition.validate();total=_number(total_area_new_m2,'new total area');wet=_number(wet_area_new_m2,'new wet area')
    if total<p.total_area_m2:raise ValueError('area contraction requires an explicit outgoing-enthalpy operator')
    if wet>total:raise ValueError('wet area exceeds total area')
    incoming=_profile(incoming_profile_C,'incoming profile');n=len(incoming)
    if p.cell_count and p.cell_count!=n:raise ValueError('incoming and existing normal grids differ')
    cap=p.capacity_J_m2_K/n;growth=total-p.total_area_m2
    aw=p.wet_area_m2;ad=p.dry_area_m2
    pw=tuple(p.wet_profile_C) if aw>0 else None
    pd=tuple(p.dry_profile_C) if ad>0 else None
    if growth>0:
        enlarged=ad+growth
        pd=tuple(((ad*pd[i] if pd is not None else 0.)+growth*incoming[i])/enlarged for i in range(n))
        ad=enlarged
    imported=tuple(cap*growth*t for t in incoming)
    area_wetting=area_drying=energy_wetting=energy_drying=0.
    if wet>aw:
        area_wetting=wet-aw
        if ad<=0 or pd is None or area_wetting>ad+1e-14*max(1.,total):raise ArithmeticError('insufficient dry donor area')
        # Fixed-grid material enthalpy, with no source or sink at reclassification.
        energy_wetting=cap*area_wetting*math.fsum(pd)
        pw=tuple(((aw*pw[i] if pw is not None else 0.)+area_wetting*pd[i])/wet for i in range(n))
        aw=wet;ad=total-wet
        if ad==0:pd=None
    elif wet<aw:
        area_drying=aw-wet
        energy_drying=cap*area_drying*math.fsum(pw)
        enlarged=total-wet
        pd=tuple(((ad*pd[i] if pd is not None else 0.)+area_drying*pw[i])/enlarged for i in range(n))
        aw=wet;ad=enlarged
        if aw==0:pw=None
    else:
        ad=total-wet
        if aw==0:pw=None
        if ad==0:pd=None
    updated=ColumnPartition(p.capacity_J_m2_K,aw,ad,pw,pd).validate()
    before=p.cell_energies_J() if p.cell_count else tuple([0.]*n)
    after=updated.cell_energies_J() if updated.cell_count else tuple([0.]*n)
    residual=tuple(b-a-e for a,b,e in zip(before,after,imported))
    energies=list(before)+list(after)+list(imported)+[energy_wetting,energy_drying]
    if not all(math.isfinite(e) for e in energies):raise ValueError('derived material enthalpy overflow')
    scale=max(1.,max((abs(e) for e in energies),default=0.))
    if max((abs(e) for e in residual),default=0.)>1e-12*scale:raise ArithmeticError('per-cell material enthalpy conservation failure')
    return RemapResult(updated,math.fsum(imported),imported,energy_wetting,energy_drying,
        math.fsum(residual),residual,growth,area_wetting,area_drying)
