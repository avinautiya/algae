"""Exact geometry for an extruded vertical step, with separate floor/face melt.

This nonconvex two-dimensional geometry is a deliberately restricted model:
planar face, horizontal stream floor, fixed upper surface, no overhang/collapse.
It does not turn a plan-view floor calculation into a full cliff simulation.
"""
from dataclasses import dataclass, asdict, replace
import math

RHO_I=917.
LATENT=334000.


@dataclass(frozen=True)
class SteppedCliff:
    upstream_edge_x_m: float
    ice_back_edge_x_m: float
    ice_bottom_z_m: float
    face_x_m: float
    floor_z_m: float
    crest_z_m: float
    extrusion_width_m: float

    def validate(self):
        for k,v in asdict(self).items():
            if isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(v):
                raise ValueError(k+' must be finite numeric')
        if not self.upstream_edge_x_m<self.face_x_m<self.ice_back_edge_x_m:
            raise ValueError('face must lie within horizontal domain')
        if not self.ice_bottom_z_m<self.floor_z_m<self.crest_z_m or self.extrusion_width_m<=0:
            raise ValueError('invalid vertical domain or width')
        return self

    @property
    def face_area_m2(self):return (self.crest_z_m-self.floor_z_m)*self.extrusion_width_m

    @property
    def floor_area_m2(self):return (self.face_x_m-self.upstream_edge_x_m)*self.extrusion_width_m

    def polygon(self):
        self.validate()
        return [(self.upstream_edge_x_m,self.ice_bottom_z_m),
                (self.ice_back_edge_x_m,self.ice_bottom_z_m),
                (self.ice_back_edge_x_m,self.crest_z_m),
                (self.face_x_m,self.crest_z_m),(self.face_x_m,self.floor_z_m),
                (self.upstream_edge_x_m,self.floor_z_m)]

    @property
    def ice_volume_m3(self):
        self.validate()
        return self.extrusion_width_m*((self.ice_back_edge_x_m-self.upstream_edge_x_m)*
            (self.floor_z_m-self.ice_bottom_z_m)+(self.ice_back_edge_x_m-self.face_x_m)*
            (self.crest_z_m-self.floor_z_m))


def advance_latent(cliff,face_latent_J,floor_latent_J,rho_ice=RHO_I,latent_heat=LATENT):
    """Apply net latent face energy first, then floor energy on the NEW floor.

    All areas are actual extruded surfaces. Two sequential removals are disjoint,
    avoiding double counting the corner. Energy outside the finite ice domain
    raises rather than silently clipping. Different stepping orders are a
    splitting uncertainty when both fluxes act continuously; refine dt.
    """
    cliff.validate()
    for k,v in {'face_latent_J':face_latent_J,'floor_latent_J':floor_latent_J,
                'rho_ice':rho_ice,'latent_heat':latent_heat}.items():
        if isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(v) or v<0:
            raise ValueError(k+' must be finite nonnegative numeric')
    if rho_ice==0 or latent_heat==0:raise ValueError('density and latent heat must be positive')
    dx=face_latent_J/(rho_ice*latent_heat*cliff.face_area_m2)
    after_face=replace(cliff,face_x_m=cliff.face_x_m+dx).validate()
    dz=floor_latent_J/(rho_ice*latent_heat*after_face.floor_area_m2)
    after=replace(after_face,floor_z_m=after_face.floor_z_m-dz).validate()
    removed_face=face_latent_J/(rho_ice*latent_heat)
    removed_floor=floor_latent_J/(rho_ice*latent_heat)
    measured=cliff.ice_volume_m3-after.ice_volume_m3
    return after,dict(face_retreat_m=dx,floor_lowering_m=dz,
                      face_melt_m3=removed_face,floor_melt_m3=removed_floor,
                      total_melt_m3=removed_face+removed_floor,
                      polygon_volume_loss_m3=measured,
                      volume_residual_m3=measured-removed_face-removed_floor,
                      phase_energy_J=face_latent_J+floor_latent_J,
                      note='prescribed uniform planar melt; not a fracture or survival model')


def face_partition(cliff,water_stage_m,debris_top_m):
    """Exclusive bare-air, bare-water, debris areas of the vertical face.

    Covered face is assigned to debris even when submerged. The caller must
    separately provide the correct submerged-debris boundary heat flux.
    """
    cliff.validate()
    for v in (water_stage_m,debris_top_m):
        if isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(v):raise ValueError('finite elevation required')
    bottom=cliff.floor_z_m;top=cliff.crest_z_m;width=cliff.extrusion_width_m
    cover=max(0.,min(top,debris_top_m)-bottom)
    bare_bottom=bottom+cover
    wet=max(0.,min(top,water_stage_m)-bare_bottom)
    air=top-bottom-cover-wet
    return dict(debris_area_m2=cover*width,bare_water_area_m2=wet*width,
                bare_air_area_m2=air*width,total_face_area_m2=(top-bottom)*width)
