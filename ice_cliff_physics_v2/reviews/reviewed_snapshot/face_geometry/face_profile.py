"""Conservative extruded height-band geometry, including geometric undercuts.

This is a prescribed latent-energy geometry operator, not a fracture, flow,
material-mask advection, or full moving-interface heat solver. Horizontal ledges
created by nonuniform retreat are reported, not silently treated as vertical face.
"""
from dataclasses import asdict, dataclass, replace
import math

RHO_ICE = 917.
LATENT = 334000.


def finite(v, name, positive=False, nonnegative=False):
    if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v):
        raise ValueError(name+' must be finite numeric')
    if positive and v <= 0 or nonnegative and v < 0:
        raise ValueError(name+' outside sign range')
    return float(v)


@dataclass(frozen=True)
class Band:
    z_low_m: float
    z_high_m: float
    front_x_m: float

    def __post_init__(self):
        for k, v in asdict(self).items():
            finite(v, k)
        if self.z_high_m <= self.z_low_m:
            raise ValueError('band needs positive height')


@dataclass(frozen=True)
class Profile:
    upstream_edge_x_m: float
    ice_back_edge_x_m: float
    ice_bottom_z_m: float
    crest_z_m: float
    extrusion_width_m: float
    bands: tuple

    def __post_init__(self):
        for k in ['upstream_edge_x_m', 'ice_back_edge_x_m', 'ice_bottom_z_m', 'crest_z_m', 'extrusion_width_m']:
            finite(getattr(self, k), k)
        if self.ice_back_edge_x_m <= self.upstream_edge_x_m or self.extrusion_width_m <= 0:
            raise ValueError('positive horizontal domain and width required')
        if not isinstance(self.bands, (tuple, list)) or not self.bands:
            raise ValueError('nonempty ordered bands required')
        object.__setattr__(self, 'bands', tuple(self.bands))
        for i, b in enumerate(self.bands):
            if not isinstance(b, Band):
                raise ValueError('bands must be Band records')
            if not self.upstream_edge_x_m < b.front_x_m < self.ice_back_edge_x_m:
                raise ValueError('front outside finite ice domain')
            if i and self.bands[i-1].z_high_m != b.z_low_m:
                raise ValueError('bands must share exact ordered interfaces; no silent gap repair')
        if self.bands[-1].z_high_m != self.crest_z_m or self.floor_z_m <= self.ice_bottom_z_m:
            raise ValueError('bands must span floor to fixed crest above ice bottom')
        finite(self.floor_area_m2,'derived floor area',positive=True)
        finite(self.vertical_face_area_m2,'derived face area',positive=True)
        finite(self.ice_volume_m3,'derived remaining volume',positive=True)

    @classmethod
    def uniform(cls, upstream_edge_x_m, ice_back_edge_x_m, ice_bottom_z_m,
                face_x_m, floor_z_m, crest_z_m, extrusion_width_m, n_bands=1):
        if isinstance(n_bands, bool) or not isinstance(n_bands, int) or n_bands <= 0:
            raise ValueError('positive integer band count required')
        for name, value in [('floor', floor_z_m), ('crest', crest_z_m), ('face', face_x_m)]:
            finite(value, name)
        if crest_z_m <= floor_z_m:
            raise ValueError('crest above floor required')
        bounds = [floor_z_m+(crest_z_m-floor_z_m)*i/n_bands for i in range(n_bands)]+[crest_z_m]
        bands = tuple(Band(a,b,face_x_m) for a,b in zip(bounds,bounds[1:]))
        return cls(upstream_edge_x_m, ice_back_edge_x_m, ice_bottom_z_m,
                   crest_z_m, extrusion_width_m, bands)

    @property
    def floor_z_m(self):
        return self.bands[0].z_low_m

    @property
    def floor_area_m2(self):
        return (self.bands[0].front_x_m-self.upstream_edge_x_m)*self.extrusion_width_m

    @property
    def vertical_face_area_m2(self):
        return (self.crest_z_m-self.floor_z_m)*self.extrusion_width_m

    @property
    def ice_volume_m3(self):
        base=(self.ice_back_edge_x_m-self.upstream_edge_x_m)*(self.floor_z_m-self.ice_bottom_z_m)
        face=math.fsum((self.ice_back_edge_x_m-b.front_x_m)*(b.z_high_m-b.z_low_m) for b in self.bands)
        return self.extrusion_width_m*(base+face)

    def polygon(self):
        points=[(self.upstream_edge_x_m,self.ice_bottom_z_m),
                (self.ice_back_edge_x_m,self.ice_bottom_z_m),
                (self.ice_back_edge_x_m,self.crest_z_m),
                (self.bands[-1].front_x_m,self.crest_z_m)]
        for i in range(len(self.bands)-1,-1,-1):
            b=self.bands[i];points.append((b.front_x_m,b.z_low_m))
            if i:
                points.append((self.bands[i-1].front_x_m,b.z_low_m))
        points.append((self.upstream_edge_x_m,self.floor_z_m))
        return points

    def ledges(self):
        rows=[]
        for lower, upper in zip(self.bands,self.bands[1:]):
            difference=lower.front_x_m-upper.front_x_m
            if difference:
                rows.append(dict(z_m=lower.z_high_m,x0_m=min(lower.front_x_m,upper.front_x_m),
                                 x1_m=max(lower.front_x_m,upper.front_x_m),
                                 area_m2=abs(difference)*self.extrusion_width_m,
                                 orientation='downward_ceiling' if difference>0 else 'upward_terrace'))
        return rows

    def undercut_diagnostic(self):
        maximum=0.; upper_min=None
        for b in reversed(self.bands):
            if upper_min is not None:
                maximum=max(maximum,b.front_x_m-upper_min)
            upper_min=b.front_x_m if upper_min is None else min(upper_min,b.front_x_m)
        return dict(maximum_horizontal_undercut_m=maximum,
                    downward_ceiling_area_m2=math.fsum(x['area_m2'] for x in self.ledges() if x['orientation']=='downward_ceiling'),
                    meaning='geometric diagnostic only; no fracture, stability or collapse prediction')


def advance(profile, band_latent_J, floor_latent_J, rho_ice=RHO_ICE, latent_heat=LATENT):
    """Apply disjoint face rectangles first, then floor on the new basal area.

    Positive x displacement means horizontal retreat; positive floor lowering is
    a negative z displacement. All inputs are NET latent-interface energies.
    The caller handles cold warming/storage and horizontal-ledge heat separately.
    """
    if not isinstance(profile, Profile):
        raise ValueError('validated Profile required')
    rho=finite(rho_ice,'ice density',positive=True);latent=finite(latent_heat,'latent heat',positive=True)
    if not isinstance(band_latent_J,(list,tuple)) or len(band_latent_J)!=len(profile.bands):
        raise ValueError('one face latent energy per existing height band required')
    energies=[finite(e,'band latent energy',nonnegative=True) for e in band_latent_J]
    floor_energy=finite(floor_latent_J,'floor latent energy',nonnegative=True)
    volumetric=finite(rho*latent,'volumetric latent heat',positive=True)
    new=[];dxs=[]
    for b,energy in zip(profile.bands,energies):
        area=(b.z_high_m-b.z_low_m)*profile.extrusion_width_m
        dx=energy/(volumetric*area)
        x=b.front_x_m+dx
        if energy>0 and x==b.front_x_m:
            raise ValueError('face change below coordinate floating-point resolution')
        new.append(replace(b,front_x_m=x));dxs.append(x-b.front_x_m)
    after_face=replace(profile,bands=tuple(new))
    floor_area=after_face.floor_area_m2
    dz=floor_energy/(volumetric*floor_area)
    floor=profile.floor_z_m-dz
    if floor_energy>0 and floor==profile.floor_z_m:
        raise ValueError('floor change below coordinate floating-point resolution')
    new[0]=replace(new[0],z_low_m=floor)
    after=replace(profile,bands=tuple(new))
    face_volume=math.fsum(energies)/volumetric;floor_volume=floor_energy/volumetric
    geometric_face=math.fsum(dx*(b.z_high_m-b.z_low_m)*profile.extrusion_width_m for dx,b in zip(dxs,profile.bands))
    geometric_floor=(profile.floor_z_m-after.floor_z_m)*floor_area
    stable_loss=geometric_face+geometric_floor
    subtraction=profile.ice_volume_m3-after.ice_volume_m3
    return after,dict(face_band_retreat_m=dxs,floor_lowering_m=profile.floor_z_m-after.floor_z_m,
        face_melt_m3=face_volume,floor_melt_m3=floor_volume,total_melt_m3=face_volume+floor_volume,
        stable_geometry_volume_loss_m3=stable_loss,
        remaining_volume_subtraction_m3=subtraction,
        volume_residual_m3=stable_loss-face_volume-floor_volume,
        remaining_subtraction_roundoff_m3=subtraction-stable_loss,
        phase_energy_J=math.fsum(energies)+floor_energy,
        floor_area_used_m2=floor_area,after_ledges=after.ledges(),
        scope='prescribed vertical-band and floor melt; newly exposed ledge heat not supplied here; no collapse')


def face_partition(profile, water_stage_m, debris_intervals=()):
    """Height-resolved exclusive vertical-face areas under a supplied water stage.

    Debris intervals are Eulerian height masks, not advected material parcels.
    Horizontal ledges are reported separately and are NOT silently assigned a
    vertical-face heat flux. Hydrostatic wetting/accessibility is an assumption.
    """
    finite(water_stage_m,'water stage')
    if not isinstance(debris_intervals,(list,tuple)):
        raise ValueError('ordered nonoverlapping debris-height intervals required')
    cover=[]
    for pair in debris_intervals:
        if not isinstance(pair,(list,tuple)) or len(pair)!=2:
            raise ValueError('two debris interval bounds required')
        lo,hi=[finite(v,'debris interval bound') for v in pair]
        if hi<=lo or cover and lo<cover[-1][1]:
            raise ValueError('debris intervals must be positive, ordered and nonoverlapping')
        cover.append((lo,hi))
    names=['bare_water_m2','bare_air_m2','debris_water_m2','debris_air_m2']
    totals={k:0. for k in names};rows=[]
    for i,b in enumerate(profile.bands):
        edges={b.z_low_m,b.z_high_m}
        for z in [water_stage_m]+[z for pair in cover for z in pair]:
            if b.z_low_m<z<b.z_high_m:edges.add(z)
        edges=sorted(edges);areas={k:0. for k in names}
        for lo,hi in zip(edges,edges[1:]):
            mid=lo+(hi-lo)/2
            covered=any(a<=mid<c for a,c in cover)
            wet=mid<water_stage_m
            name=('debris' if covered else 'bare')+('_water_m2' if wet else '_air_m2')
            areas[name]+=(hi-lo)*profile.extrusion_width_m
        for k in names:totals[k]+=areas[k]
        rows.append(dict(index=i,z_low_m=b.z_low_m,z_high_m=b.z_high_m,front_x_m=b.front_x_m,
                         **areas,total_vertical_area_m2=(b.z_high_m-b.z_low_m)*profile.extrusion_width_m))
    return dict(bands=rows,totals=totals,total_vertical_area_m2=profile.vertical_face_area_m2,
                area_residual_m2=math.fsum(totals.values())-profile.vertical_face_area_m2,
                horizontal_ledges=profile.ledges(),
                scope='conditional vertical-face height masks and hydrostatic stage; no material advection/air trapping/ledge heat')


def front_at(profile,z_m):
    finite(z_m,'query height')
    for b in profile.bands:
        if b.z_low_m<=z_m<b.z_high_m or z_m==profile.crest_z_m and b is profile.bands[-1]:
            return b.front_x_m
    raise ValueError('height outside exposed vertical profile')
