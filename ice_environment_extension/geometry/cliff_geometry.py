"""Conservative geometry operator, not a cliff evolution or collapse closure.

Coordinates are (horizontal metres, elevation metres). Ice is a convex polygon.
A supporting face retreats inward by a parallel cut; supplied *net latent* energy
determines the removed area. Width is perpendicular to the cross section.
"""
from dataclasses import dataclass
import math


def signed_area(vertices):
    p = tuple(tuple(map(float, v)) for v in vertices)
    if len(p) < 3 or any(len(v) != 2 or not all(map(math.isfinite, v)) for v in p):
        raise ValueError("Need at least three finite 2-D vertices")
    # Shift origin to reduce cancellation at projected-map coordinates.
    x0, z0 = p[0]
    q = [(x-x0, z-z0) for x, z in p]
    return .5 * math.fsum(a[0]*b[1]-b[0]*a[1] for a,b in zip(q,q[1:]+q[:1]))


def convex_polygon(vertices):
    p = tuple(tuple(map(float, v)) for v in vertices)
    area = signed_area(p)
    if area == 0:
        raise ValueError("Zero area")
    if area < 0:
        p = p[::-1]
    scale = max(math.dist(a,b) for a,b in zip(p,p[1:]+p[:1]))
    tol = 1e-12 * scale**2
    for i in range(len(p)):
        a,b,c = p[i-1],p[i],p[(i+1)%len(p)]
        if math.dist(a,b) == 0:
            raise ValueError("Repeated adjacent vertex")
        if (b[0]-a[0])*(c[1]-b[1])-(b[1]-a[1])*(c[0]-b[0]) < -tol:
            raise ValueError("Only simple convex polygons are supported")
    # Convexity requires every vertex to lie on every edge's interior side.
    for a,b in zip(p,p[1:]+p[:1]):
        if any((b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]) < -tol for c in p):
            raise ValueError("Non-convex or self-intersecting polygon")
    return p


def _normal(n):
    n = tuple(map(float,n))
    if len(n)!=2 or not all(map(math.isfinite,n)) or math.hypot(*n)==0:
        raise ValueError("Finite nonzero normal required")
    norm=math.hypot(*n)
    return n[0]/norm,n[1]/norm


def _clip(p,n,c):
    out=[]
    for a,b in zip(p,p[1:]+p[:1]):
        da=sum(x*y for x,y in zip(a,n))-c
        db=sum(x*y for x,y in zip(b,n))-c
        if da>=0:
            out.append(a)
        if (da<0)!=(db<0):
            t=da/(da-db)
            out.append(tuple(x+t*(y-x) for x,y in zip(a,b)))
    return tuple(out)


@dataclass(frozen=True)
class RetreatResult:
    vertices: tuple
    normal_retreat_m: float
    removed_ice_m3: float
    latent_energy_J: float
    energy_residual_J: float


def retreat_by_energy(vertices, inward_normal, energy_J, width_m=1., rho_ice=917., latent_heat=334000.):
    """Single planar face update; energy is already net of all other budgets.

    Refuses ice exhaustion. This cannot infer a face direction from floor melt.
    A non-edge-aligned normal creates a corner cut, so callers must explicitly
    supply a physical planar boundary normal (see integration contract).
    """
    p=convex_polygon(vertices); n=_normal(inward_normal)
    if not all(math.isfinite(v) and v>0 for v in (width_m,rho_ice,latent_heat)):
        raise ValueError("Positive finite width, density and latent heat required")
    if not math.isfinite(energy_J) or energy_J<0:
        raise ValueError("Finite nonnegative net latent energy required; no refreezing")
    area=signed_area(p); target=energy_J/(width_m*rho_ice*latent_heat)
    if target>=area:
        raise ValueError("Ice exhaustion or geometry-domain change requires external handling")
    if target==0:
        return RetreatResult(p,0.,0.,0.,0.)
    proj=[sum(x*y for x,y in zip(v,n)) for v in p]
    base=min(proj); lo=base; hi=max(proj)
    for _ in range(90):
        mid=(lo+hi)/2
        q=_clip(p,n,mid)
        removed=area-(signed_area(q) if len(q)>=3 else 0.)
        if removed<target: lo=mid
        else: hi=mid
    cut=(lo+hi)/2; q=_clip(p,n,cut)
    volume=(area-signed_area(q))*width_m
    latent=volume*rho_ice*latent_heat
    return RetreatResult(q,cut-base,volume,latent,latent-energy_J)


def segment_exposure(start, end, water_level_m, debris_top_m):
    """Exact line lengths by elevation for a prescribed toe mantle.

    Debris covers everything at/below debris_top; remaining face partitions into
    submerged and exposed. Thus debris may itself be wet but is counted once.
    This height mask is an input, not a prediction of grain transport/coverage.
    """
    a=tuple(map(float,start)); b=tuple(map(float,end))
    if len(a)!=2 or len(b)!=2 or not all(math.isfinite(v) for v in (*a,*b,water_level_m,debris_top_m)):
        raise ValueError("Finite 2-D endpoints and levels required")
    length=math.dist(a,b)
    def below(level):
        low,high=sorted((a[1],b[1]))
        if high==low: return length if low<=level else 0.
        return length*max(0.,min(1.,(level-low)/(high-low)))
    debris=below(debris_top_m)
    submerged=max(0.,below(water_level_m)-debris)
    return {"debris_covered_m":debris,"submerged_bare_m":submerged,
            "exposed_bare_m":length-debris-submerged,"total_m":length}
