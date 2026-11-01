"""Optimistic latent-energy lower bounds for opening a resolved ice route.

These are necessary conditions, not hydraulic forecasts. Elevation must be the
removable ice surface after any assumed perfect debris clearing. No collapse,
uplift, refreezing, or time-dependent vertical datum is represented.
"""
import heapq
import math
import numpy as np


def opening_bound(elevation, *, level_upper_m, dx_m, dy_m,
                  minimum_cross_section_width_m=0., resolved_cell_route=False, blocked=None,
                  rho_i=917., latent_J_kg=334000.):
    """Inlet = first row; target = last row; four-neighbor grid connectivity.

    The path bound is exact for an arbitrarily allocated single-cell-wide route
    on this grid. The width bound relaxes connectivity and permits disjoint
    cross-sectional pieces, so it is optimistic for every corridor occupying at
    least the declared transverse width at EACH row. Their maximum remains a
    lower bound only if whole-cell routes are required. The default uses ONLY
    the relaxed width bound, permitting sub-cell corridors. Set
    resolved_cell_route=True to explicitly require full cells along a route.
    It is not an estimate of the actual heat delivered to ice.
    """
    z=np.asarray(elevation,dtype=float)
    if z.ndim!=2 or min(z.shape)<1 or not np.isfinite(z).all():
        raise ValueError('finite nonempty 2-D elevation required')
    scalars=(level_upper_m,dx_m,dy_m,minimum_cross_section_width_m,rho_i,latent_J_kg)
    if not all(not isinstance(v,(bool,np.bool_)) and math.isfinite(v) for v in scalars) or min(dx_m,dy_m,rho_i,latent_J_kg)<=0 or minimum_cross_section_width_m<0:
        raise ValueError('invalid geometry or material inputs')
    walls=np.zeros(z.shape,dtype=bool) if blocked is None else np.asarray(blocked,dtype=bool)
    if walls.shape!=z.shape:raise ValueError('wall shape mismatch')
    depth=np.maximum(z-level_upper_m,0.)
    nr,nc=z.shape
    cost=depth*dx_m*dy_m
    dist=np.full(z.shape,np.inf);parent={};queue=[]
    for j in range(nc):
        if not walls[0,j]:
            dist[0,j]=cost[0,j];heapq.heappush(queue,(float(cost[0,j]),0,j))
    goal=None
    while queue:
        value,i,j=heapq.heappop(queue)
        if value!=dist[i,j]:continue
        if i==nr-1:goal=(i,j);break
        for ii,jj in ((i-1,j),(i+1,j),(i,j-1),(i,j+1)):
            if not (0<=ii<nr and 0<=jj<nc) or walls[ii,jj]:continue
            candidate=value+cost[ii,jj]
            if candidate<dist[ii,jj]:
                dist[ii,jj]=candidate;parent[ii,jj]=(i,j)
                heapq.heappush(queue,(candidate,ii,jj))
    path=[]
    if goal is not None:
        path=[goal]
        while path[-1] in parent:path.append(parent[path[-1]])
        path.reverse()
    path_volume=float(dist[goal]) if goal is not None else math.inf
    width_volume=0.
    for row,mask in zip(depth,walls):
        remaining=minimum_cross_section_width_m
        available=np.sort(row[~mask])
        if remaining>len(available)*dy_m:
            width_volume=math.inf;break
        for d in available:
            take=min(remaining,dy_m)
            width_volume+=float(d)*take*dx_m
            remaining=max(0.,remaining-take)
            if remaining==0:break
    if not isinstance(resolved_cell_route,bool):raise ValueError('resolved_cell_route must be bool')
    volume=max(path_volume,width_volume) if resolved_cell_route else width_volume
    return dict(minimum_ice_volume_lower_bound_m3=volume,
                latent_energy_lower_bound_J=volume*rho_i*latent_J_kg,
                path_only_volume_m3=path_volume,
                relaxed_width_volume_m3=width_volume,
                path_only_optimizer=[list(p) for p in path],
                level_upper_m=float(level_upper_m),
                minimum_cross_section_width_m=float(minimum_cross_section_width_m),
                resolved_cell_route=resolved_cell_route,
                scope='necessary static geometric condition; no discharge or sustained-contact guarantee')


def exclusion(bound, available_energy_upper_J):
    if isinstance(available_energy_upper_J,bool) or not math.isfinite(available_energy_upper_J) or available_energy_upper_J<0:
        raise ValueError('finite nonnegative upper heat budget required')
    required=bound['latent_energy_lower_bound_J']
    if isinstance(required,bool) or math.isnan(required) or required<0:
        raise ValueError('nonnegative required energy (possibly infinite) required')
    return dict(reconnection_excluded=available_energy_upper_J<required,
                available_energy_upper_J=float(available_energy_upper_J),
                required_energy_lower_J=required,
                warning='failure to exclude does not predict reconnection; heat allocation is deliberately optimistic')
