"""Minimum-barrier path through a plan-view elevation grid (four-neighbor).

The result is a geometric, uniform-level reachability threshold, NOT a discharge
prediction. Momentum, losses, time-dependent level and porous flow need the
hydrodynamic model. Do not mark porous material blocked by default.
"""
import heapq
import numpy as np


def minimax_crest(elevation,blocked=None):
    z=np.asarray(elevation,dtype=float)
    if z.ndim!=2 or z.shape[0]<2 or z.shape[1]<1 or not np.all(np.isfinite(z)):
        raise ValueError('finite 2-D grid with at least two streamwise rows required')
    walls=np.zeros(z.shape,dtype=bool) if blocked is None else np.asarray(blocked,dtype=bool)
    if walls.shape!=z.shape:raise ValueError('wall mask shape mismatch')
    nr,nc=z.shape;cost=np.full(z.shape,np.inf);parent={};heap=[]
    for j in range(nc):
        if not walls[0,j]:cost[0,j]=z[0,j];heapq.heappush(heap,(float(z[0,j]),0,j))
    goal=None
    while heap:
        value,i,j=heapq.heappop(heap)
        if value!=cost[i,j]:continue
        if i==nr-1:goal=(i,j);break
        for ni,nj in ((i-1,j),(i+1,j),(i,j-1),(i,j+1)):
            if not 0<=ni<nr or not 0<=nj<nc or walls[ni,nj]:continue
            candidate=max(value,float(z[ni,nj]))
            if candidate<cost[ni,nj]:
                cost[ni,nj]=candidate;parent[(ni,nj)]=(i,j);heapq.heappush(heap,(candidate,ni,nj))
    if goal is None:return dict(crest=None,path=[],status='no face-connected path')
    path=[goal]
    while path[-1] in parent:path.append(parent[path[-1]])
    path.reverse()
    straight=[]
    for j in range(nc):
        if not walls[:,j].any():straight.append(float(z[:,j].max()))
    return dict(crest=float(cost[goal]),path=[list(v) for v in path],
                straight_column_crest=min(straight) if straight else None,
                status='geometric path found',
                scope='lowest max elevation along a four-neighbor inlet-to-outlet path; not a flow solution')
