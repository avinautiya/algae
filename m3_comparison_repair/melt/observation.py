"""Identical-planform melt observation operators for the pinned M3 campaign.

Pure analysis/export functions. They never run or patch a solver. Volume means
vertical ice lowering integrated over planform area, not face-normal retreat.
"""
import math
import numpy as np


def positive(value,name):
    if isinstance(value,bool) or not math.isfinite(value) or value<=0:raise ValueError(name+' must be finite positive')
    return float(value)


def rectangle_overlap(x_edges,y_edges,rectangle):
    """Exact intersection area with piecewise-constant rectangular grid cells."""
    x=np.asarray(x_edges,float);y=np.asarray(y_edges,float);r=np.asarray(rectangle,float)
    if x.ndim!=1 or y.ndim!=1 or len(x)<2 or len(y)<2 or r.shape!=(4,):raise ValueError('edge/rectangle shape')
    if not all(np.isfinite(v).all() for v in (x,y,r)) or np.any(np.diff(x)<=0) or np.any(np.diff(y)<=0):raise ValueError('finite increasing edges required')
    x0,x1,y0,y1=r
    if x1<=x0 or y1<=y0:raise ValueError('rectangle must have positive extent')
    wx=np.maximum(0.,np.minimum(x[1:],x1)-np.maximum(x[:-1],x0))
    wy=np.maximum(0.,np.minimum(y[1:],y1)-np.maximum(y[:-1],y0))
    return wx[:,None]*wy[None,:]


def geometry_melt(before,after,x_edges,y_edges,rectangle,ice_mask=None,negative_tolerance_m=1e-12):
    """Observe actual before-minus-after bed change; refuse resolved refreezing.

    Ghost cells must be removed by the caller. Tolerance-level negative changes
    remain signed and are reported; they are not silently clipped into melt.
    """
    a=np.asarray(before,float);b=np.asarray(after,float);weights=rectangle_overlap(x_edges,y_edges,rectangle)
    if a.shape!=weights.shape or b.shape!=a.shape or not np.isfinite(a).all() or not np.isfinite(b).all():raise ValueError('finite field shapes must match grid')
    if not math.isfinite(negative_tolerance_m) or negative_tolerance_m<0:raise ValueError('invalid tolerance')
    mask=np.ones(a.shape,bool) if ice_mask is None else np.asarray(ice_mask)
    if mask.shape!=a.shape or mask.dtype!=np.bool_:raise ValueError('Boolean ice mask shape mismatch')
    dz=a-b
    if np.any((dz < -negative_tolerance_m)&mask&(weights>0)):raise ValueError('resolved bed rise is outside melt-only scope')
    return dict(volume_m3=float(np.sum(dz*weights*mask)),observed_ice_area_m2=float(np.sum(weights*mask)),
                most_negative_lowering_m=float(np.min(dz[mask&(weights>0)])) if np.any(mask&(weights>0)) else None,
                basis='actual vertical ice-bed lowering over declared planform rectangle')


def latent_energy(volume_m3,rho_i=917.,latent_J_kg=334000.):
    if not math.isfinite(volume_m3) or volume_m3<0:raise ValueError('finite nonnegative volume required')
    return volume_m3*positive(rho_i,'rho_i')*positive(latent_J_kg,'latent')


def volume_from_latent(energy_J,rho_i=917.,latent_J_kg=334000.):
    if not math.isfinite(energy_J) or energy_J<0:raise ValueError('finite nonnegative energy required')
    return energy_J/(positive(rho_i,'rho_i')*positive(latent_J_kg,'latent'))


def profile_melt(records,*,dx_m,toe_width_m,reach_start_m,roi_start_m,roi_end_m,rounding_decimals=5):
    """Existing CFD toe_floor_x exports suffice for a full-width axial ROI.

    Profile is the full transverse mean, not an arbitrary local cross-section.
    Initial/final errors from independently rounded profiles are bounded by
    10**(-decimals) times the observed area. This is a deterministic storage
    precision bound, not numerical model uncertainty or a confidence interval.
    """
    dx=positive(dx_m,'dx');width=positive(toe_width_m,'toe width')
    if not isinstance(records,list) or len(records)<2:raise ValueError('at least initial and final records required')
    first,last=records[0],records[-1]
    if not all(not isinstance(r['day'],bool) and math.isfinite(r['day']) for r in (first,last)) or last['day']<=first['day']:
        raise ValueError('finite strictly increasing profile endpoint times required')
    a=np.asarray(first['toe_floor_x'],float);b=np.asarray(last['toe_floor_x'],float)
    if a.ndim!=1 or not len(a) or b.shape!=a.shape or not np.isfinite(a).all() or not np.isfinite(b).all():raise ValueError('finite equal-length mean profiles required')
    if type(rounding_decimals) is not int or not 0<=rounding_decimals<=15:raise ValueError('invalid rounding precision')
    edges=reach_start_m+np.arange(len(a)+1)*dx
    if roi_start_m<edges[0]-1e-12 or roi_end_m>edges[-1]+1e-12 or roi_end_m<=roi_start_m:raise ValueError('ROI must lie inside exported reach')
    weights=rectangle_overlap(edges,[0.,width],[roi_start_m,roi_end_m,0.,width])[:,0]
    volume=float(np.dot(a-b,weights));precision=float(weights.sum()*10**(-rounding_decimals))
    if volume < -precision:raise ValueError('net bed rise exceeds rounding allowance')
    return dict(volume_m3=volume,rounding_volume_bound_m3=precision,area_m2=float(weights.sum()),
                first_day=float(first['day']),last_day=float(last['day']),
                roi=dict(x_start_m=roi_start_m,x_end_m=roi_end_m,toe_width_m=width),
                basis='initial/final transverse-mean ice elevation; no simulation rerun',
                warning='rounding bound excludes discretization, physical and closure uncertainty')


def reduced_axial_roi(total_toe_volume_m3,*,method,layout,reach_length_m=2.,
                      cover_fraction=.5,roi_start_from_reach_m=.5,roi_end_from_reach_m=1.5):
    """Map pinned reduced scalar to FULL-WIDTH axial ROI under its uniform cells.

    Connectivity strip represents ONLY a central cover_fraction*reach control.
    Corridor branches and fractional/fixed models cover the whole reach and are
    axially uniform. No unmodelled open-reach melt is invented for strip models.
    """
    if not math.isfinite(total_toe_volume_m3) or total_toe_volume_m3<0:raise ValueError('nonnegative finite reduced volume required')
    L=positive(reach_length_m,'reach length');a=roi_start_from_reach_m;b=roi_end_from_reach_m
    if not all(math.isfinite(v) for v in (a,b,cover_fraction)) or not 0<=a<b<=L or not 0<cover_fraction<1:raise ValueError('invalid ROI or cover')
    if layout=='porous_strip':raise ValueError('pinned reduced methods have no porous prediction')
    if method=='connectivity_reduced':
        if layout=='strip':start=(1-cover_fraction)*L/2;end=(1+cover_fraction)*L/2
        elif layout=='corridor':start,end=0.,L
        else:raise ValueError('unknown connectivity layout')
    elif method in ('fractional_cover','fixed_geometry','bare_reduced(evolving)','bare_fixed_geometry'):
        start,end=0.,L
    else:raise ValueError('unknown reduced method')
    if a<start-1e-12 or b>end+1e-12:raise ValueError('requested ROI includes unmodelled reduced reach')
    fraction=(b-a)/(end-start)
    return dict(volume_m3=float(total_toe_volume_m3*fraction),fraction_of_model_support=float(fraction),
                model_support_start_m=start,model_support_end_m=end,
                basis='exact axial restriction of declared uniform reduced controls; full toe width')


def export_sim_interval(sim,before_zb,*,roi_rectangle,dt_s):
    """Optional future pure hook immediately AFTER sim.lower_bed, with saved bed.

    Uses actual bed change, excludes ghosts and wall/nonice cells. Compares heat
    implied volume with actual lowering so floor clipping cannot hide a defect.
    Physical-time dt is supplied explicitly; hydraulic spin-up time is invalid.
    """
    dt=positive(dt_s,'physical morphology interval')
    edges_x=np.arange(sim.nx+1)*sim.dx;edges_y=np.arange(sim.ny+1)*sim.dx
    mask=sim.ice[1:-1,1:-1]&~sim.wall[1:-1,1:-1]
    observation=geometry_melt(np.asarray(before_zb)[1:-1,1:-1],sim.zb[1:-1,1:-1],
                              edges_x,edges_y,roi_rectangle,mask)
    weights=rectangle_overlap(edges_x,edges_y,roi_rectangle)
    heat_J=float(np.sum(sim.last_qacc[1:-1,1:-1]*weights*mask)*dt)
    observation.update(heat_removed_J=heat_J,heat_implied_volume_m3=volume_from_latent(heat_J),
                       geometric_minus_heat_volume_m3=observation['volume_m3']-volume_from_latent(heat_J))
    return observation
