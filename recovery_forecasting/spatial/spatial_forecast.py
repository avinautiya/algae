"""Conditional spatial latent-energy exclusion CLI, not a flow/survival model."""
import argparse
import hashlib
import heapq
import json
import math
from pathlib import Path
import numpy as np

ASSUMPTIONS = (
    'ice_surface_after_clearing', 'global_local_free_surface_bounds',
    'melting_only_opening', 'all_usable_heat_bounded',
    'terrain_is_cellwise_bound', 'target_width_is_declared',
    'walls_are_immutable_and_no_unmodeled_routes',
)


def number(value, name, *, positive=False, nonnegative=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(name + ' must be a finite number, not a bool')
    if positive and value <= 0 or nonnegative and value < 0:
        raise ValueError(name + ' is outside its allowed range')
    return float(value)


def matrix(value, name):
    if not isinstance(value, list) or not value or not all(isinstance(row,list) and row for row in value):
        raise ValueError(name + ' must be a nonempty rectangular numeric matrix')
    if len({len(row) for row in value}) != 1:
        raise ValueError(name + ' must be rectangular')
    return np.array([[number(v,name) for v in row] for row in value],dtype=float)


def interval(value, name, *, positive=False, nonnegative=False):
    if not isinstance(value,dict) or set(value)!={'lower','upper'}:
        raise ValueError(name + ' requires exactly lower and upper')
    lo=number(value['lower'],name+'.lower',positive=positive,nonnegative=nonnegative)
    hi=number(value['upper'],name+'.upper',positive=positive,nonnegative=nonnegative)
    if lo>hi:raise ValueError(name + ' lower exceeds upper')
    return lo,hi


def validate(raw):
    if not isinstance(raw,dict):raise ValueError('input must be a JSON object')
    allowed={'schema_version','elevation_m','elevation_bounds_m','dx_m','dy_m',
             'minimum_cross_section_width_m','level_bounds_m','usable_heat_bounds_J',
             'ice_density_bounds_kg_m3','latent_heat_bounds_J_kg','blocked',
             'resolved_cell_route','assumptions','provenance'}
    if set(raw)-allowed:raise ValueError('unknown input fields: '+', '.join(sorted(set(raw)-allowed)))
    if type(raw.get('schema_version')) is not int or raw['schema_version']!=1:
        raise ValueError('schema_version must be integer 1')
    z=matrix(raw.get('elevation_m'),'elevation_m')
    eb=raw.get('elevation_bounds_m',{'lower':raw['elevation_m'],'upper':raw['elevation_m']})
    if not isinstance(eb,dict) or set(eb)!={'lower','upper'}:raise ValueError('elevation_bounds_m requires lower and upper')
    zl,zu=(matrix(eb[k],'elevation_bounds_m.'+k) for k in ('lower','upper'))
    if zl.shape!=z.shape or zu.shape!=z.shape or not np.all((zl<=z)&(z<=zu)):
        raise ValueError('elevation bounds must enclose nominal matrix with identical shapes')
    dx=number(raw.get('dx_m'),'dx_m',positive=True);dy=number(raw.get('dy_m'),'dy_m',positive=True)
    width=number(raw.get('minimum_cross_section_width_m',0.),'minimum_cross_section_width_m',nonnegative=True)
    level=interval(raw.get('level_bounds_m'),'level_bounds_m')
    energy=interval(raw.get('usable_heat_bounds_J'),'usable_heat_bounds_J',nonnegative=True)
    rho=interval(raw.get('ice_density_bounds_kg_m3',{'lower':917.,'upper':917.}),'ice_density_bounds_kg_m3',positive=True)
    latent=interval(raw.get('latent_heat_bounds_J_kg',{'lower':334000.,'upper':334000.}),'latent_heat_bounds_J_kg',positive=True)
    walls=raw.get('blocked',[[False]*z.shape[1] for _ in range(z.shape[0])])
    if not isinstance(walls,list) or len(walls)!=z.shape[0] or any(not isinstance(r,list) or len(r)!=z.shape[1] or any(type(v) is not bool for v in r) for r in walls):
        raise ValueError('blocked must be a Boolean matrix matching elevation')
    discrete=raw.get('resolved_cell_route',False)
    if type(discrete) is not bool:raise ValueError('resolved_cell_route must be Boolean')
    claims=raw.get('assumptions')
    if not isinstance(claims,dict) or set(claims)!=set(ASSUMPTIONS) or any(claims[k] is not True for k in ASSUMPTIONS):
        raise ValueError('assumptions must explicitly assert every documented prerequisite as true')
    if 'provenance' in raw and not isinstance(raw['provenance'],dict):raise ValueError('provenance must be an object')
    return dict(z=z,zlower=zl,zupper=zu,dx=dx,dy=dy,width=width,level=level,
                energy=energy,rho=rho,latent=latent,walls=np.array(walls,dtype=bool),discrete=discrete)


def volume_bound(z, level, dx, dy, width, walls, discrete):
    """Fractional row-width relaxation; optional WHOLE-CELL graph constraint."""
    with np.errstate(over='raise',invalid='raise'):
        try:depth=np.maximum(z-level,0.)
        except FloatingPointError:raise ValueError('numeric range exceeded in elevation difference')
    row_bound=0.
    for row,wall in zip(depth,walls):
        ordered=np.sort(row[~wall]);remaining=width
        if width>len(ordered)*dy:return dict(volume=math.inf,width=math.inf,path=math.inf,path_cells=[])
        for d in ordered:
            take=min(remaining,dy);row_bound+=float(d)*take*dx
            if not math.isfinite(row_bound):raise ValueError('numeric range exceeded in width cost')
            remaining=max(0.,remaining-take)
            if remaining==0:break
    with np.errstate(over='raise',invalid='raise'):
        try:cost=depth*dx*dy
        except FloatingPointError:raise ValueError('numeric range exceeded in cell cost')
    nr,nc=z.shape;dist=np.full(z.shape,np.inf);parents={};heap=[]
    for j in range(nc):
        if not walls[0,j]:dist[0,j]=cost[0,j];heapq.heappush(heap,(float(cost[0,j]),0,j))
    goal=None
    while heap:
        value,i,j=heapq.heappop(heap)
        if value!=dist[i,j]:continue
        if i==nr-1:goal=(i,j);break
        for ii,jj in ((i-1,j),(i+1,j),(i,j-1),(i,j+1)):
            if 0<=ii<nr and 0<=jj<nc and not walls[ii,jj]:
                candidate=value+float(cost[ii,jj])
                if not math.isfinite(candidate):raise ValueError('numeric range exceeded in path cost')
                if candidate<dist[ii,jj]:
                    dist[ii,jj]=candidate;parents[ii,jj]=(i,j);heapq.heappush(heap,(candidate,ii,jj))
    path=[]
    if goal is not None:
        path=[goal]
        while path[-1] in parents:path.append(parents[path[-1]])
        path.reverse()
    path_bound=float(dist[goal]) if goal is not None else math.inf
    return dict(volume=max(row_bound,path_bound) if discrete else row_bound,
                width=row_bound,path=path_bound,path_cells=[list(p) for p in path])


def finite_record(value):
    if isinstance(value,float) and math.isinf(value):return {'value':None,'infinite':True}
    return {'value':float(value),'infinite':False}


def forecast(raw):
    p=validate(raw)
    optimistic=volume_bound(p['zlower'],p['level'][1],p['dx'],p['dy'],p['width'],p['walls'],p['discrete'])
    pessimistic=volume_bound(p['zupper'],p['level'][0],p['dx'],p['dy'],p['width'],p['walls'],p['discrete'])
    elo=optimistic['volume']*p['rho'][0]*p['latent'][0]
    ehi=pessimistic['volume']*p['rho'][1]*p['latent'][1]
    if any(math.isfinite(v) and not math.isfinite(e) for v,e in ((optimistic['volume'],elo),(pessimistic['volume'],ehi))):
        raise ValueError('numeric range exceeded in latent energy')
    robust=p['energy'][1]<elo
    any_exclusion=p['energy'][0]<ehi
    status='excluded_for_all_declared_box_values' if robust else ('exclusion_depends_on_declared_values' if any_exclusion else 'not_excluded_for_any_declared_box_value')
    payload=json.dumps(raw,sort_keys=True,separators=(',',':'),allow_nan=False).encode()
    return dict(schema_version=1,model='conditional_spatial_latent_energy_exclusion_v1',
        input_sha256=hashlib.sha256(payload).hexdigest(),
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        classification=status,reconnection_excluded_under_all_declared_bounds=robust,
        assumption_status='user_asserted_not_independently_verified',
        feasible_class='whole_cell_graph_and_declared_width' if p['discrete'] else 'fractional_row_width_relaxation',
        required_latent_energy_enclosure_J={'lower':finite_record(elo),'upper':finite_record(ehi)},
        usable_heat_bounds_J={'lower':p['energy'][0],'upper':p['energy'][1]},
        optimistic_ice_volume_bound_m3=finite_record(optimistic['volume']),
        optimistic_width_only_volume_m3=finite_record(optimistic['width']),
        diagnostic_whole_cell_path_volume_m3=finite_record(optimistic['path']),
        diagnostic_path_only_optimizer=optimistic['path_cells'],
        minimum_cross_section_width_m=p['width'],
        warnings=['No positive discharge, sustained contact or cliff-survival prediction follows from non-exclusion.',
                  'The parameter box is a deterministic enclosure, not a confidence interval.',
                  'The optional graph path is a whole-cell diagnostic and is not a continuum trench bound.',
                  'Upstream stage alone is not a global local-free-surface upper bound.',
                  'Terrain bounds must bound subcell terrain; cell-centre samples alone are insufficient.'],
        provenance=raw.get('provenance',{}))


def strict_json(path):
    def invalid(token):raise ValueError('nonstandard nonfinite JSON constant: '+token)
    def unique(pairs):
        out={}
        for k,v in pairs:
            if k in out:raise ValueError('duplicate JSON field: '+k)
            out[k]=v
        return out
    return json.loads(Path(path).read_text(),parse_constant=invalid,object_pairs_hook=unique)


def plot(raw, record, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator, FormatStrFormatter
    p=validate(raw);levels=np.linspace(p['level'][0],p['level'][1],101)
    costs=[volume_bound(p['zlower'],h,p['dx'],p['dy'],p['width'],p['walls'],p['discrete'])['volume']*p['rho'][0]*p['latent'][0] for h in levels]
    fig,axs=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    im=axs[0].imshow(np.ma.array(p['zlower'],mask=p['walls']),origin='lower',aspect='auto',extent=(0,p['z'].shape[1]*p['dy'],0,p['z'].shape[0]*p['dx']))
    fig.colorbar(im,ax=axs[0],label='Optimistic ice elevation bound (m)')
    axs[0].set(xlabel='Transverse distance (m)',ylabel='Streamwise distance (m)',title='Post-clearance ice surface')
    axs[1].plot(levels,np.array(costs)/1e6,label='Optimistic required latent energy')
    axs[1].axhline(p['energy'][1]/1e6,ls=':',color='black',label='Usable heat upper bound')
    axs[1].set(xlabel='Global local-water-level upper bound (m)',ylabel='Energy (MJ)',title='Necessary opening condition')
    axs[1].xaxis.set_major_locator(MaxNLocator(5))
    axs[1].xaxis.set_major_formatter(FormatStrFormatter('%.3g'))
    axs[1].legend(fontsize=8)
    label='Synthetic method check' if 'SYNTHETIC' in str(raw.get('provenance',{}).get('evidence','')).upper() else 'Declared-input diagnostic'
    fig.suptitle(label+' — energy exclusion only, not a flow forecast',fontsize=11)
    fig.savefig(path,dpi=180);plt.close(fig)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input');parser.add_argument('--output');parser.add_argument('--plot')
    a=parser.parse_args(argv)
    try:
        raw=strict_json(a.input);result=forecast(raw)
        content=json.dumps(result,indent=2,allow_nan=False)+'\n'
        if a.output:Path(a.output).write_text(content)
        else:print(content,end='')
        if a.plot:plot(raw,result,a.plot)
        return 0
    except (ValueError,TypeError,KeyError,OSError,OverflowError) as exc:
        parser.exit(2,json.dumps({'error':str(exc),'status':'invalid_input_or_output'})+'\n')


if __name__=='__main__':raise SystemExit(main())
