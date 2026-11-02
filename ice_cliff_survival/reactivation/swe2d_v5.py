"""Version5: explicit meltwater source and auditable split-time accounting.

The inherited hydraulic kernel is the exact pinned v4 kernel. New
step_physical(dt) advances water/heat for dt then converts that SAME removed
heat to ice melt and liquid mass. It is first-order operator splitting, not
VOF, cold ice, or a full mechanical-energy/ALE solver. The legacy hydraulics
+ lower_bed API remains an explicitly quasi-steady accelerated approximation;
its hydraulic spin-up time and morphological time are NEVER conflated.
"""
import math
import inspect
import numpy as np
if __package__:
    from . import swe2d_v4_pinned as base
else:
    import swe2d_v4_pinned as base

G,RHO_W,CP,RHO_I,LF,HDRY=(base.G,base.RHO_W,base.CP,base.RHO_I,base.LF,base.HDRY)
SOLVER_VERSION='v5-phase-source-physical-split-1'

def positive(x,name):
    if isinstance(x,(bool,np.bool_)) or not np.isscalar(x) or not math.isfinite(x) or x<=0:
        raise ValueError(name+' must be finite and positive')
    return float(x)

class Sim(base.Sim):
    def __init__(self,*args,**kwargs):
        supplied=inspect.signature(base.Sim.__init__).bind(self,*args,**kwargs)
        supplied.apply_defaults();cfg=supplied.arguments
        for name in ('dx','Lpool','Lr','Ltail','Wl','heff','k_deb','cfl'):
            positive(cfg[name],name)
        for name in ('Q','Tw','manning'):
            value=cfg[name]
            if isinstance(value,(bool,np.bool_)) or not np.isscalar(value) or not np.isfinite(value) or value<0:
                raise ValueError(name+' must be finite nonnegative')
        if cfg['cfl']>.5:raise ValueError('this2D explicit kernel requires cfl<=0.5')
        if int(round(cfg['Wl']/cfg['dx']))<1:raise ValueError('lane width must contain a grid cell')
        self.source_log=[];self.hydraulic_window_log=[];self.reset_log=[]
        self.physical_time_s=0.;self.morphological_time_s=0.;self.hydraulic_time_s=0.
        self._heat_epoch=0;self._heat_ready=False;self._physical_sources_only=True
        self._phase_water=0.;self._phase_ice=0.;self._phase_latent=0.
        self._boundary_water=0.;self._positivity_water=0.;self._boundary_heat=0.;self._sink_heat=0.
        self._reset_water=0.;self._reset_heat=0.;self._reference_volume=None;self._reference_heat=None
        super().__init__(*args,**kwargs)

    @property
    def interior(self):return np.s_[1:-1,1:-1]

    def liquid_volume(self):return float(self.h[self.interior].sum()*self.dx**2)
    def sensible_heat(self):return float(RHO_W*CP*self.hT[self.interior].sum()*self.dx**2)

    def reset_state(self,eta):
        if isinstance(eta,(bool,np.bool_)) or not np.isscalar(eta) or not np.isfinite(eta):
            raise ValueError('eta must be finite')
        old_v=self.liquid_volume() if hasattr(self,'h') else None
        old_e=self.sensible_heat() if hasattr(self,'hT') else None
        super().reset_state(eta)
        v=self.liquid_volume();e=self.sensible_heat()
        if self._reference_volume is None:
            self._reference_volume=v;self._reference_heat=e
        elif not self.hydraulic_window_log and not self.source_log:
            # A pre-evolution reset defines the initial condition, not a
            # fictitious physical inflow. No completed evolution is erased.
            self._reference_volume=v;self._reference_heat=e
            self.reset_log.append(dict(kind='initial_condition_definition',volume_m3=v,sensible_J=e))
        elif old_v is not None:
            self._reset_water+=v-old_v;self._reset_heat+=e-old_e
            self.reset_log.append(dict(water_delta_m3=v-old_v,sensible_delta_J=e-old_e,
                meaning='explicit state replacement, not physical boundary flow'))
            self._physical_sources_only=False
        self._heat_ready=False

    def _window(self,duration,cool):
        duration=positive(duration,'hydraulic window')
        h0=self.liquid_volume();e0=self.sensible_heat()
        qacc=np.zeros_like(self.zb);fa=np.zeros(5)
        pk,phi,D=(1.,self.porous[1],self.porous[2]) if self.porous is not None else (0.,.4,.2)
        elapsed,steps,smax=base.advance(self.h,self.hu,self.hv,self.hT,self.zs(),self.pdep,
            self.heff_field(),self.ice,self.manning,self.Q/(self.ny*self.dx),self.Tw,
            self.dx,self.nx,self.ny,duration,self.cfl,pk,phi,D,1e-6,cool,qacc,fa,self.closed_out)
        v=self.liquid_volume();e=self.sensible_heat()
        active=np.zeros_like(self.ice);active[self.interior]=self.ice[self.interior]&~self.wall[self.interior]
        sink=float((qacc*active).sum()*self.dx**2)
        record=dict(elapsed_s=float(elapsed),inlet_m3=float(fa[0]),outlet_m3=float(fa[1]),
            positivity_added_m3=float(fa[4]),heat_in_J=float(fa[2]),heat_out_J=float(fa[3]),
            ice_heat_removed_J=sink,volume_before_m3=h0,volume_after_m3=v,
            sensible_before_J=e0,sensible_after_J=e,
            water_residual_m3=v-h0-fa[0]+fa[1]-fa[4],
            thermal_residual_J=e-e0-fa[2]+fa[3]+sink)
        self.hydraulic_window_log.append(record)
        self._boundary_water+=fa[0]-fa[1];self._positivity_water+=fa[4]
        self._boundary_heat+=fa[2]-fa[3];self._sink_heat+=sink
        self.hydraulic_time_s+=elapsed
        self.last_qacc=qacc/elapsed;self.last_fa=fa/elapsed
        self._dV=(v-h0)/elapsed;self._dE=(e-e0)/elapsed;self._win=elapsed
        self._heat_epoch+=1;self._heat_ready=True;self._source_window=record
        self.hyd_steps=steps;self.smax=smax
        return record,steps,smax

    def hydraulics(self,tmax=120.,tol=2e-4,window=4.,min_t=8.,cool=True):
        tmax=positive(tmax,'tmax');window=positive(window,'window');tol=positive(tol,'tol')
        if isinstance(min_t,(bool,np.bool_)) or not np.isscalar(min_t) or not np.isfinite(min_t) or min_t<0:raise ValueError('min_t must be finite nonnegative')
        # Merely executing a settling window creates pseudo-time heat/flow
        # history; subsequent physical steps cannot relabel that as physical.
        self._physical_sources_only=False
        total=0.;steps=0;smax=0.;last=None;converged=False
        while total<tmax:
            r,n,sm=self._window(min(window,tmax-total),cool)
            total+=r['elapsed_s'];steps+=n;smax=max(smax,sm)
            cur=np.array([self.Qsec(self.toe_cols()),self.Qsec(self.byp_cols()),
                self.last_fa[1],self.last_fa[3]/max(abs(self.last_fa[2]),1e-300)])
            qscale=max(abs(self.Q),1e-12)
            if last is not None and total>=min_t and abs(cur[0]-last[0])<tol*qscale and abs(cur[1]-last[1])<tol*qscale and abs(cur[3]-last[3])<1e-4 and abs(self._dV)<1e-3*qscale:
                converged=True;break
            last=cur
        self.hyd_time=total;self.hyd_steps=steps;self.smax=smax
        self._last_settlement=dict(converged=converged,hyd_time=total,steps=steps,
            dt_limit=self.cfl*self.dx/max(smax,1e-300),time_semantics='hydraulic pseudo-time; not morphological interval')
        return dict(self._last_settlement)

    def lower_bed(self,dtm):
        return self._lower_bed(dtm,mode='quasi_steady_morphological')

    def _lower_bed(self,dtm,mode):
        dtm=positive(dtm,'morphological interval')
        if not self._heat_ready:raise RuntimeError('fresh hydraulic heat window required; source already consumed or state reset')
        if self.incision not in ('full','uniform'):raise ValueError('unknown incision mode')
        active=np.zeros_like(self.ice);active[self.interior]=self.ice[self.interior]&~self.wall[self.interior]
        flux=np.where(active,self.last_qacc,0.)
        if not np.isfinite(flux).all() or (flux<0).any():raise ValueError('negative/nonfinite melt flux unsupported')
        dz=flux*dtm/(RHO_I*LF)
        if self.incision=='uniform':dz=active*(dz.sum()/max(active.sum(),1))
        if not np.isfinite(dz).all():raise ValueError('melt increment overflows')
        covered=active&((self.solid&(self.zd>0))|(self.pdep>0))
        covered_volume=float((dz*covered).sum()*self.dx**2)
        if mode=='physical_split' and covered_volume>0:
            raise ValueError('physical phase transfer beneath solid/porous debris requires a justified drainage/storage closure; nativev5 supports bare connected ice only')
        old_z=self.zb.copy();old_v=self.liquid_volume();old_e=self.sensible_heat()
        self.zb=self.zb-dz
        liquid=RHO_I/RHO_W*dz;self.h=self.h+liquid
        # New water is atTm and has zero horizontal momentum: hu,hv,hT unchanged.
        vi=float(dz.sum()*self.dx**2);vw=float(liquid.sum()*self.dx**2)
        volume_error=self.liquid_volume()-old_v-vw
        self.ledger_gap=float(abs((old_z-self.zb)[self.interior].sum()-dz[self.interior].sum())*self.dx**2)
        latent=RHO_I*LF*vi
        self._phase_water+=vw;self._phase_ice+=vi;self._phase_latent+=latent
        if mode!='physical_split':self._physical_sources_only=False
        self.morphological_time_s+=dtm
        source=dict(mode=mode,heat_epoch=self._heat_epoch,interval_s=dtm,
            ice_volume_m3=vi,water_added_m3=vw,
            liquid_volume_change_m3=self.liquid_volume()-old_v,
            volume_source_residual_m3=volume_error,
            phase_mass_residual_kg=RHO_W*vw-RHO_I*vi,
            latent_implied_J=latent,sensible_delta_J=self.sensible_heat()-old_e,
            matched_window_heat_removed_J=self._source_window['ice_heat_removed_J'],
            latent_minus_window_heat_J=latent-self._source_window['ice_heat_removed_J'],
            heat_time_scaling=dtm/self._source_window['elapsed_s'],
            covered_ice_melt_m3=covered_volume,
            under_debris_source_transport=(
                'assumed instantaneous transfer to the depth-averaged water; not physically verified'
                if covered_volume>0 else 'no covered phase source'),
            maximum_surface_drop_m=float(((1-RHO_I/RHO_W)*dz).max()),
            limitations='stationary temperate ice; no birth horizontal momentum/sensible heat; pressure-work/ALE mechanical energy not modelled')
        self.source_log.append(source);self.last_phase_source=source
        self._heat_ready=False
        rm=np.zeros_like(active);rm[self.i0:self.i1+1]=True
        self.last_reach=(float((dz*(self.toe&rm)).sum()*self.dx**2),float((dz*(self.byp&rm)).sum()*self.dx**2))
        for name,mask in self.regions.items():self.region_vol[name]=self.region_vol.get(name,0.)+float((dz*mask).sum()*self.dx**2)
        toe=float((dz*self.toe).sum()*self.dx**2);byp=float((dz*self.byp).sum()*self.dx**2)
        source.update(toe_ice_volume_m3=toe,bypass_ice_volume_m3=byp)
        source['route_source_diagnostics']={}
        for name,volume,cols in [('toe',toe,self.toe_cols()),('bypass',byp,self.byp_cols())]:
            q=self.Qsec(cols);water=RHO_I/RHO_W*volume
            source['route_source_diagnostics'][name]=dict(water_added_m3=water,
                sampled_section_discharge_m3_s=q,
                phase_to_absolute_section_throughflow=(water/(dtm*abs(q)) if q!=0 else None),
                zero_section_flow_with_generated_water=(q==0 and water>0),
                interpretation='diagnostic of a sampled section, not proof of local negligible feedback or uniform route discharge')
        return vi,toe,byp

    def step_physical(self,dt,cool=True):
        """One physical-time split step; uses identical heat and morph duration."""
        dt=positive(dt,'physical step')
        if self.incision!='full':
            raise ValueError('native physical mode requires local full incision; uniform redistribution is a nonlocal counterfactual')
        active=np.zeros_like(self.ice);active[self.interior]=self.ice[self.interior]&~self.wall[self.interior]
        if cool and np.any(active&((self.solid&(self.zd>0))|(self.pdep>0))):
            # Refuse BEFORE removing heat; do not leave a half-applied step.
            raise ValueError('native physical mode supports bare connected ice only; covered meltwater needs a drainage/storage closure')
        self._window(dt,cool)
        vols=self._lower_bed(dt,mode='physical_split')
        self.physical_time_s+=dt
        self.t+=dt
        return dict(dt_s=dt,physical_time_s=self.physical_time_s,ice_melt_m3=vols[0],
            source=dict(self.last_phase_source),ledgers=self.combined_ledgers())

    def combined_ledgers(self):
        mass_volume_res=self.liquid_volume()-self._reference_volume-self._boundary_water-self._positivity_water-self._phase_water-self._reset_water
        sensible_res=self.sensible_heat()-self._reference_heat-self._boundary_heat+self._sink_heat-self._reset_heat
        return dict(solver_version=SOLVER_VERSION,executed_state_water_volume_residual_m3=mass_volume_res,
            executed_state_sensible_residual_J=sensible_res,
            combined_ice_water_mass_residual_kg=RHO_W*(mass_volume_res+self._positivity_water),
            ice_volume_lost_m3=self._phase_ice,generated_water_m3=self._phase_water,
            hydraulic_boundary_net_water_m3=self._boundary_water,positivity_added_water_m3=self._positivity_water,
            explicit_reset_water_m3=self._reset_water,explicit_reset_heat_J=self._reset_heat,
            hydraulic_ice_heat_removed_J=self._sink_heat,morphological_latent_J=self._phase_latent,
            physical_split_thermal_latent_residual_J=(sensible_res+self._phase_latent-self._sink_heat if self._physical_sources_only else None),
            physical_split_only=self._physical_sources_only,
            scope='algorithmic control-volume ledger across executed hydraulic windows; pseudo-time throughflow must not be reported as environmental elapsed time')

    def diagnostics(self):
        if not hasattr(self,'last_fa'):raise RuntimeError('no hydraulic window available')
        qt=self.Qsec(self.toe_cols());qb=self.Qsec(self.byp_cols());fa=self.last_fa
        scale=max(abs(float(fa[0])),abs(float(fa[1])),abs(self.Q),1e-300)
        active=np.zeros_like(self.ice);active[self.interior]=self.ice[self.interior]&~self.wall[self.interior]
        melt=float((self.last_qacc*active).sum()*self.dx**2)
        return dict(Q_toe=qt,Q_byp=qb,f_toe=(qt/self.Q if self.Q!=0 else None),Q_in=fa[0],Q_out=fa[1],
            water_resid=(fa[0]-fa[1]+fa[4]-self._dV)/scale,water_storage_rate=self._dV,
            heat_in_W=fa[2],heat_out_W=fa[3],melt_heat_W=melt,
            energy_resid_rel=(fa[2]-fa[3]-melt-self._dE)/max(abs(fa[2]),abs(fa[3]),abs(melt),1.),storage_rate_W=self._dE,
            melt_heat_toe_W=float((self.last_qacc*self.toe).sum()*self.dx**2),
            melt_heat_byp_W=float((self.last_qacc*self.byp).sum()*self.dx**2),
            wet_ice_area_toe=float(((self.h>1e-3)&self.ice&self.toe).sum()*self.dx**2),
            wet_ice_area_byp=float(((self.h>1e-3)&self.ice&self.byp).sum()*self.dx**2),
            hydraulic_window_state_is_current=self._heat_ready,
            phase_source=self.last_phase_source if hasattr(self,'last_phase_source') else None,
            combined_ledgers=self.combined_ledgers())

    def region_heat(self):
        active=np.zeros_like(self.ice);active[self.interior]=self.ice[self.interior]&~self.wall[self.interior]
        return {name:float((self.last_qacc*mask*active).sum()*self.dx**2) for name,mask in self.regions.items()}

    def set_deposit(self,*args,**kwargs):
        result=super().set_deposit(*args,**kwargs)
        self._heat_ready=False
        return result

    def remove_debris(self):
        result=super().remove_debris()
        self._heat_ready=False
        return result
