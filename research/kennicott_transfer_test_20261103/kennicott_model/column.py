"""Passive heterogeneous-debris/cold-ice fixed-interface latent-sink column.

Backward-Euler FV; interfaceTb<=0,m>=0,Tb*m=0. KnownC reference is NOT fitted
by this solver. No target values, no fictitious below-ice thermistor forcing.
"""
import math

def factor(diag,lower,upper):
 d=list(diag);u=list(upper);l=list(lower);multipliers=[]
 for j in range(1,len(d)):
  if d[j-1]<=0:raise ValueError('nonpositive passive matrix pivot')
  m=l[j-1]/d[j-1];multipliers.append(m);d[j]-=m*u[j-1]
 if d[-1]<=0:raise ValueError('nonpositive pivot')
 def solve(rhs):
  b=list(rhs)
  for j,m in enumerate(multipliers,1):b[j]-=m*b[j-1]
  b[-1]/=d[-1]
  for j in range(len(d)-2,-1,-1):b[j]=(b[j]-u[j]*b[j+1])/d[j]
  return b
 return solve

def system(storage,faces,leftG,rightG,h):
 diag=[s/h+(faces[j-1]if j else leftG)+(faces[j]if j<len(storage)-1 else rightG)for j,s in enumerate(storage)]
 return factor(diag,[-v for v in faces],[-v for v in faces])

def simulate(*,surface_temperature_C,interval_seconds,debris_depth_m,debris_alpha_cells_m2_s,debris_capacity_J_m3K,
             initial_debris_temperature_C,initial_ice_temperature_C,ice_depth_m=2.,ice_k_W_mK=2.1,ice_C_J_m3K=917.*2100.,latent_J_kg=333700.,record_intervals=True):
 n=len(debris_alpha_cells_m2_s);m=len(initial_ice_temperature_C)
 if len(initial_debris_temperature_C)!=n or min(n,m)<1 or not surface_temperature_C:raise ValueError('nonempty matching explicit states')
 vals=list(surface_temperature_C)+list(initial_debris_temperature_C)+list(initial_ice_temperature_C)+list(debris_alpha_cells_m2_s)+[interval_seconds,debris_depth_m,debris_capacity_J_m3K,ice_depth_m,ice_k_W_mK,ice_C_J_m3K,latent_J_kg]
 if any(not math.isfinite(x)for x in vals)or min(interval_seconds,debris_depth_m,debris_capacity_J_m3K,ice_depth_m,ice_k_W_mK,ice_C_J_m3K,latent_J_kg)<=0 or min(debris_alpha_cells_m2_s)<=0:raise ValueError('finite positive material/time/geometry')
 if any(x>0 for x in initial_ice_temperature_C):raise ValueError('initialiceabove0unsupported')
 dz=debris_depth_m/n;dzi=ice_depth_m/m;kd=[a*debris_capacity_J_m3K for a in debris_alpha_cells_m2_s];Sd=[debris_capacity_J_m3K*dz]*n;Si=[ice_C_J_m3K*dzi]*m
 fd=[1/(dz/(2*a)+dz/(2*b))for a,b in zip(kd[:-1],kd[1:])];fi=[ice_k_W_mK/dzi]*(m-1)
 Gtop=2*kd[0]/dz;Gd=2*kd[-1]/dz;Gi=2*ice_k_W_mK/dzi;Gint=1/(1/Gd+1/Gi)
 h=interval_seconds;warmD=system(Sd,fd,Gtop,Gd,h);warmI=system(Si,fi,Gi,Gi,h);cold=system(Sd+Si,fd+[Gint]+fi,Gtop,Gi,h)
 Td=list(initial_debris_temperature_C);Ti=list(initial_ice_temperature_C);E0=sum(s*t for s,t in zip(Sd+Si,Td+Ti));Ein=Edeep=latent=0.;maxdef=0.;coldsteps=0;trace=[]
 for j,Ts in enumerate(surface_temperature_C):
  rd=[s*t/h for s,t in zip(Sd,Td)];rd[0]+=Gtop*Ts;ri=[s*t/h for s,t in zip(Si,Ti)];D=warmD(rd);I=warmI(ri);net=Gd*D[-1]+Gi*I[0]
  if net>=0:Tb=0.;melt=net*h/latent_J_kg
  else:
   allT=cold(rd+ri);D,I=allT[:n],allT[n:];Tb=(Gd*D[-1]+Gi*I[0])/(Gd+Gi);melt=0.;coldsteps+=1
   if Tb>1e-12:raise ValueError('coldbranchpositiveinterface')
  if any(t>1e-12 for t in I):raise ValueError('iceexceedsmeltingpoint')
  supply=h*Gtop*(Ts-D[0]);deep=h*Gi*I[-1];change=sum(s*(v-u)for s,u,v in zip(Sd+Si,Td+Ti,D+I));L=latent_J_kg*melt;defect=change+L-supply+deep
  Ein+=supply;Edeep+=deep;latent+=L;maxdef=max(maxdef,abs(defect));Td,Ti=D,I
  if record_intervals:trace.append({'start_seconds':j*h,'end_seconds':(j+1)*h,'melt_mmwe':melt,'interface_T_C':Tb,'surface_heat_J_m2':supply,'deep_heat_out_J_m2':deep,'sensible_storage_change_J_m2':change,'latent_heat_J_m2':L,'ledger_residual_J_m2':defect})
 Ef=sum(s*t for s,t in zip(Sd+Si,Td+Ti));res=Ef-E0+latent-Ein+Edeep;scale=max(abs(E0),abs(Ef),abs(latent),abs(Ein),abs(Edeep),1.)
 return {'physical_melt_mmwe':latent/latent_J_kg,'initial_storage_J_m2':E0,'final_storage_J_m2':Ef,'surface_heat_J_m2':Ein,'deep_heat_out_J_m2':Edeep,'latent_heat_J_m2':latent,'energy_residual_J_m2':res,'relative_energy_residual':res/scale,'maximum_step_ledger_residual_J_m2':maxdef,'cold_interface_steps':coldsteps,'trace':trace,'final_debris_temperature_C':Td,'final_ice_temperature_C':Ti,'fixed_geometry':True,'moving_interface_or_ice_advection':False,'wet_debris_phase_or_advection':False,'physical_capacity_fitted':False,'field_validation_established':False}
