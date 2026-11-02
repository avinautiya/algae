from run import *
r=json.loads(Path(__file__).with_name('results.json').read_text());rows=[]
for col in r['columns']:
 i=col['profile'];d=(PROFILE if col['geometry']=='profile' else STAKE)[IDS.index(i)];t,z,T=read(i);inside=(z>z[0])&(z<d);zz=z[inside];initz=np.r_[z[0],zz,d];init=np.r_[T[0,0],T[0,inside],0.]
 p=PreparedSlab(t-t[0],T[:,0],np.zeros(len(t)),zz,initz,init,modes=128,domain_provenance=PROVENANCE);a=col['alpha']['pde'];exact=p.bottom_interval_gradient(a);direct=p.bottom_downward_gradient(a).gradient_K_m
 Q=float(np.sum(exact.integral_K_s_m));trap=float(np.trapz(direct,t));qs=float(np.trapz(T[:,0]/(d-z[0]),t));M=exact.temperature_moment_K_m
 rows.append(dict(profile=i,geometry=col['geometry'],exact_signed_unit_integral=Q,instantaneous_trapezoid=trap,relative_operator_difference=(trap-Q)/Q,quasi_steady_signed_integral=qs,storage_correction=(float(M[-1])-float(M[0]))/a,moment_identity_residual=Q-qs+(float(M[-1])-float(M[0]))/a,negative_sample_count=int(np.sum(direct<0))))
Path(__file__).with_name('SIGNED_ENERGY_AUDIT.json').write_text(json.dumps(rows,indent=2));print(json.dumps(rows,indent=2))
