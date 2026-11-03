"""TRAIN-only supplemental diagnostics; never refits on or scores Tower2 melt."""
import json,math
from pathlib import Path
import numpy as np
from layered_column import predict
from run_training_diagnosis import load_training,DEPTHS,DOMAIN,LAYER

ROOT=Path(__file__).resolve().parent


def run():
    out=ROOT/'ADDITIONAL_IDENTIFIABILITY.json'
    if out.exists():raise ValueError('new diagnostic output required')
    source=json.loads((ROOT/'TRAINING_DIAGNOSIS.json').read_text());data,capacities,_=load_training()
    C=capacities['finite_sensible_capacity_scenes']['training_valid_probe_proxy_not_measured_depth_mean']['C_J_m3K']
    observed=data[:-1,1:];grad=[]
    for j in range(1,5):
        a,b=(j-1,j+1) if j<4 else (3,4)
        grad.append((data[:-1,b]-data[:-1,a])/(DEPTHS[b]-DEPTHS[a]))
    eps=.4+.02*np.abs(np.array(grad).T)
    def model(ku,kl,c=C,shift=0.):
        z=DEPTHS.copy();z[-1]+=shift
        return predict(k_upper=ku,k_lower=kl,capacity=c,depth=DOMAIN,layer_depth=LAYER,
            initial_depths=list(z)+[DOMAIN],initial_temperatures=list(data[0])+[0.],
            surface_samples=list(data[:,0]),sample_depths=list(z[1:]))
    chosen=next(x for x in source['winners'] if x['family']=='two_layer' and x['objective']=='unweighted_temperature_MSE_C2')
    ku=chosen['k_upper_W_mK'];kl=chosen['k_lower_W_mK'];base=model(ku,kl)
    sensitivity=[]
    for h in [.03,.015,.0075]:
        cols=[]
        for j in range(4):
            pp=[]
            for sign in [-1,1]:
                u,l,c=ku,kl,C;shift=0.
                if j==0:u*=math.exp(sign*h)
                elif j==1:l*=math.exp(sign*h)
                elif j==2:c*=math.exp(sign*h)
                else:shift=sign*h*.02
                pp.append(np.array(model(u,l,c,shift)['predicted_temperature_C'][:-1]))
            cols.append(((pp[1]-pp[0])/(2*h)/eps).ravel())
        J=np.array(cols).T/np.sqrt(observed.size);s=np.linalg.svd(J,compute_uv=False)
        sensitivity.append(dict(log_parameter_difference_step=h,singular_values=s.tolist(),
            scaling_null_direction_norm=float(np.linalg.norm(J@np.array([1.,1.,1.,0.])/math.sqrt(3)))))
    scales=[]
    for family,first,second in [('homogeneous',4.5,4.5),('two_layer',ku,kl)]:
        reference=model(first,second);p=np.array(reference['predicted_temperature_C'])
        for label,scene in capacities['finite_sensible_capacity_scenes'].items():
            ratio=scene['C_J_m3K']/C;r=model(first*ratio,second*ratio,scene['C_J_m3K'])
            delta=np.array(r['predicted_temperature_C'])-p
            scales.append(dict(family=family,capacity_scene=label,scale=ratio,
                cold_interface_steps=r['cold_interface_steps'],
                max_temperature_change_C=float(np.max(np.abs(delta))),
                RMS_temperature_change_C=float(np.sqrt(np.mean(delta*delta))),
                raw_TRAIN_melt_mmwe=r['melt_mmwe']))
    local=[]
    for j in [1,2,3]:
        Tt=(data[2:,j]-data[:-2,j])/3600
        a=DEPTHS[j]-DEPTHS[j-1];b=DEPTHS[j+1]-DEPTHS[j]
        Tzz=2*(data[1:-1,j-1]/(a*(a+b))-data[1:-1,j]/(a*b)+data[1:-1,j+1]/(b*(a+b)))
        alpha=float(Tzz@Tt/(Tzz@Tzz))
        local.append(dict(sensor='DT'+str(j+1),derivative_OLS_alpha_m2_s=alpha,
            conditional_k_W_mK=alpha*C,derivative_residual_RMS_K_s=float(np.sqrt(np.mean((Tt-alpha*Tzz)**2))),
            interpretation='noise-sensitive local PDE regression, not independent conductivity measurement'))
    output=dict(two_layer_unweighted_reference_cold_interface_steps=base['cold_interface_steps'],
        finite_difference_step_sensitivity=sensitivity,known_ice_absolute_scale_tests=scales,
        sensor_means_C=np.mean(data[:-1],axis=0).tolist(),local_derivative_diagnostics=local,
        note='Any nearzero singular value scaling as h² is differentiation truncation, not evidence of identifiable absolute heat scale.',
        Tower2_camera_opened=False,production_validation_run=False)
    out.write_text(json.dumps(output,indent=2)+'\n');print(json.dumps(output,indent=2))


if __name__=='__main__':run()
