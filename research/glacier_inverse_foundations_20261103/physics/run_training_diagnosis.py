"""Tower1 TRAIN development only: loss/objective/layering/identifiability audit.

No Tower2 camera or new validation scoring. No model is selected using melt.
The already inspected Tower1 TRAIN mass is a posterior diagnostic only.
"""
import csv,datetime as dt,hashlib,json,math,sys
from pathlib import Path
import numpy as np
from layered_column import predict

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from paths import require_piramide
OLD=require_piramide()
START=dt.datetime(2022,12,2,3);END=dt.datetime(2022,12,16,3)
GRID=[.2,.35,.6,1.,1.6,2.5,3.5,4.5]
DEPTHS=np.array([0.,.105,.23,.355,.48]);DOMAIN=.5;LAYER=.25


def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def load_training():
    filename=OLD/'data/decoded/debtemp_Tower1.csv';manifest=json.loads((OLD/'data/decoded/DECODE_VERIFICATION.json').read_text())
    h=next(r['sha256'] for r in manifest['checks'] if r['family']=='debtemp' and r['tower']=='Tower1')
    if digest(filename)!=h:raise ValueError('TRAIN source hash mismatch')
    rows=list(csv.DictReader(filename.open()));rows=[r for r in rows if START<=dt.datetime.fromisoformat(r['SOURCE_LOCAL_LABEL'])<=END]
    t=[dt.datetime.fromisoformat(r['SOURCE_LOCAL_LABEL']) for r in rows]
    if t[0]!=START or t[-1]!=END or any((b-a).total_seconds()!=1800 for a,b in zip(t,t[1:])):raise ValueError('exact contiguous TRAIN support required')
    temperatures=np.array([[float(r['DT'+str(i)]) for i in range(1,6)] for r in rows])
    if not np.isfinite(temperatures).all():raise ValueError('TRAIN missing data refused')
    capacity=json.loads((OLD/'data/TRAINING_CAPACITY_AUDIT.json').read_text())
    return temperatures,capacity,dict(temperature_csv_sha256=h,capacity_audit_sha256=digest(OLD/'data/TRAINING_CAPACITY_AUDIT.json'))


def losses(prediction,observed,epsilon):
    r=prediction-observed
    return dict(unweighted_temperature_MSE_C2=float(np.mean(r*r)),
        depth_uncertainty_standardized_interval_loss=float(np.mean(np.maximum(np.abs(r)-epsilon,0.)**2/epsilon**2)),
        basal_temperature_and_gradient_loss=float(np.mean((r[:,3]/.4)**2+((r[:,2]-r[:,3])/.8)**2)+.1*np.mean((r[:,:2]/.4)**2)))


def run():
    output=ROOT/'TRAINING_DIAGNOSIS.json'
    if output.exists():raise ValueError('refusing to overwrite development diagnosis')
    data,capacity,provenance=load_training();C=capacity['finite_sensible_capacity_scenes']['training_valid_probe_proxy_not_measured_depth_mean']['C_J_m3K']
    observed=data[:-1,1:];gradients=[]
    for j in range(1,5):
        a,b=(j-1,j+1) if j<4 else (3,4)
        gradients.append((data[:-1,b]-data[:-1,a])/(DEPTHS[b]-DEPTHS[a]))
    epsilon=.4+.02*np.abs(np.array(gradients).T)
    plan=dict(scope='Tower1 TRAIN DEVELOPMENT only; not a new validation result',start=START.isoformat(),end=END.isoformat(),
        k_grid_W_mK=GRID,capacity_J_m3K=C,domain_depth_m=DOMAIN,layer_interface_m=LAYER,
        source_geometry='same source-bounded middle20mm scene; assumed static geometry',
        loss_definitions=['plain equal-sensor MSE','engineering .4C + abs(observed vertical TRAIN gradient)*.02m standardized interval loss',
                          'T5/.4 and (T4-T5)/.8 squared, plus0.1 upper-sensor standardized MSE; predeclared heuristic, not claimed optimal'],
        no_parameter_selection_uses_melt=True,no_Tower2_camera_access=True,source_hashes=provenance)
    (ROOT/'STUDY_PLAN.json').write_text(json.dumps(plan,indent=2)+'\n')
    def model(ku,kl,capacity=C,deep_shift=0.):
        z=DEPTHS.copy();z[-1]+=deep_shift
        return predict(k_upper=float(ku),k_lower=float(kl),capacity=float(capacity),depth=DOMAIN,layer_depth=LAYER,
            initial_depths=list(z)+[DOMAIN],initial_temperatures=list(data[0])+[0.],
            surface_samples=list(data[:,0]),sample_depths=list(z[1:]))
    cache={};records=[]
    for ku in GRID:
        for kl in GRID:
            result=model(ku,kl);p=np.array(result['predicted_temperature_C'][:-1]);score=losses(p,observed,epsilon)
            row=dict(k_upper_W_mK=ku,k_lower_W_mK=kl,alpha_upper_m2_s=ku/C,alpha_lower_m2_s=kl/C,
                     raw_TRAIN_melt_mmwe=result['melt_mmwe'],sensor_RMSE_C=list(np.sqrt(np.mean((p-observed)**2,axis=0))),
                     max_standardized_absolute_residual=float(np.max(np.abs(p-observed)/epsilon)),
                     engineering_box_all_rows_compatible=bool(np.all(np.abs(p-observed)<=epsilon)),
                     relative_energy_residual=result['relative_energy_residual'],**score)
            records.append(row);cache[(ku,kl)]=(row,p)
    names=list(losses(observed,observed,epsilon));winners=[]
    for family in ['homogeneous','two_layer']:
        candidates=[r for r in records if family=='two_layer' or r['k_upper_W_mK']==r['k_lower_W_mK']]
        for loss in names:
            win=min(candidates,key=lambda x:(x[loss],x['k_upper_W_mK'],x['k_lower_W_mK']))
            winners.append(dict(family=family,objective=loss,**win,
                selected_at_grid_boundary=win['k_upper_W_mK'] in [GRID[0],GRID[-1]] or win['k_lower_W_mK'] in [GRID[0],GRID[-1]]))
    # All hypotheses chosen by TRAIN temperature objectives. Only now inspect the
    # already known original TRAIN camera scalar as a development diagnostic.
    known=json.loads((OLD/'model/frozen_piramide_v3/melt_comparison/PREDICTIONS.json').read_text())
    observed_train=float(known['calibration']['observed_ice_ablation_mmwe'])
    for w in winners:
        w['already_inspected_TRAIN_camera_mmwe']=observed_train
        w['TRAIN_output_gain_needed_diagnostic_only']=observed_train/w['raw_TRAIN_melt_mmwe']
        w['TRAIN_absolute_melt_discrepancy_mmwe_not_validation']=abs(w['raw_TRAIN_melt_mmwe']-observed_train)
    base=next(w for w in winners if w['family']=='two_layer' and w['objective']==names[0]);ku=base['k_upper_W_mK'];kl=base['k_lower_W_mK'];h=.03
    jac=[];goal=[]
    for coordinate in range(4):
        predictions=[];masses=[]
        for sign in [-1,1]:
            u,l,c=ku,kl,C;shift=0.
            if coordinate==0:u*=math.exp(sign*h)
            elif coordinate==1:l*=math.exp(sign*h)
            elif coordinate==2:c*=math.exp(sign*h)
            else:shift=sign*h*.02
            r=model(u,l,c,shift);predictions.append(np.array(r['predicted_temperature_C'][:-1]));masses.append(r['melt_mmwe'])
        jac.append(((predictions[1]-predictions[0])/(2*h)/epsilon).ravel())
        goal.append((math.log(masses[1])-math.log(masses[0]))/(2*h))
    J=np.array(jac).T/math.sqrt(observed.size);U,s,Vt=np.linalg.svd(J,full_matrices=False);b=np.array(goal)
    sensitivity=dict(reference='unweighted two-layer TRAIN optimum, no melt-selected reference',parameters=['log k_upper','log k_lower','log C','deep probe shift/.02m'],
        normalized_J_singular_values=list(s),condition_number=float(s[0]/s[-1]),
        goal_log_melt_gradient=list(b),right_singular_vectors=Vt.tolist(),
        local_goal_error_amplification=float(np.linalg.norm((Vt@b)/s)),
        interpretation='engineering-error-scaled local RMS sensitivities; not Fisher information, confidence intervals or a global uncertainty bound')
    scales=[]
    for label,scene in capacity['finite_sensible_capacity_scenes'].items():
        c=scene['C_J_m3K'];a=c/C;r=model(ku*a,kl*a,c);p=np.array(r['predicted_temperature_C'][:-1]);reference=cache[(ku,kl)][1]
        scales.append(dict(capacity_scene=label,scale=a,capacity_J_m3K=c,
            alpha_upper_unchanged=True,alpha_lower_unchanged=True,
            max_temperature_difference_C=float(np.max(np.abs(p-reference))),
            RMS_temperature_difference_C=float(np.sqrt(np.mean((p-reference)**2))),
            raw_TRAIN_melt_mmwe=r['melt_mmwe']))
    report=dict(plan_sha256=digest(ROOT/'STUDY_PLAN.json'),all_candidates=records,winners=winners,
        source_sensor_error_C=.4,depth_error_m=.02,number_engineering_box_compatible_candidates=sum(r['engineering_box_all_rows_compatible'] for r in records),
        identifiability=sensitivity,cold_ice_absolute_scale_test=scales,
        Tower2_camera_values_opened=False,new_validation_or_originality_claim=False,
        previous_failed_field_verdicts_unchanged=True)
    output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(output=str(output),candidates=len(records),winners=winners,identifiability=sensitivity,coldscale=scales),indent=2))


if __name__=='__main__':run()
