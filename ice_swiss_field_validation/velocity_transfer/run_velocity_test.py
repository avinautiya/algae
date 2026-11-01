"""Small retrospective transfer test; externally imposed geometry, not melt validation."""
import csv
from datetime import datetime
import hashlib
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parent

def diameter(width, up_depth, down_depth):
    if not all(math.isfinite(v) and v>0 for v in [width, up_depth, down_depth]):
        raise ValueError('Finite positive width and depths required')
    depth=(up_depth+down_depth)/2
    return 4*width*depth/(width+2*depth)

def velocity_basis(head, hydraulic_diameter, length):
    if not all(math.isfinite(v) and v>0 for v in [head,hydraulic_diameter,length]):
        raise ValueError('Finite positive head, diameter, length required')
    return math.sqrt(2*9.81*hydraulic_diameter*head/length)

def fit_amplitude(basis, training_targets):
    if len(basis)!=len(training_targets) or not basis:
        raise ValueError('Matched nonempty training arrays required')
    if not all(math.isfinite(v) and v>0 for v in basis+training_targets):
        raise ValueError('Positive finite training values required')
    return math.fsum(x*y for x,y in zip(basis,training_targets))/math.fsum(x*x for x in basis)

def metrics(observed,predicted):
    errors=[p-y for p,y in zip(predicted,observed)]
    return dict(MAE_m_s=math.fsum(abs(e) for e in errors)/len(errors),
                RMSE_m_s=math.sqrt(math.fsum(e*e for e in errors)/len(errors)),
                bias_m_s=math.fsum(errors)/len(errors), cases=len(errors))

def source_date(text):
    # fromisoformat accepts source unpadded hours inconsistently across Python;
    # only date is used for the predeclared same-day floor matching.
    return text[:10]

def main():
    protocol_bytes=(ROOT/'PROTOCOL.json').read_bytes()
    protocol=json.loads(protocol_bytes)
    original=ROOT.parent/'original_mirror'
    floors={}
    for station in [3,5]:
        path=original/f'data/Channel_elevation/P{station}_channel_elevation.csv'
        with path.open() as f:
            floors[station]={r[0][:10]:float(r[1]) for r in list(csv.reader(f,delimiter=';'))[1:]}
    with (ROOT.parent/'P5_P3_tracer_velocity_printed.csv').open() as f:
        velocities={r['time_at_downstream_station_source_naive']:float(r['velocity_m_s']) for r in csv.DictReader(f)}
    times=protocol['train_source_times']+protocol['test_source_times']
    records=[]
    for i,t in enumerate(times):
        d=t[:10]
        zu,zd=floors[5][d],floors[3][d]
        hu,hd=protocol['P5_depths_m_in_order'][i],protocol['P3_depths_m_in_order'][i]
        records.append(dict(time_source_naive=t,split='train' if i<3 else 'test',observed_velocity_m_s=velocities[t],
                            z5_m=zu,z3_m=zd,h5_m=hu,h3_m=hd,piezometric_head_m=(zu+hu)-(zd+hd),
                            archived_plus_head_m=zu+hu-zd+hd,
                            P5_depth_status='author deduced' if d=='2019-07-11' else 'manual observed',
                            P3_depth_status='pressure-sensor observation' if d=='2019-07-31' else 'nearest same-day manual observation'))
    obs=[r['observed_velocity_m_s'] for r in records]
    constant=math.fsum(obs[:3])/3
    results={};flat=[]
    for width in protocol['width_scenarios_m']:
        candidates={};parameters={}
        for expression in ['piezometric','archived_plus']:
            basis=[velocity_basis(r[expression+'_head_m'],diameter(width,r['h5_m'],r['h3_m']),protocol['distance_m']) for r in records]
            amp=fit_amplitude(basis[:3],obs[:3])
            parameters[expression]=dict(amplitude=amp,Darcy_f=1/amp**2,training_objective_velocity_SSE=math.fsum((amp*x-y)**2 for x,y in zip(basis[:3],obs[:3])))
            candidates[expression+'_fit']=[amp*x for x in basis]
            candidates[expression+'_fixed_f03']=[x/math.sqrt(.3) for x in basis]
        candidates['training_constant_velocity']=[constant]*len(records)
        scored={key:metrics(obs[3:],pred[3:]) for key,pred in candidates.items()}
        primary=scored['piezometric_fit'];comparator=scored['archived_plus_fit']
        passes=primary['MAE_m_s']<comparator['MAE_m_s'] and primary['RMSE_m_s']<comparator['RMSE_m_s']
        passes=passes and all(abs(p-y)<abs(q-y) for p,q,y in zip(candidates['piezometric_fit'][3:],candidates['archived_plus_fit'][3:],obs[3:]))
        results[str(width)]=dict(fitted_parameters=parameters,training_constant_velocity_m_s=constant,test_metrics=scored,
                                 paired_MAE_reduction_piezometric_vs_archived_plus=comparator['MAE_m_s']-primary['MAE_m_s'],consistent_superiority_within_width=passes)
        for i,record in enumerate(records):
            for key,pred in candidates.items():
                flat.append(dict(record,width_m=width,model=key,predicted_velocity_m_s=pred[i],signed_error_m_s=pred[i]-obs[i],absolute_error_m_s=abs(pred[i]-obs[i])))
    outcome=dict(protocol_sha256=hashlib.sha256(protocol_bytes).hexdigest(),code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                 scope=protocol['scope'],status='Retrospective conditional comparison; source targets previously inspected',
                 all_width_consistent_superiority=all(v['consistent_superiority_within_width'] for v in results.values()),widths=results)
    (ROOT/'RESULTS.json').write_text(json.dumps(outcome,indent=2)+'\n')
    with (ROOT/'CASE_PREDICTIONS.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(flat[0]));writer.writeheader();writer.writerows(flat)
    print(json.dumps(outcome,indent=2))

if __name__=='__main__':main()
