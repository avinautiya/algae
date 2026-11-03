"""Decode process observations only; never loads camera ablation/calibration files.

Uses public mat-io0.4.2, verifies its converted arrays against raw MCOS
property maps, raw epoch millisecond values, native dimensions/column names.
No timezone is invented; naive native labels remain naive.
"""
import argparse,hashlib,json,importlib.metadata
from pathlib import Path
import numpy as np
from matio import load_from_mat

FAMILIES=('debtemp','sm','trh','wind')
def digest(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def decode(source,out):
    out=Path(out);out.mkdir(exist_ok=True,parents=True)
    converted=load_from_mat(source,raw_data=False,add_table_attrs=True,variable_names=list(FAMILIES))
    raw=load_from_mat(source,raw_data=True,variable_names=list(FAMILIES))
    report={'source_sha256':digest(source),'parser':{'mat-io':importlib.metadata.version('mat-io')},
            'scope':'Process observations only; no melt targets, model fits or scores',
            'native_clock':'No native timezone; do not label timestamps UTC. Same-source relative local-clock analysis only until external clock authority established.',
            'checks':[]}
    for fam in FAMILIES:
        assert converted[fam].dtype.names==('Tower1','Tower2','Tower3')
        for tower in ('Tower1','Tower2','Tower3'):
            frame=converted[fam][tower][0,0]
            obj=raw[fam][tower][0,0];assert obj.classname=='timetable'
            props=obj.properties['any']
            assert props.shape==(1,1)
            def prop(k):return props[k][0,0]
            names=[str(x.item()) for x in prop('varNames').ravel()]
            assert names==list(frame.columns)
            assert int(prop('numRows').item())==len(frame)
            assert int(prop('numVars').item())==len(frame.columns)
            columns=prop('data').ravel()
            for j,c in enumerate(columns):np.testing.assert_equal(c.ravel(),frame.iloc[:,j].to_numpy())
            rt=prop('rowTimes');assert rt.classname=='datetime'
            tz=rt.properties.get('tz',rt.properties.get('tmz'))
            assert np.asarray(tz).size==0,'Nonempty clock zone requires explicitly audited handling'
            epoch=np.asarray(rt.properties['data']).ravel()
            assert not np.iscomplexobj(epoch) and np.isfinite(epoch).all()
            assert np.equal(epoch,np.rint(epoch)).all(),'Do not round fractional milliseconds'
            native_times=epoch.astype('int64').astype('datetime64[ms]')
            np.testing.assert_equal(native_times,frame.index.to_numpy().astype('datetime64[ms]'))
            assert frame.index.is_unique and frame.index.is_monotonic_increasing
            steps=np.diff(epoch)/1000
            name=f'{fam}_{tower}.csv';path=out/name
            frame.to_csv(path,index_label='SOURCE_LOCAL_LABEL',date_format='%Y-%m-%dT%H:%M:%S')
            report['checks'].append({'family':fam,'tower':tower,'rows':len(frame),'columns':names,
                'first_source_local_label':str(frame.index[0]),'last_source_local_label':str(frame.index[-1]),
                'step_seconds_counts':{str(k):int(v) for k,v in zip(*np.unique(steps,return_counts=True))},
                'missing_per_column':{k:int(v) for k,v in frame.isna().sum().items()},
                'minimum_per_column':{k:float(v) for k,v in frame.min().items()},
                'maximum_per_column':{k:float(v) for k,v in frame.max().items()},
                'raw_column_values_equal':True,'raw_epoch_times_equal':True,
                'file':name,'sha256':digest(path)})
    (out/'DECODE_VERIFICATION.json').write_text(json.dumps(report,indent=2)+'\n')
    return report
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('source');a.add_argument('output');args=a.parse_args()
    r=decode(args.source,args.output)
    for c in r['checks']:print(c['family'],c['tower'],c['rows'],c['first_source_local_label'],c['last_source_local_label'])
