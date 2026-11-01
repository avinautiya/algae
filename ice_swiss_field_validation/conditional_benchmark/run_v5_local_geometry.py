"""Frozen native-gap forcing reanalysis. Thermal model definitions unchanged."""
import argparse
import datetime as dt
from functools import lru_cache
import json
from pathlib import Path
import sys
import run_local_geometry_development as comparison

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'temperature_digitization'))
from gap_reconstructed_temperature import GapReconstructedTemperatureSeries

METHODS={'linear':'linear_endpoint_bridge','left':'constant_previous_endpoint','right':'constant_next_endpoint'}

class Adapter:
    def __init__(self,series,method):self.series=series;self.method=method
    def integrate_positive(self,a,b,extra_offset_C=0):
        return self.series.integrate_positive(a,b,extra_offset_C=extra_offset_C,method=self.method)

def main(source,out):
    freeze=json.loads((ROOT/'FREEZE_V5.json').read_text())
    if freeze['protocol_sha256']!=comparison.filehash(ROOT/'PROTOCOL_V5_LOCAL_GEOMETRY.md'):
        raise ValueError('V4 protocol changed')
    for relative,expected in freeze['new_inputs_verified_sha256'].items():
        if comparison.filehash(source/relative)!=expected:raise ValueError('V5 source input changed')
    series={s:GapReconstructedTemperatureSeries.from_archive(source,s,max_observed_gap_seconds=60,max_reconstruction_gap_seconds=86400) for s in [3,5]}
    membership={};summaries={}
    for label,method in METHODS.items():
        destination=out/label
        comparison.run(source,destination,temperature_factory=lambda root,s:Adapter(series[s],method))
        path=destination/'CONDITIONAL_RESULTS.json';result=json.loads(path.read_text())
        result.update(v5_freeze=freeze,v5_execution_code_sha256=comparison.filehash(Path(__file__)),
            gap_method=method,selection='linear primary from temperature-only mask experiment; left/right sensitivity; no melt-based choice',
            reconstruction='local rectangular D/Re from imposed daily station depth and common Q; native-time Tw plus <=24h endpoint gaps; target-exposed development reanalysis',
            gap_limits=dict(native_max_segment_seconds=60,reconstruction_max_gap_seconds=86400))
        for name,scenario in result['scenarios'].items():
            membership[label,name]=tuple(r['id'] for r in scenario.get('accepted_intervals',[]))
            for row in scenario.get('accepted_intervals',[]):
                a=comparison.timestamp(row['start'])+dt.timedelta(hours=scenario['clock_shift_h'])
                b=comparison.timestamp(row['end'])+dt.timedelta(hours=scenario['clock_shift_h'])
                row['temperature_support']=series[int(row['station'][1:])].integrate_positive(a,b,extra_offset_C=scenario['temperature_shift_C'],method=method)
        path.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
        nominal=result['scenarios']['clock+0h_T+0.00C_W2.00m']
        summaries[label]=dict(status=nominal['status'],calibration=nominal.get('calibration'),
                              test_metrics=nominal.get('test_metrics'),excluded_intervals=nominal.get('excluded_intervals'))
    for name in result['scenarios']:
        if not membership['linear',name]==membership['left',name]==membership['right',name]:
            raise AssertionError('method-dependent availability requires common-case reanalysis before comparison')
    out.mkdir(parents=True,exist_ok=True)
    (out/'V5_SUMMARY.json').write_text(json.dumps(dict(freeze=freeze,availability_identical_across_imputers=True,nominal=summaries),indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(availability_identical_across_imputers=True,nominal=summaries),indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('source',type=Path);p.add_argument('--out',type=Path,required=True);a=p.parse_args();main(a.source,a.out)
