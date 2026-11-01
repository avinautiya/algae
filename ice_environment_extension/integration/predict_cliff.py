"""Executable conditional cliff-volume predictor from net latent-energy packets.

Run: python predict_cliff.py INPUT.json OUTPUT.json
Input energies must already be net of cold-ice warming and debris heat storage.
This geometry driver does not synthesize missing Himalayan flow or weather.
"""
import argparse
from dataclasses import asdict
import hashlib
import json
import math
from pathlib import Path
from stepped_cliff import SteppedCliff, advance_latent, face_partition


def predict(config):
    if not isinstance(config,dict) or set(config)!={'geometry','intervals','energy_provenance'}:
        raise ValueError('exact geometry, intervals and energy_provenance fields required')
    if not isinstance(config['energy_provenance'],str) or not config['energy_provenance'].strip():
        raise ValueError('describe energy source and cold/storage deductions')
    cliff=SteppedCliff(**config['geometry']).validate();initial=cliff
    if not isinstance(config['intervals'],list) or not config['intervals']:raise ValueError('nonempty interval list required')
    rows=[];t=0.;face=0.;floor=0.;energy=0.
    for interval in config['intervals']:
        required={'end_time_s','face_net_latent_J','floor_net_latent_J','water_stage_m','debris_top_m'}
        if not isinstance(interval,dict) or set(interval)!=required:raise ValueError('interval keys invalid')
        end=interval['end_time_s']
        if isinstance(end,bool) or not isinstance(end,(int,float)) or not math.isfinite(end) or end<=t:
            raise ValueError('strictly increasing finite end times required')
        cliff,ledger=advance_latent(cliff,interval['face_net_latent_J'],interval['floor_net_latent_J'])
        partition=face_partition(cliff,interval['water_stage_m'],interval['debris_top_m'])
        face+=ledger['face_melt_m3'];floor+=ledger['floor_melt_m3'];energy+=ledger['phase_energy_J']
        rows.append(dict(time_s=end,geometry=asdict(cliff),**ledger,**partition,
                         cumulative_face_melt_m3=face,cumulative_floor_melt_m3=floor))
        t=end
    actual=initial.ice_volume_m3-cliff.ice_volume_m3
    return dict(input_sha256=hashlib.sha256(json.dumps(config,sort_keys=True,allow_nan=False).encode()).hexdigest(),
                evidence='conditional planar-step forecast; not validated Himalayan cliff survival',
                energy_provenance=config['energy_provenance'],series=rows,
                cumulative=dict(face_melt_m3=face,floor_melt_m3=floor,total_melt_m3=face+floor,
                                polygon_volume_loss_m3=actual,volume_residual_m3=actual-face-floor,
                                phase_energy_J=energy),
                scope=['uniform normal face retreat and horizontal floor lowering',
                       'fixed crest and extrusion width; no overhang or collapse',
                       'prescribed stage, debris elevation and net latent energies',
                       'no atmospheric face energy assigned to floor incision'])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('input',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();result=predict(json.loads(a.input.read_text()))
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result['cumulative'],indent=2))
