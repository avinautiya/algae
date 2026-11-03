"""Predict calibrated heat/melt on an explicitly supplied supported linear branch.
Usage: python predict.py input.json --output NEW_FILE.json
This does not infer geometry, units or phase validity from absent metadata.
"""
import argparse,json
from pathlib import Path
from ratio_enclosure import enclosure,coefficient_box_enclosure

def predict(config):
 if config.get('units')!='mm_water_equivalent' or config.get('linear_phase_branch_justified')is not True:
  raise ValueError('Require explicit water-equivalent units and justified linear phase branch')
 if config.get('source_feasible_set_justified')is not True:raise ValueError('Unsupported feasible set: abstain before computing')
 keys=['A','b','bounds','training_target','tolerance'];kwargs={k:config[k]for k in keys}
 if config.get('coefficient_boxes'):
  kwargs.update(config['coefficient_boxes']);result=coefficient_box_enclosure(**kwargs)
 else:
  kwargs.update(train=config['train'],future=config['future']);result=enclosure(**kwargs)
 result['prediction_type']='Conditional calibrated melt/heat inference; not automatically an advance environmental forecast'
 return result
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('input',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 config=json.loads(a.input.read_text());result=predict(config)
 with a.output.open('x')as out:json.dump(result,out,indent=2,allow_nan=False);out.write('\n')
