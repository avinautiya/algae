import argparse
import json
from pathlib import Path
from .evidence import readiness
from .evaluate import evaluate

p=argparse.ArgumentParser(description='Self-reported evidence checklist or descriptive forecast comparison')
p.add_argument('mode',choices=('readiness','evaluate'));p.add_argument('input',type=Path);p.add_argument('--output',type=Path)
a=p.parse_args()
try:
    record=json.loads(a.input.read_text(),parse_constant=lambda x: (_ for _ in ()).throw(ValueError('nonstandard JSON constant '+x)))
    result=(readiness if a.mode=='readiness' else evaluate)(record)
except (ValueError,TypeError,KeyError,OSError) as e:
    p.exit(2,'Invalid input: '+str(e)+'\n')
text=json.dumps(result,indent=2,allow_nan=False)+'\n'
if a.output:a.output.write_text(text)
else:print(text,end='')
