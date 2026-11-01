"""Run: python3 -m dynamic.cli input.json --output forecast.json"""
import argparse
import json
import sys
import os
import tempfile
from pathlib import Path
from .model import forecast


def unique_pairs(pairs):
    result={}
    for key,value in pairs:
        if key in result:raise ValueError('duplicate JSON key: '+key)
        result[key]=value
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('input',type=Path);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    try:
        if a.input.resolve()==a.output.resolve():raise ValueError('output must not overwrite input')
        data=json.loads(a.input.read_text(),object_pairs_hook=unique_pairs,
                        parse_constant=lambda x:(_ for _ in ()).throw(ValueError('nonstandard JSON '+x)))
        result=forecast(data)
        content=json.dumps(result,indent=2,allow_nan=False)+'\n'
        with tempfile.NamedTemporaryFile(mode='w',dir=a.output.parent,delete=False) as temporary:
            temporary.write(content);name=temporary.name
        try:os.replace(name,a.output)
        finally:
            if os.path.exists(name):os.unlink(name)
    except (ValueError,RuntimeError,OSError,TypeError,OverflowError) as exc:
        print('forecast refused: '+str(exc),file=sys.stderr);return 2
    print(json.dumps({'output':str(a.output),'samples':len(result['series']),'checks':result['checks']}));return 0


if __name__=='__main__':sys.exit(main())
