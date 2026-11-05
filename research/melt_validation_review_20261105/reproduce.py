"""Recompute in an isolated temporary directory, leaving released files intact."""
from pathlib import Path
import shutil,tempfile,subprocess,sys,json
R=Path(__file__).resolve().parent

def main():
    source=R/'followup/Debris_properties_estimations.zip'
    if not source.is_file():raise SystemExit('Run download_source.py first')
    with tempfile.TemporaryDirectory(prefix='melt-release-reproduce-') as d:
        copy=Path(d)/'followup';shutil.copytree(R/'followup',copy)
        for name in ['ELIGIBILITY.json','RESULTS.json','PUBLISHED_RESULTS.json','AUDIT.json']:(copy/'piramide_test'/name).unlink()
        for script in ['run_test.py','published_comparison.py','audit.py']:
            result=subprocess.run([sys.executable,str(copy/'piramide_test'/script)],text=True,capture_output=True)
            if result.returncode:raise RuntimeError(result.stdout+result.stderr)
        for name in ['RESULTS.json','PUBLISHED_RESULTS.json']:
            fresh=json.loads((copy/'piramide_test'/name).read_text());saved=json.loads((R/'followup/piramide_test'/name).read_text())
            assert fresh['scores']==saved['scores'], 'Different score tables: '+name
        print('PASS: raw-source scoring, published comparisons and audit reproduce both saved score tables; released files unchanged.')

if __name__=='__main__':main()
