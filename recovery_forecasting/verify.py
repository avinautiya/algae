"""Run all delivery tests in isolated subprocesses, fail on any failure."""
from pathlib import Path
import subprocess
import sys

root=Path(__file__).resolve().parent
commands=[
    ([sys.executable,'-m','unittest','discover','-s','dynamic','-t','.','-p','test_*.py','-v'],root),
    ([sys.executable,'-m','unittest','discover','-s','audit','-t','.','-p','test_*.py','-v'],root),
    ([sys.executable,'-m','unittest','-v','test_forecast.py'],root/'predictor'),
    ([sys.executable,'-m','unittest','-v','test_spatial_forecast.py'],root/'spatial'),
]
if __name__=='__main__':
    for command,cwd in commands:
        subprocess.run(command,cwd=cwd,check=True)
    print('All four delivery test suites passed.')
