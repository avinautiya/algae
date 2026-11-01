"""Run the four additive repair suites; never launch the finest campaign."""
from pathlib import Path
import subprocess
import sys


if __name__=='__main__':
    root=Path(__file__).resolve().parent
    commands=[(root,[sys.executable,'-m','unittest','discover','-s','melt','-t','.','-p','test_*.py','-v']),
              (root/'classification',[sys.executable,'-m','unittest','-v','test_recovery.py']),
              (root/'numerics',[sys.executable,'-m','unittest','-v','test_uncertainty.py','test_finest_check.py']),
              (root,[sys.executable,'-m','unittest','discover','-s','bed_audit','-t','.','-p','test_*.py','-v'])]
    for cwd,cmd in commands:subprocess.run(cmd,cwd=cwd,check=True)
    print('All four repair suites passed; no production or finest-grid case launched.')
