"""CLI: python -m ice_cliff_survival <command>
  data-download   pinned, checksum-verified download of the two Zenodo archives (needs network)
  data-audit      audit the archives -> data_audit/AUDIT.json (needs the raw archives)
  verify          run the component verification suite (no data, no CFD)
  demo            SYNTHETIC toe-stock demonstration with illustrative coefficients (never an environmental forecast)
  run             environmental run: NOT SUPPORTED in milestone 0 (exits with an explicit error)"""
import json, os, subprocess, sys, unittest
HERE = os.path.dirname(os.path.abspath(__file__))
def main(argv):
    if len(argv) < 2: print(__doc__); return 2
    c = argv[1]
    if c == "data-download": return subprocess.call([sys.executable, os.path.join(HERE, "data_audit", "download.py")])
    if c == "data-audit": return subprocess.call([sys.executable, os.path.join(HERE, "data_audit", "audit.py")])
    if c == "verify":
        suite = unittest.defaultTestLoader.discover(os.path.join(HERE, "tests"), top_level_dir=os.path.dirname(HERE)); r = unittest.TextTestRunner(verbosity=1).run(suite); return 0 if r.wasSuccessful() else 1
    if c == "demo":
        from .demo import main as dm; return dm(argv[2:])
    if c == "run":
        from . import UnsupportedRun
        raise UnsupportedRun("environmental ice-cliff survival runs are not supported: the coupled flow-heat / moving-ice / grain transport model is not implemented or verified (see ice_cliff_survival/INTEGRATION_CONTRACT.md)")
    print(__doc__); return 2
if __name__ == "__main__": sys.exit(main(sys.argv))
