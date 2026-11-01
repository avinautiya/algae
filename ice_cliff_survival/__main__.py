"""CLI: python -m ice_cliff_survival <command>
  data-download   pinned, checksum-verified download of the two Zenodo archives (needs network)
  data-audit      audit the archives -> data_audit/AUDIT.json (needs the raw archives)
  verify          run the component verification suite (no data, no CFD)
  demo            SYNTHETIC toe-stock demonstration with illustrative coefficients (never an environmental forecast)
  heat-ice-2d     2-D heat + moving-ice component run (needs the built OpenFOAM app; NOT environmental)
  flow-stefan     flow + heat + moving-ice coupled channel case (needs the built flowStefanFoam; COMPONENT run, NOT environmental)
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
    if c == "heat-ice-2d":
        from .cfd.cli2d import main as m2; return m2(argv[2:])
    if c == "flow-stefan":
        from .cfd.cli_flow import main as mf; return mf(argv[2:])
    if c == "run":
        sys.stderr.write("UnsupportedRun: environmental ice-cliff survival runs are not supported. Implemented/verified components: 2-D heat-moving-ice with static shielding and the laminar rigid-lid flow-heat-moving-ice channel (ice_cliff_survival/cfd); NOT implemented: free-surface flow coupling, mobile-grain transport, burial dynamics, environmental forcing. See INTEGRATION_CONTRACT.md.\n"); return 2
    print(__doc__); return 2
if __name__ == "__main__": sys.exit(main(sys.argv))
