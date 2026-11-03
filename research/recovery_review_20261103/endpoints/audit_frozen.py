"""Audit pinned historical labels without retroactively changing their rules."""
import hashlib,json
from collections import Counter
from pathlib import Path
ROOT=Path(__file__).resolve().parent
def main():
 manifest=json.loads((ROOT/'SOURCE_MANIFEST.json').read_text())
 for f in manifest['files']:
  b=(ROOT/'source'/f['local']).read_bytes();actual=hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()
  if actual!=f['git_blob_sha']:raise ValueError('pinned source mismatch')
 results=json.loads((ROOT/'source/analytic_lane__analytic_lane_results.json').read_text())
 cases=results['cases'];evolving=[c for c in cases if c['control']=='evolving'];labelmatches=0
 for c in evolving:
  expected=['initially_disconnected'] if c['pulse_multiple']==1 else ['temporary_or_incomplete_reconnection','initially_disconnected']
  labelmatches+=c['classification'] in expected
 report=dict(pinned_commit=manifest['commit'],verified_source_blobs=4,n_evolving=len(evolving),
  historical_classes=dict(Counter(c['classification'] for c in evolving)),
  frozen_primary_NOT_sustained_matches=sum(not c['sustained_postpulse'] for c in evolving),
  frozen_finer_label_matches=labelmatches,frozen_finer_label_mismatches=len(evolving)-labelmatches,
  stored_history_digest_nulls=sum(c['prefix_sha256'] is None for c in cases),stored_history_digest_total=len(cases),
  runner_wrong_history_keys=['preloss_prefix_sha256','fingerprint'],actual_v1_key='sampled_preloss_history_sha256',
  no_retroactive_wins=True,
  interpretation='Primary no-sustained prediction succeeds; finer label prediction mismatches27/36. Historical v1 is internally consistent, but its labels are coarse for initially disconnected pulse access. Null history IDs cannot establish independence. V2 is prospective only.')
 if report['n_evolving']!=36 or report['frozen_primary_NOT_sustained_matches']!=36 or report['frozen_finer_label_mismatches']!=27 or report['stored_history_digest_nulls']!=72:raise ValueError('unexpected historical evidence; review new source, never force expected result')
 (ROOT/'AUDIT_RESULT.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
