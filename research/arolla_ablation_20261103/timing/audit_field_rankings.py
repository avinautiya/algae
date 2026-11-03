"""Describe rank identifiability for ALL frozen native H1D timing scenarios."""
import argparse,collections,itertools,json
from pathlib import Path
from ranking_bounds import rank_bound
MODELS=['constant_rate','positive_degree_day','net_shortwave_dose']
def audit(data):
 groups={}
 for r in data['rows']:
  if r['stake']=='H1D' and r['model'] in MODELS:
   key=(r['mean_label_end_offset_h'],tuple(r['visit_hours']))
   if r['model'] in groups.setdefault(key,{}):raise ValueError('duplicate model on fixed timing support')
   groups[key][r['model']]=r
 if len(groups)!=81 or any(set(g)!=set(MODELS) for g in groups.values()):raise ValueError('incomplete fixed native H1D support')
 rows=[]
 for key,g in sorted(groups.items()):
  for a,b in itertools.combinations(MODELS,2):
   ra,rb=g[a],g[b]
   if ra['calibration_observed_mmwe']!=rb['calibration_observed_mmwe'] or ra['heldout_observed_mmwe']!=rb['heldout_observed_mmwe']:raise ValueError('unmatched native observations')
   for association in ['each_shared_reading','each_interval']:
    # Source0.1m × native900kg/m³ ->90mmwater equivalent. Interpretation
    # as absolute bounds is a sensitivity, not an author-certified interval.
    result=rank_bound(ra['calibration_observed_mmwe'],ra['heldout_observed_mmwe'],ra['test_train_dose_ratio'],rb['test_train_dose_ratio'],90.,association=association)
    rows.append(dict(timing=[key[0],list(key[1])],model_a=a,model_b=b,**result))
 counts={}
 for association in ['each_shared_reading','each_interval']:
  counts[association]=dict(collections.Counter(r['status'] for r in rows if r['uncertainty_association']==association))
 return dict(rows=rows,status_counts=counts,source_uncertainty_association_unspecified=True,probabilistic_confidence_claim=False,interpretation='Conditional bounded-error ranking audit. Nominal field errors remain reported; this does not manufacture precision or certify independent glaciers.')
def main():
 p=argparse.ArgumentParser();p.add_argument('results',type=Path);p.add_argument('out',type=Path);a=p.parse_args()
 if a.out.exists():raise ValueError('output exists')
 result=audit(json.loads(a.results.read_text()));a.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result['status_counts'],indent=2))
if __name__=='__main__':main()
