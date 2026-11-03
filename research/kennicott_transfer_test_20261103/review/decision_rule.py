"""Pre-release engineering rule, not a significance test. No file/target access."""
import math

def evaluate(metrics,primary='physical_layered_local',closest='homogeneous_empirical_scale'):
 if primary not in metrics or closest not in metrics:return {'flag':False,'reason':'Primary or closest-standard comparator unscorable; no meaningful-gain claim'}
 def check(m):
  es=m['absolute_errors_mmwe'];v=m['MAE_mmwe']
  if len(es)!=3 or any(not math.isfinite(x)or x<0 for x in es)or not math.isfinite(v)or v<0:raise ValueError('Three finite nonnegative errors and MAE required')
  if abs(sum(es)/3-v)>1e-6*max(1,v):raise ValueError('MAE does not match all three frozen stakes')
 for m in metrics.values():check(m)
 p=metrics[primary];others={n:x for n,x in metrics.items()if n!=primary}
 if not others:return {'flag':False,'reason':'No scorable baselines; refusals are not wins'}
 checks={n:{'MAE_mmwe':x['MAE_mmwe'],'at_least_10percent_and_strictly_lower':p['MAE_mmwe']<=.9*x['MAE_mmwe']and p['MAE_mmwe']<x['MAE_mmwe']}for n,x in others.items()}
 wins=sum(a<b for a,b in zip(p['absolute_errors_mmwe'],metrics[closest]['absolute_errors_mmwe']))
 return {'flag':all(x['at_least_10percent_and_strictly_lower']for x in checks.values())and wins>=2,'all_scorable_baseline_checks':checks,'strict_stake_wins_vs_closest':wins,'required_stake_wins':2,'n_stakes':3,'interpretation':'Predeclared engineering decision only; one glacier/year, no statistical confidence claim','unscorable_methods':'Must be listed separately; never counted as wins'}
