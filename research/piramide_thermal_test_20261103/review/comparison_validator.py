"""Fail-closed comparison specification and conservative clustered win audit.

Validates declarations, not their historical truth. No model fitting or target
loading. A sign test concerns independent block-wise wins, not universal skill.
"""
import hashlib,json,math
from math import comb
SHA=lambda s:isinstance(s,str) and len(s)==64 and all(c in '0123456789abcdef' for c in s)
def validate_protocol(p):
 errors=[]
 def require(ok,message):
  if not ok:errors.append(message)
 require(p.get('version')==1,'version1 required')
 candidate=p.get('candidate',{});cid=candidate.get('id')
 require(isinstance(cid,str) and bool(cid),'one candidate id required')
 require(candidate.get('selected_using')=='development_only','candidate selection must use development only')
 require(isinstance(candidate.get('all_considered_candidates'),list) and cid in candidate.get('all_considered_candidates',[]),'complete hypothesis ledger required')
 require(bool(candidate.get('selection_rule')),'development selection rule required')
 require(SHA(candidate.get('code_sha256')),'candidate source hash required')
 splits=p.get('splits',{});dev=set(splits.get('development_groups',[]));ev=set(splits.get('evaluation_groups',[]))
 require(bool(dev) and bool(ev) and not dev&ev,'nonempty disjoint development/evaluation groups required')
 require(not ev&set(p.get('targets_inspected_before_freeze',[])),'previously inspected targets cannot be confirmatory evaluation')
 require(p.get('prediction_freeze_before_target_scoring') is True,'predictions/rules must be frozen before scoring')
 for par in p.get('parameters',[]):
  uses=set(par.get('all_transitive_melt_dependencies',[]))
  require(par.get('dependency_audit_complete') is True,'parameter '+str(par.get('id'))+': full dependency audit required')
  require(not uses&ev,'parameter '+str(par.get('id'))+': evaluation melt leakage, including pooled published constants')
  require(uses<=dev,'parameter '+str(par.get('id'))+': unknown melt dependency')
 require(bool(p.get('parameters')),'parameter provenance ledger required')
 for d in p.get('datasets',[]):require(bool(d.get('version')) and SHA(d.get('sha256')) and bool(d.get('source_url')),'versioned hashed dataset source required')
 require(bool(p.get('datasets')),'dataset manifest required')
 comparison=p.get('comparison',{})
 require(comparison.get('target')=='native_ice_ablation_mmwe','native independent ice target required')
 for key in ['same_spatial_footprint','same_time_support','same_target_quantity','observation_operator_verified','measurement_uncertainty_declared','all_prespecified_baselines_reported']:
  require(comparison.get(key) is True,key+' must be verified')
 baselines=comparison.get('baselines',[])
 require(bool(baselines) and len(baselines)==len(set(baselines)) and cid not in baselines,'distinct fixed baselines required')
 require(comparison.get('primary_loss')=='block_mean_absolute_error','fixed blockwise MAE required')
 require(type(comparison.get('minimum_practical_improvement_mmwe')) in (int,float) and math.isfinite(comparison.get('minimum_practical_improvement_mmwe',float('nan'))) and comparison.get('minimum_practical_improvement_mmwe',-1)>=0,'prespecified practical improvement margin required')
 require(comparison.get('same_information_budget') is True,'equal observed-input budget required; otherwise separate conditional-reconstruction comparison')
 baseline_sources=comparison.get('baseline_sources',{})
 for b in baselines:
  r=baseline_sources.get(b,{})
  require(bool(r.get('source')) and SHA(r.get('code_or_output_sha256')) and r.get('native_ice_component_verified') is True,'baseline '+str(b)+': source/hash/native-ice mapping required')
 require(comparison.get('multiple_comparison_rule')=='holm_one_sided_sign','fixed Holm family required')
 blocks=comparison.get('independence_blocks',{})
 require(set(blocks)==ev and all(isinstance(v,str) and v for v in blocks.values()),'every evaluation group requires an independence block')
 if p.get('claim_scope')=='cross_glacier_generalization':
  require(comparison.get('block_definition')=='glacier','cross-glacier claims require glacier-level aggregation')
  require(len(set(blocks.values()))>=2,'multiple independent glaciers required')
 elif p.get('claim_scope')!='retrospective_event':require(False,'supported explicit claim scope required')
 require(comparison.get('independence_justification_reviewed') is True,'independent block assumptions must be justified')
 return dict(ready=not errors,errors=errors,declaration_verification_only=True,guaranteed_success=False)

def sign_tail(wins,n):
 if type(wins) is not int or type(n) is not int or not 0<=wins<=n or n<1:raise ValueError('valid nonempty block counts required')
 return sum(comb(n,k) for k in range(wins,n+1))/2**n

def holm(pvalues,alpha=.05):
 if not 0<alpha<1 or not pvalues or any(not math.isfinite(p) or not 0<=p<=1 for p in pvalues.values()):raise ValueError('invalid test family')
 ordered=sorted(pvalues,key=pvalues.get);out={};running=0.
 for j,key in enumerate(ordered):
  running=max(running,min(1.,(len(ordered)-j)*pvalues[key]))
  out[key]=dict(raw_p=pvalues[key],holm_adjusted_p=running,reject=running<=alpha)
 return out

def conservative_block_test(block_rows,baselines,alpha=.05):
 """Bounds must already include joint clocks/reading errors on matched support.

 Independent blocks, not individual camera frames, are the units. Every fixed
 block remains in n; ties/uncertain differences count as non-wins. An upper bound
 <0 certifies the candidate's lower loss on that block under the declared bounds.
 """
 if not block_rows or not baselines:raise ValueError('nonempty fixed evaluation family required')
 ids=[r['block'] for r in block_rows]
 if len(ids)!=len(set(ids)):raise ValueError('aggregate correlated intervals to one row per block first')
 pvalues={};wins={}
 for b in baselines:
  count=0
  for row in block_rows:
   bounds=row['candidate_minus_baseline_loss_bounds_mmwe'][b]
   if len(bounds)!=2 or not all(type(v) in (float,int) and math.isfinite(v) for v in bounds) or bounds[0]>bounds[1]:raise ValueError('finite ordered joint-support bounds required')
   count+=bounds[1]<0
  wins[b]=count;pvalues[b]=sign_tail(count,len(block_rows))
 tests=holm(pvalues,alpha)
 return dict(independent_blocks=len(block_rows),certified_block_wins=wins,tests=tests,
  all_baselines_pass=all(t['reject'] for t in tests.values()),
  statistical_assumptions='independent prespecified blocks; under each null conditional block-win chance is at most1/2',
  uncertainty_semantics='external bounded-support assumptions; not inferred from these p values',
  universal_superiority_established=False)

def write_freeze(protocol,path):
 from pathlib import Path
 path=Path(path)
 if path.exists():raise ValueError('refuse overwriting freeze')
 result=validate_protocol(protocol)
 if not result['ready']:raise ValueError('; '.join(result['errors']))
 canonical=json.dumps(protocol,sort_keys=True,separators=(',',':'),allow_nan=False).encode()
 path.write_text(json.dumps(dict(protocol=protocol,protocol_sha256=hashlib.sha256(canonical).hexdigest()),indent=2)+'\n')
 return hashlib.sha256(canonical).hexdigest()
