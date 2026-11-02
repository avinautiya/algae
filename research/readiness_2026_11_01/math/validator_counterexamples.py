"""Check the exact committed validator in isolation; no repository mutation."""
import json,copy
from pathlib import Path
from original_validity_contract import check,REQUIRED
base={g:{k:True for k in ks} for g,ks in REQUIRED.items()}
cases={}
def add(label,group,key,value):
 x=copy.deepcopy(base);x[group][key]=value;cases[label]=dict(result=check(x),edited_field=f'{group}.{key}',edited_value_repr=repr(value))
add('false_string_flag','measurement','independent_of_calibration_data','false')
add('support_fraction_above_one','temperature','observed_support_fraction',2.)
add('negative_uncertainty','geometry','floor_survey_sigma_m',-.1)
add('nan_uncertainty','discharge','discharge_uncertainty',float('nan'))
add('empty_source_dictionary','geometry','cross_section_width_source',{})
add('empty_survey_dates','geometry','floor_elevation_survey_dates',[])
cases['all_boolean_placeholders_no_numeric_or_source_evidence']=dict(result=check(base))
assert all(x['result']['may_issue_environmental_prediction'] for x in cases.values())
out=dict(pinned_commit='05e7312e3bbf248acf9f646267a5b8ea06a8b32e',scope='demonstrated schema/attestation gate failures; no environmental claims',numerics_group_in_code='numerics' in REQUIRED,cases=cases)
Path(__file__).with_name('VALIDATOR_COUNTEREXAMPLES.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
print(json.dumps(out,indent=2))
