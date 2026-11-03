"""Process features only; never imports or opens camera objects."""
import json,hashlib,sys
from pathlib import Path
import pandas as pd,numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from paths import require_piramide
sys.path.insert(0,str(require_piramide()/'data'))
from read_verified_process import read_process
root=require_piramide()/'data/decoded';t=read_process(root,'debtemp','Tower1')
train=t.loc[(t.index>='2022-12-02T03:00')&(t.index<'2022-12-16T03:00')]
out={'scope':'Tower1 process features only; no camera/melt observations', 'source_sha256':hashlib.sha256((root/'debtemp_Tower1.csv').read_bytes()).hexdigest(),
 'sensor_absolute_accuracy_C':.4,'not_sigma':True,'independent_offset_calibration_found':False,
 'deepest_sensor_is_not_independently_verified_ice_probe':True,'full_record_rows':len(t),'training_rows':len(train),
 'records':{}}
for name,d in [('predeclared_train',train),('full_released_process_record',t)]:
 grad=d.DT4-d.DT5
 out['records'][name]={'DT5_min_C':float(d.DT5.min()),'DT5_max_C':float(d.DT5.max()),'DT5_below_minus_accuracy_count':int((d.DT5<-.4).sum()),'DT5_above_plus_accuracy_count':int((d.DT5>.4).sum()),
 'DT5_constant_offset_interval_for_all_T0_lower_C':float(d.DT5.max()-.4),'DT5_constant_offset_interval_for_all_T0_upper_C':float(d.DT5.min()+.4),'DT5_all_T0_constant_offset_interval_nonempty':bool(d.DT5.max()-.4<=d.DT5.min()+.4),
 'DT4_minus_DT5_min_C':float(grad.min()),'DT4_minus_DT5_max_C':float(grad.max()),'DT4_minus_DT5_mean_C':float(grad.mean()),
 'gradient_nominal_spacing_m':.125,'difference_bounded_accuracy_C':.8,'gradient_bounded_accuracy_K_m_ignoring_depth_uncertainty':6.4,
 'finite_difference_temporal_derivative_worstcase_error_C_hour':1.6,'interpretation':'Separately bounded per-reading errors; no independence assumption; sensor accuracy is not a probability distribution. Stable offsets would cancel in temporal differences but stability is not independently calibrated.'}
(Path(__file__).resolve().parent/'PIRAMIDE_COLD_SIGNAL_AUDIT.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
