"""Plots saved field predictions without fitting or selecting a winner."""
import os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/field-thermal-mpl')
import csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
r=json.loads((ROOT/'results/FIELD_TEMPERATURE_RESULTS.json').read_text());out=ROOT/'figures';out.mkdir(exist_ok=True)
models=['diffusion','quasisteady','lag_filter','ridge'];labels=['Homogeneous diffusion','Quasi-steady','Lag filter','Train-only ridge']
fig,ax=plt.subplots(figsize=(12,5),constrained_layout=True);x=np.arange(len(r['profile_results']))
for j,(m,l) in enumerate(zip(models,labels)):
 ax.bar(x+(j-1.5)*.19,[p['test_metrics'][m]['MAE_C'] for p in r['profile_results']],width=.19,label=l)
ax.set_xticks(x,[p['profile_id'].replace('KDT_','') for p in r['profile_results']],rotation=35,ha='right')
ax.set_ylabel('Held-out interior-temperature MAE (°C)');ax.set_title('Kennicott observed-boundary thermal test — every eligible profile\n13-day calibration objective; next 14 days withheld; one glacier, dependent histories')
ax.legend(ncol=2,fontsize=9);ax.grid(axis='y',alpha=.2)
for ext in ['png','pdf']:fig.savefig(out/f'field_temperature_all_profiles.{ext}',dpi=180)
plt.close(fig)
# Deterministic illustrative series: first sorted profile, deepest interior probe,
# first three held-out days. The complete campaign is displayed above.
p=sorted(r['profile_results'],key=lambda p:p['profile_id'])[0]
with (ROOT/'results'/p['prediction_csv']).open() as f:rr=list(csv.DictReader(f))
z=max(float(q['depth_m']) for q in rr)
rr=[q for q in rr if q['role']=='test' and float(q['depth_m'])==z]
rr=rr[:73];t=np.arange(len(rr))/24
fig,ax=plt.subplots(figsize=(11,4),constrained_layout=True)
ax.plot(t,[float(q['observed_C']) for q in rr],color='black',marker='.',linewidth=1.7,label='Measured interior')
for m,l in zip(models,labels):ax.plot(t,[float(q[m+'_C']) for q in rr],label=l,linewidth=1.1)
ax.set(xlabel='Elapsed days from held-out start',ylabel='Temperature (°C)',title=f"Deterministic illustration: {p['profile_id']}, depth {z*100:g} cm\nMeasured upper/lower probes supplied to every predictor")
ax.legend(ncol=3,fontsize=8);ax.grid(alpha=.2)
for ext in ['png','pdf']:fig.savefig(out/f'field_temperature_example.{ext}',dpi=180)
plt.close(fig)
print('Saved two figure pairs')
