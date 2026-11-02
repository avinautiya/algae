import os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/layered-mpl')
import json
from pathlib import Path
import numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parent;r=json.loads((P/'results/RESULTS.json').read_text());p=r['profiles'];models=['layered','homogeneous','quasisteady','lag_filter','ridge'];labels=['Two-layer diffusion','Homogeneous diffusion','Quasi-steady','Lag filter','Train-only ridge'];fig,ax=plt.subplots(figsize=(14,5));x=np.arange(len(p));w=.16
for j,(m,label) in enumerate(zip(models,labels)):ax.bar(x+(j-2)*w,[q['metrics'][m]['MAE_C'] for q in p],w,label=label)
ax.set_xticks(x);ax.set_xticklabels([q['profile_id'].replace('KDT_','') for q in p]);ax.set_ylabel('Reserved interior-temperature MAE (°C)');ax.set_title('Fixed-midpoint heterogeneity on reserved days28–42\nEarly13-day training only; every eligible profile; observed endpoint forcing');ax.legend(ncol=3,fontsize=9);ax.grid(axis='y',alpha=.2);fig.tight_layout();(P/'figures').mkdir(exist_ok=True)
for ext in ['png','pdf']:fig.savefig(P/'figures'/f'reserved_layered_all_profiles.{ext}',dpi=160)
