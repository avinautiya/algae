"""Plot the frozen retrospective comparison; no fitting or selection here."""
import csv
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parent
rows=list(csv.DictReader((ROOT/'CASE_PREDICTIONS.csv').open()))
results=json.loads((ROOT/'RESULTS.json').read_text())
labels={'piezometric_fit':'Piezometric head, fitted f','archived_plus_fit':'Archived plus-head, fitted f',
        'training_constant_velocity':'Training constant velocity','piezometric_fixed_f03':'Piezometric head, fixed f=0.3',
        'archived_plus_fixed_f03':'Archived plus-head, fixed f=0.3'}
fig,axes=plt.subplots(1,2,figsize=(12,4.6),constrained_layout=True)
x=np.arange(4)
for key,label in labels.items():
    case=[r for r in rows if float(r['width_m'])==2 and r['split']=='test' and r['model']==key]
    axes[0].plot(x,[float(r['predicted_velocity_m_s']) for r in case],'o-',label=label,alpha=.85)
axes[0].scatter(x,[float(r['observed_velocity_m_s']) for r in case],marker='x',s=75,color='black',label='Measured tracer mean',zorder=5)
axes[0].set_xticks(x,['25 Jul','30 Jul 14h','30 Jul 18h','31 Jul'],rotation=25)
axes[0].set(ylabel='Velocity (m/s)',title='Temporal test cases, width 2 m')
axes[0].legend(fontsize=7)
widths=[1.5,2.,2.5]
for key,label in labels.items():
    axes[1].plot(widths,[results['widths'][str(w)]['test_metrics'][key]['MAE_m_s'] for w in widths],'o-',label=label)
axes[1].set(xlabel='Prespecified channel width (m)',ylabel='Test MAE (m/s)',title='Every width sensitivity retained')
axes[1].set_xticks(widths)
for ax in axes:ax.grid(alpha=.2)
fig.suptitle('Correcting the head expression does not establish transferable constant-friction skill\nRetrospective conditional test: 3 training and 4 dependent test experiments',fontsize=12)
for suffix in ['png','pdf']:fig.savefig(ROOT/f'VELOCITY_TRANSFER.{suffix}',dpi=180)
