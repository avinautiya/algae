"""Plot audit outputs without fitting or changing archived author products."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parent
rows=json.loads((ROOT/'HEAD_DROP_NOMINAL_AUDIT.json').read_text())['cases']
x=np.arange(len(rows))
fig,axes=plt.subplots(1,2,figsize=(11,4),constrained_layout=True)
axes[0].plot(x,[r['source_expression_drop_m'] for r in rows],'o-',label='Archived source expression')
axes[0].plot(x,[r['piezometric_drop_m'] for r in rows],'s-',label='Difference of piezometric heads')
axes[0].set_ylabel('Nominal water-head drop (m)')
axes[1].plot(x,[r['source_expression_nominal_f'] for r in rows],'o-',label='Archived expression, nominal inputs')
axes[1].plot(x,[r['piezometric_nominal_f'] for r in rows],'s-',label='Piezometric expression, nominal inputs')
axes[1].scatter(x,[r['source_exported_MC_mean_f'] for r in rows],marker='x',color='black',label='Original exported Monte Carlo mean')
axes[1].set_ylabel('Width-fixed Darcy friction factor')
labels=['11 Jul 09h','11 Jul 12h','16 Jul','25 Jul','30 Jul 14h','30 Jul 18h','31 Jul']
for ax in axes:
    ax.set_xticks(x,labels,rotation=35,ha='right');ax.grid(alpha=.2);ax.legend(fontsize=7)
fig.suptitle('Author-notebook head-expression audit\nNominal sensitivity only; not a corrected field-calibration product',fontsize=12)
for suffix in ['png','pdf']:fig.savefig(ROOT/f'HEAD_DROP_AUDIT.{suffix}',dpi=180)
