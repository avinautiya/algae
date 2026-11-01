"""Own plots of graph-digitized factual forcing; no bed-target data."""
import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
rows=list(csv.DictReader((ROOT/'FORCING_GRAPH_SAMPLES.csv').open()))
t=[float(r['elapsed_hours']) for r in rows];temp=[float(r['temperature_graph_C']) for r in rows];q=[float(r['discharge_graph_m3_s']) for r in rows]
fig,axes=plt.subplots(2,1,figsize=(9,5.5),sharex=True,constrained_layout=True)
axes[0].plot(t,temp,'o-',ms=3,label='Graph-digitized lake temperature');axes[0].fill_between(t,[v-.1 for v in temp],[v+.1 for v in temp],alpha=.15,label='Source stated accuracy +/-0.1 C; not CI');axes[0].set_ylabel('Lake temperature (C)')
axes[1].plot(t,q,'o-',ms=3,label='Graph-digitized discharge from lake level');axes[1].fill_between(t,[max(0,v-.1) for v in q],[v+.1 for v in q],alpha=.15,label='Source stated accuracy +/-0.1 m3/s; not CI');axes[1].set_ylabel('Discharge (m3/s)');axes[1].set_xlabel('Elapsed hours from 25 August 2005 at 19:00, source clock')
for ax in axes:
 ax.axvline(12,color='gray',linestyle='--');ax.axvline(108,color='gray',linestyle='--');ax.grid(alpha=.2);ax.legend(fontsize=8)
fig.suptitle('Rochemelon external forcing only\nSamples of published curves; no measured bed-change targets included')
for suffix in ['png','pdf']:fig.savefig(ROOT/f'FORCING_INPUTS.{suffix}',dpi=180)
