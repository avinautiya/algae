from pathlib import Path
import json
from dataclasses import replace
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from model import Case,evolve
root=Path(__file__).parent;rec=root/'records';out=root/'figures';out.mkdir(exist_ok=True)
base=Case(discharge=.07,initial_debris_depth=.2,allow_supply=False,source_fractions=(.1,.1,.8))
h=[]
for label,pulse in [('no flood',None),('eightfold flood on day 1',dict(start_day=1.,multiplier=8.,duration_day=.25)),('eightfold flood on day 14',dict(start_day=14.,multiplier=8.,duration_day=.25))]:
 r=evolve([base],days=45,dt=112.5,pulse=pulse)
 h.append(dict(label=label,pulse=pulse,run=r))
(rec/'accepted_spatial_histories.json').write_text(json.dumps(h,indent=2))
fig,ax=plt.subplots(2,2,figsize=(11,8),layout='constrained')
for x in h:
 hist=x['run']['history'];t=[v['day'] for v in hist]
 ax[0,0].plot(t,[v['toe_flow'][0] for v in hist],label=x['label'])
 ax[0,1].plot(t,[v['toe_ice_z'][0]-v['bypass_z'][0] for v in hist],label=x['label'])
K=(2/3)**1.5*np.sqrt(9.81)*.5
ax[0,1].axhline((.07/K)**(2/3),ls='--',color='black',label='bare-route activation limit')
ax[0,0].set(xlabel='Day',ylabel='Cliff-side flow (m³/s)',title='Identical floods at different times')
ax[0,1].set(xlabel='Day',ylabel='Toe ice − bypass ice elevation (m)',title='The routing barrier outlasts the debris')
ax[0,0].legend(fontsize=8);ax[0,1].legend(fontsize=8)
r=json.loads((rec/'spatial_coverage_check.json').read_text())
for mode,label in [(False,'uniform-layer assumption'),(True,'grain-footprint coverage')]:
 inds=[i for i,c in enumerate(r['cases']) if c['footprint_shield']==mode]
 ax[1,0].plot([r['cases'][i]['source_fractions'][-1] for i in inds],[r['wet_days'][i] for i in inds],'o-',label=label)
ax[1,0].set(xlabel='Initial coarse-grain volume fraction',ylabel='Days with water contact (45-day run)',title='Sparse cover changes the verdict')
ax[1,0].legend(fontsize=8)
rr=json.loads((rec/'spatial_pulse_check.json').read_text())
sub=[x for x in rr if x['pulse']['multiplier']==8.]
ax[1,1].bar([str(int(x['pulse']['start_day'])) for x in sub],[x['toe_melt_m3'] for x in sub],color='steelblue')
ax[1,1].set(xlabel='Day of the identical six-hour flood',ylabel='45-day local toe melt (m³)',title='All debris exported in each case')
fig.suptitle('Finite-stock model with patchy grain coverage — exploratory, not field validation',fontsize=12)
fig.savefig(out/'spatial_coverage_and_reactivation.png',dpi=180);fig.savefig(out/'spatial_coverage_and_reactivation.pdf');plt.close(fig)
print([(x['label'],x['run']['wet_days'],x['run']['toe_melt_m3']) for x in h])
