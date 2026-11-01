"""Regenerate the saved coarse operator diagnostic; no solver call."""
import json
from pathlib import Path
import numpy as np

def main():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    folder=Path(__file__).parent
    fig,axs=plt.subplots(1,2,figsize=(10,3.8),layout='constrained')
    for mode,label in [('legacy','Pinned bed update'),('ice_only','Ice-only unbounded update')]:
        data=json.loads((folder/('paired_'+mode+'.json')).read_text())
        rows=data['run']['records'];updates=data['operator_updates']
        axs[0].plot([r['day'] for r in rows],[r['f_toe'] for r in rows],label=label)
        axs[1].plot([u['time_day']+.25 for u in updates],np.cumsum([u['declared_melt_m3']-u['actual_ice_lowering_m3'] for u in updates]),label=label)
    axs[0].set(xlabel='Time (days)',ylabel='Toe-flow fraction',title='Early contact-loss sensitivity')
    axs[1].set(xlabel='Time (days)',ylabel='Counted minus geometric ice melt (m3)',title='Whole-domain phase bookkeeping')
    axs[0].legend(fontsize=8)
    fig.suptitle('Exploratory coarse dx = 0.1 m, 3.5 days; not production data',fontsize=11)
    for ext in ('png','pdf'):fig.savefig(folder/('paired_operator_diagnostic.'+ext),dpi=180)

if __name__=='__main__':main()
