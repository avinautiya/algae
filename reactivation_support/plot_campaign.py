"""Plot completed run JSON without inventing absent observables."""
import argparse,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def plot(inputs,out,*,label):
    fig,ax=plt.subplots(2,2,figsize=(11,7),layout='constrained')
    for f in inputs:
        run=json.loads(Path(f).read_text());records=run['records']
        t=np.array([r['day'] for r in records]);name=Path(f).stem
        if all('f_toe' in r for r in records):ax[0,0].plot(t,[r['f_toe'] for r in records],label=name)
        if all('cum_melt_toe' in r for r in records):ax[0,1].plot(t,[r['cum_melt_toe'] for r in records],label=name)
        if all('crest_toe' in r and 'crest_byp' in r for r in records):ax[1,0].plot(t,[r['crest_toe']-r['crest_byp'] for r in records],label=name)
        if all('energy_resid_rel' in r for r in records):ax[1,1].plot(t,[r['energy_resid_rel'] for r in records],label=name)
    titles=['Toe discharge / prescribed inflow','Cumulative local toe melt','Reported toe crest − bypass crest','Hydraulic-window thermal ledger residual']
    labels=['Fraction','Ice volume (m³)','Elevation difference (m)','Dimensionless residual']
    for a,title,y in zip(ax.ravel(),titles,labels):
        a.set(title=title,xlabel='Physical day',ylabel=y)
        if a.lines:a.legend(fontsize=7)
        else:a.text(.5,.5,'Required observable not recorded',ha='center',va='center',transform=a.transAxes)
    ax[1,1].set_yscale('symlog',linthresh=1e-12)
    fig.suptitle(label,fontsize=12)
    out=Path(out);out.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(out.with_suffix('.png'),dpi=180);fig.savefig(out.with_suffix('.pdf'));plt.close(fig)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('inputs',nargs='+',type=Path)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--label',required=True,help='Explicit evidence class, e.g. frozen evaluation or development only')
    a=p.parse_args();plot(a.inputs,a.output,label=a.label)
