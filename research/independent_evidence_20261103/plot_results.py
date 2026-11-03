"""Reproducible deterministic views; never fits or selects model parameters."""
import json,os,tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parent
CACHE=Path(tempfile.gettempdir())/'glacier-plot-cache-20261103';CACHE.mkdir(exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR',str(CACHE))
os.environ.setdefault('XDG_CACHE_HOME',str(CACHE))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
def main():
 out=ROOT/'figures';out.mkdir(exist_ok=True)
 folder=ROOT/'extension/results/frozen_point_diagnostic'
 data=json.loads((folder/'RESULTS.json').read_text());rows=json.loads((folder/'PREDICTIONS.json').read_text())
 names=['diffusion','instantaneous_gradient','training_median','causal_boundary_ridge'];labels=['Diffusion','Instantaneous\ngradient','Training\nmedian','Causal ridge']
 plt.rcParams.update({'axes.spines.top':False,'axes.spines.right':False,'font.size':10})
 f,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
 for ax,r in zip(axes,data['results']):
  vals=[r['scores'][n]['mae_W_m2'] for n in names];ax.bar(range(4),vals,color=['#386cb0','#7fc97f','#beaed4','#fdc086']);ax.set_xticks(range(4),labels);ax.set_yscale('log');ax.set_ylim(.4,20);ax.set_ylabel('Held-out MAE (raw W/m²; log scale)');ax.set_title('Site '+str(r['site'])+' — all336 records')
  for i,v in enumerate(vals):ax.text(i,v*1.08,f'{v:.3f}',ha='center',fontsize=9)
 f.suptitle('Conditional sensor-signal test — declared superiority rule FAILED')
 for ext in ['png','pdf']:f.savefig(out/('plate_errors.'+ext),dpi=180)
 plt.close(f)
 f,axes=plt.subplots(2,1,figsize=(10,6),layout='constrained')
 for ax,site in zip(axes,[1,2]):
  selected=[r for r in rows if r['site']==site][:96];hours=range(len(selected))
  ax.plot(hours,[r['observed_raw_plate_W_m2'] for r in selected],color='black',label='Observed raw plate')
  for key,col in [('diffusion','#386cb0'),('causal_boundary_ridge','#e17c05')]:ax.plot(hours,[r['predictions'][key] for r in selected],color=col,label=key)
  ax.set_title(f'Site{site}: first96 evaluation hours (fixed chronological selection)');ax.set_ylabel('Raw W/m²');ax.legend(fontsize=8,ncol=3)
 axes[-1].set_xlabel('Hours from first evaluation record')
 f.suptitle('Representative-point diagnostic; not native interval-flux or melt validation')
 for ext in ['png','pdf']:f.savefig(out/('plate_dynamics.'+ext),dpi=180)
 plt.close(f)
if __name__=='__main__':main()
