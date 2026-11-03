"""Scientific development figures; all data already inspected/released."""
import json,csv,datetime as dt
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).parent;r=json.loads((P/'results_v1/DEVELOPMENT_RESULTS.json').read_text())
fig,axs=plt.subplots(1,2,figsize=(10,4),sharey=True)
for ax,tower in zip(axs,['Tower1','Tower2']):
 rows=[x for x in r['all_trials'] if x['tower']==tower and x['supported']]
 ax.plot([x['timescale_reference_days']for x in rows],[x['DEVELOPMENT_MAE_mmwe']for x in rows],'o-',label='Radiation reservoir')
 for method,color,label in [('air_pdd','tab:orange','Locked air PDD'),('calendar_exp_air','tab:green','Locked calendar')]:
  ax.axhline(r['locked_old_baselines'][tower][method]['DEVELOPMENT_MAE_mmwe'],color=color,label=label,ls='--')
 ax.set_title(tower+' — DEVELOPMENT');ax.set_xlabel('Reference reservoir timescale (days)');ax.grid(alpha=.2)
axs[0].set_ylabel('Mean absolute period ablation error (mm w.e.)');axs[1].legend(fontsize=8)
fig.suptitle('Finite radiation/storage family fails strong development baselines',fontsize=12);fig.tight_layout()
fig.savefig(P/'development_family_errors.png',dpi=180);fig.savefig(P/'development_family_errors.svg');plt.close(fig)
with (P/'results_v1/Tower1_tau0.25_trace.csv').open()as f:rows=list(csv.DictReader(f))
rows=[x for x in rows if '2022-12-17'<=x['SOURCE_LOCAL_LABEL']<'2022-12-21']
t=[dt.datetime.fromisoformat(x['SOURCE_LOCAL_LABEL'])for x in rows]
fig,axs=plt.subplots(3,1,figsize=(10,6),sharex=True)
axs[0].plot(t,[float(x['net_radiation_proxy_W_m2'])for x in rows],color='tab:orange');axs[0].axhline(0,c='k',lw=.5);axs[0].set_ylabel('Q proxy (W m⁻²)')
axs[1].plot(t,[float(x['effective_E_start_J_m2'])/1e6 for x in rows],color='tab:blue');axs[1].axhline(0,c='k',lw=.5);axs[1].set_ylabel('Effective E (MJ m⁻²)')
axs[2].plot(t,[float(x['delivered_proxy_J_m2'])/1800 for x in rows],color='tab:green');axs[2].set_ylabel('Proxy delivery (W m⁻²)');axs[2].set_xlabel('Source local time (UTC offset not established)')
for ax in axs:ax.grid(alpha=.2)
fig.suptitle('Illustrative uncalibrated cooling/storage dynamics — not measured ice state',fontsize=11);fig.autofmt_xdate();fig.tight_layout();fig.savefig(P/'effective_reservoir_dynamics.png',dpi=180);fig.savefig(P/'effective_reservoir_dynamics.svg');plt.close(fig)
print('Two development figures rendered.')
