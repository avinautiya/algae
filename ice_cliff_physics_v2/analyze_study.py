"""Recompute all frozen contrasts and publication figures; no parameter fitting."""
import json,math,csv
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def main():
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 out=ROOT/'results';fig=ROOT/'figures';fig.mkdir(exist_ok=True)
 records=[json.loads(p.read_text()) for p in sorted((out/'coupled').glob('*.json')) if p.name!='RUN_INDEX.json']
 done={d['id']:d for d in records if d['status']=='completed'}
 comparisons=[]
 for cap in [20000,100000]:
  for q in [5,20]:
   stem=f'C{cap}_Q{q}_';ref=done.get(stem+'fv32dt60')
   if ref is None:continue
   rt=ref['result']['totals']
   for name in ['rc','fv16','fv32','fv32face64','no_exchange_control']:
    if stem+name not in done:continue
    r=done[stem+name]['result']['totals']
    comparisons.append(dict(capacity_J_m2_K=cap,source_W_m2=q,comparator=name,
     covered_melt_difference_fraction=r['covered_face_melt_m3']/rt['covered_face_melt_m3']-1,
     total_face_melt_difference_fraction=r['face_melt_m3']/rt['face_melt_m3']-1,
     interpretation='controlled model comparison; reference is refined model, not environmental truth'))
 boundaries=[]
 for p in sorted((out/'boundary').glob('W*.json')):
  d=json.loads(p.read_text());a=d['analytic_robin'];fine=d['finite_volume'][-1]
  for t in ['0.1','1.0','10.0']:
   exact=a[t]['ice_heat_J']
   boundaries.append(dict(case=d['id'],time_over_diffusion_time=float(t),exact_ice_J=exact,
    rc_relative_error=d['rc'][t]['error_J']/exact,
    refined_fv_relative_error=fine['samples'][t]['ice_heat_J']/exact-1,
    refined_fv_cells=fine['cells'],refined_fv_steps_per_diffusion_time=fine['steps_per_diffusion_time']))
 geometry=[]
 for tw in ['0.5','1']:
  candidates=[d for k,d in done.items() if k.startswith('G') and float(d['configuration']['periods'][0]['water_superheat_K'])==float(tw)]
  for d in candidates:
   geometry.append(dict(id=d['id'],water_temperature_C=float(tw),totals=d['result']['totals'],final_undercut=d['result']['snapshots'][-1]['undercut']))
 summary=dict(freeze_commit='b9d6f39d09f7ead306ea87c852da8f14c53e9e89',case_count=len(records),completed=len(done),
  failures=[{k:d.get(k) for k in ['id','status','reason']} for d in records if d['status']!='completed'],
  maximum_relative_energy_residual=max(d['result']['max_time_series_relative_energy_residual'] for d in done.values()),
  comparisons=comparisons,analytic_boundary_errors=boundaries,geometry=geometry,
  limitations=['The refined coupled model is not an independent physical truth.',
  'Equal wet/dry net source is a controlled boundary, not identical environmental weather.',
  'Geometric undercut is not fracture or collapse.', 'Field temperature and ablation tests are separate component tests.'])
 (out/'ANALYSIS.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
 with (out/'COMPARISONS.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(comparisons[0]));w.writeheader();w.writerows(comparisons)
 plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
 fig1,ax=plt.subplots(2,2,figsize=(11,7),layout='constrained')
 for j,t in enumerate([.1,1.]):
  rows=[r for r in boundaries if r['time_over_diffusion_time']==t];x=list(range(len(rows)))
  for col,key,label in [(0,'rc_relative_error','Single-temperature debris'),(1,'refined_fv_relative_error','Resolved diffusion: 64 cells, 400 steps / diffusion time')]:
   a=ax[j,col];yy=[100*r[key] for r in rows];a.plot(x,yy,'o-')
   a.axhline(0,color='k',lw=.7);a.set_xticks(x,[r['case'] for r in rows],rotation=60)
   a.set_title(f'{label}\nt / diffusion time = {t:g}',fontsize=10);a.set_ylabel('Ice-heat error vs analytic (%)')
   a.text(.03,.96,f'Range: {min(yy):.3g} to {max(yy):.3g}%',transform=a.transAxes,va='top',bbox=dict(facecolor='white',alpha=.8,edgecolor='none'))
 fig1.suptitle('Independent analytic verification — early-time error remains timestep-limited')
 fig1.savefig(fig/'01_analytic_heat_errors.png',dpi=180);fig1.savefig(fig/'01_analytic_heat_errors.pdf');plt.close(fig1)
 fig2,ax=plt.subplots(figsize=(8,4),layout='constrained')
 for name,marker in [('rc','o'),('fv16','s'),('fv32','^'),('no_exchange_control','x')]:
  rows=[r for r in comparisons if r['comparator']==name]
  ax.plot(range(len(rows)),[100*r['covered_melt_difference_fraction'] for r in rows],marker+'-',label={'rc':'Single-temperature debris','fv16':'Resolved diffusion: 16 cells','fv32':'Resolved diffusion: 32 cells','no_exchange_control':'Water–debris exchange disabled (intervention)'}[name])
 ax.set_xticks(range(4),['C20k / q5','C20k / q20','C100k / q5','C100k / q20']);ax.set_xlabel('Areal heat capacity (J m$^{-2}$ K$^{-1}$) / imposed net source (W m$^{-2}$)');ax.set_ylabel('Covered-footprint melt difference (%)');ax.set_title('Conditional contrasts against FV32 / 60-s model');ax.axhline(0,color='k',lw=.7);ax.legend(fontsize=8)
 fig2.savefig(fig/'02_matched_covered_melt.png',dpi=180);fig2.savefig(fig/'02_matched_covered_melt.pdf');plt.close(fig2)
 fig3,ax=plt.subplots(figsize=(6,5),layout='constrained')
 for name in ['C20000_Q5_fv32dt60','C20000_Q20_fv32dt60']:
  p=done[name]['result']['snapshots'][-1]['face_profile'];bands=p['bands'];xx=[];zz=[]
  for b in bands:xx.extend([1000*b['front_x_m']]*2);zz.extend([b['z_low_m'],b['z_high_m']])
  ax.plot(xx,zz,label=name.replace('C20000_','').replace('_fv32dt60',''))
 ax.set_ylim(-.09,.15);ax.set_xlabel('Horizontal retreat (mm)');ax.set_ylabel('Elevation (m)');ax.set_title('Height-resolved near-toe geometry\nHorizontal ledges are adiabatic; no collapse prediction');ax.legend()
 fig3.savefig(fig/'03_toe_profiles.png',dpi=180);fig3.savefig(fig/'03_toe_profiles.pdf');plt.close(fig3)
 print(json.dumps({k:summary[k] for k in ['case_count','completed','failures','maximum_relative_energy_residual']},indent=2))
if __name__=='__main__':main()
