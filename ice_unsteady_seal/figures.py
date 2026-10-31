"""Plot the recorded verification checks; creates no new physical experiment."""

import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record',required=True); parser.add_argument('--output',required=True)
    args=parser.parse_args(); record=Path(args.record); folder=Path(args.output)
    if folder.exists(): raise FileExistsError(folder)
    data=json.loads(record.read_text()); folder.mkdir(parents=True)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                         'svg.fonttype':'none','savefig.facecolor':'white'})
    fig,axes=plt.subplots(2,2,figsize=(11,7.6),constrained_layout=True)
    forecast=data['forecast']; base=forecast['base']; x=np.array(base['x']); v=np.array(base['velocity'])
    bed=np.array(base['bed']); q=base['configuration']['discharge']; beta=base['configuration']['beta']
    axes[0,0].plot(x,bed,label='Ice bed',color='#234e70')
    axes[0,0].plot(x,bed+(q/v)**beta,label='Water surface',color='#159b99')
    axes[0,0].set(xlabel='Inner coordinate X',ylabel='Elevation',title='A   Nonuniform seal base')
    axes[0,0].legend(frameon=False)
    rows=data['uniform_shape_quadrature']; spacing=np.array([r['spacing'] for r in rows])
    errors=np.array([r['absolute_error'] for r in rows])
    axes[0,1].loglog(spacing,errors,'o-',color='#a24437',label='Independent quadrature error')
    axes[0,1].loglog(spacing,errors[-1]*(spacing/spacing[-1])**2,'--',color='.5',label='Second-order guide')
    axes[0,1].set(xlabel='Station spacing',ylabel='Discharge-derivative error',title='B   Hydraulic spatial verification')
    axes[0,1].legend(frameon=False,fontsize=8)
    rows=data['canonical_hydraulic_matching']; cutoff=np.array([r['left_cutoff'] for r in rows])
    tail=np.array([r['reference_difference'] for r in rows]); leading=np.array([r['leading_asymptotic_term'] for r in rows])
    axes[1,0].loglog(cutoff,tail,'o-',color='#234e70',label='Tail-corrected finite reference')
    axes[1,0].loglog(cutoff,leading,'--',color='#a24437',label=r'$\frac{1}{2} L^{-4}$ asymptote')
    axes[1,0].set(xlabel='Upstream cutoff L',ylabel='Hydraulic coefficient difference',title='C   Pond-limit verification')
    axes[1,0].legend(frameon=False,fontsize=8)
    samples=forecast['samples']; times=np.array([r['time'] for r in samples])
    axes[1,1].plot(times,[r['discharge_perturbation'] for r in samples],color='#159b99',label='Discharge perturbation')
    axes[1,1].plot(times,[r['lake_head_perturbation'] for r in samples],color='#a24437',label='Total-head perturbation')
    axes[1,1].set(xlabel='Inner time T',ylabel='Dimensionless perturbation',title='D   Manufactured linear scenario')
    axes[1,1].legend(frameon=False,fontsize=8)
    fig.suptitle('Unsteady ice-seal model verification\nDimensionless checks — no environmental skill claim',fontsize=13)
    for ext in ['png','svg']: fig.savefig(folder/('verification_panels.'+ext),dpi=180)
    plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4),constrained_layout=True)
    for key,label,color in [('total_head_ledger_residual','Surrogate total-head bookkeeping','#234e70'),
                            ('surface_stage_storage_residual','Finite-cutoff surface-stage storage residual','#a24437')]:
        values=np.abs([r[key] for r in samples]); ax.semilogy(times,np.maximum(values,1e-22),label=label,color=color)
    ax.set(xlabel='Inner time T',ylabel='Absolute residual (dimensionless)',
           title='Bookkeeping and finite-cutoff stage are different quantities')
    ax.legend(frameon=False,fontsize=8)
    for ext in ['png','svg']: fig.savefig(folder/('head_and_stage_ledger.'+ext),dpi=180)
    plt.close(fig)
    if 'matched_outer_gauge' in data:
        matched=data['matched_outer_gauge']; fig,axes=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
        axes[0].plot(matched['kernel_times'],matched['renewal_kernel'],color='#234e70')
        axes[0].set(xlabel='Inner time T',ylabel='Renewal-kernel value',title='A   Derived outer-gauge shape response')
        rows=matched['forecast']['samples']; times=np.array([r['time'] for r in rows])
        axes[1].plot(times,[r['outer_seal_relative_head'] for r in rows],label='Head above outer seal',color='#a24437')
        axes[1].plot(times,[r['inner_crest_relative_head'] for r in rows],label='Head above inner crest',color='#159b99')
        axes[1].set(xlabel='Inner time T',ylabel='Dimensionless head perturbation',title='B   Two control-height definitions')
        axes[1].legend(frameon=False,fontsize=8)
        fig.suptitle('Canonical linear closure verification — finite discretization',fontsize=12)
        for ext in ['png','svg']: fig.savefig(folder/('matched_shape_response.'+ext),dpi=180)
        plt.close(fig)
    hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(folder.iterdir())}
    manifest={'record':str(record),'record_sha256':hashlib.sha256(record.read_bytes()).hexdigest(),
              'plotter_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'files':hashes,
              'caption':'Equation and numerical verification of a dimensionless truncated model. Panel C uses a tail-corrected finite reference, not a certified exact error. Panel D is a manufactured scenario. These figures establish neither a new glacier discovery nor forecast skill.'}
    (folder/'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({'figures':str(folder.resolve()),'count':len(hashes)}))


if __name__=='__main__': main()
