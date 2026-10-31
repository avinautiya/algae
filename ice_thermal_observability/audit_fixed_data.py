"""Supplemental inverse-grid audit with one fixed observation vector per flow.

The main predeclared study regenerates observations from a uniform physical
profile on each grid. This additional audit holds the data vector fixed, so
discretization and changing observations are not conflated.
"""
import argparse
from pathlib import Path
import numpy as np
from .channel import Channel,CrossStream,plug_cell_bulk_kernel
from .forecast import dump
from .recovery import sharp_interval


def audit(output):
    path=Path(output)
    if path.exists():raise FileExistsError("Choose a new audit output")
    stations=[.02,.04];target=[.06,.12];rows=[];references={}
    for flow in ("plug","shear_free"):
        ch=Channel(velocity=flow)
        if flow=="plug":
            edges=np.array([0,ch.depth_m])
            data=np.array([float(plug_cell_bulk_kernel(ch,edges,x,1024)@np.array([.1])) for x in stations])
            difference=np.zeros(2)
            label="Analytical uniform-inlet sine series, 1024 terms; 2048-term agreement checked"
            more=np.array([float(plug_cell_bulk_kernel(ch,edges,x,2048)@np.array([.1])) for x in stations])
            if np.max(np.abs(data-more))>1e-14:raise ArithmeticError("Analytical truncation unresolved")
        else:
            coarse=CrossStream(ch,384);fine=CrossStream(ch,768)
            data=np.array([fine.bulk(np.full(768,.1),x) for x in stations])
            difference=data-np.array([coarse.bulk(np.full(384,.1),x) for x in stations])
            label="768-cell synthetic reference; 384-to-768 difference is an error indicator, not a certified bound"
        references[flow]=dict(bulk_data_K=data.tolist(),reference=label,
                              reference_refinement_difference_K=difference.tolist())
        for cells in (24,48,96,192):
            m=CrossStream(ch,cells);e=np.array([m.bulk_kernel(x) for x in stations])
            c=m.heat_kernel_W_K(*target)
            r=sharp_interval(c,e,data,np.zeros(2),np.zeros(cells),np.full(cells,.2))
            rows.append(dict(velocity=flow,cells=cells,fixed_bulk_data_K=data.tolist(),
                             lower_W=r["lower"],upper_W=r["upper"],minimax_W=r["minimax"],
                             max_dual_gap_W=max(abs(r["upper_certificate"]["dual_gap"]),abs(r["negative_lower_certificate"]["dual_gap"]))))
    result=dict(scope="Supplemental numerical-validation audit; no changes to primary protocol/test cases",
                references=references,rows=rows,continuum_error_certified=False)
    dump(path,result)
    print(path)
    return result


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--output",required=True);audit(p.parse_args().output)
