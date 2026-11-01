"""Digitize Fig6 forcing panels only; never inspect/extract panel-b target pixels."""
from pathlib import Path
from datetime import datetime,timedelta
import csv,json,hashlib,os
import numpy as np
from PIL import Image,ImageDraw
ROOT=Path(__file__).resolve().parent
SOURCE_ROOT=Path(os.environ.get('ROCHE_SOURCE_DIR',str(ROOT)))
SOURCE=SOURCE_ROOT/'figure6_source.jpeg'
X0,X120=59.,782.
T_Y12,T_Y06=53.,191.
Q_Y0,Q_Y12=655.,449.

def runs(indices):
    if not len(indices):raise ValueError('no plotted forcing pixels')
    breaks=np.flatnonzero(np.diff(indices)>1)+1
    groups=np.split(indices,breaks)
    groups.sort(key=len,reverse=True)
    if len(groups)>1 and len(groups[0])==len(groups[1]):raise ValueError('ambiguous equal-length traces')
    return groups[0]

def main():
    source=Image.open(SOURCE).convert('RGB');gray=np.array(source.convert('L'))
    rows=[];start=datetime(2005,8,25,19)
    for hour in range(2,112,2):
        x=round(X0+(X120-X0)*hour/120)
        # Deliberately disjoint forcing ROIs: panel b is never sampled.
        ty=runs(np.flatnonzero(gray[57:188,x]<135)+57)
        qy=runs(np.flatnonzero(gray[430:640,x]<135)+430)
        yT=float(np.median(ty));yQ=float(np.median(qy))
        tw=1.2-(yT-T_Y12)*.6/(T_Y06-T_Y12)
        q=(Q_Y0-yQ)*1.2/(Q_Y0-Q_Y12)
        rows.append(dict(elapsed_hours=hour,source_clock_naive=(start+timedelta(hours=hour)).isoformat(),
                         temperature_graph_C=tw,discharge_graph_m3_s=q,x_source_px=x,
                         temperature_y_source_px=yT,discharge_y_source_px=yQ,
                         temperature_curve_run_px=len(ty),discharge_curve_run_px=len(qy),
                         preferred_regular_breach_window=12<=hour<=108))
    with (ROOT/'FORCING_GRAPH_SAMPLES.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    # Proof exposes only temperature/discharge input panels, never bed-change data.
    top=source.crop((0,0,799,214));bottom=source.crop((0,414,799,708))
    proof=Image.new('RGB',(799,508),'white');proof.paste(top,(0,0));proof.paste(bottom,(0,214));draw=ImageDraw.Draw(proof)
    for row in rows:
        x=row['x_source_px'];yt=row['temperature_y_source_px'];yq=row['discharge_y_source_px']-414+214
        for y in [yt,yq]:draw.ellipse((x-2,y-2,x+2,y+2),outline=(210,0,0),width=1)
    proof.save(ROOT/'FORCING_SELECTION_PROOF.png')
    manifest=dict(source_pdf_url='https://www.cambridge.org/core/services/aop-cambridge-core/content/view/C183D075F0E9ABCD7CE5251D6096FCA4/S0022143000206357a.pdf/outburst-flood-hazard-for-glacier-dammed-lac-de-rochemelon-france.pdf',
                  pdf_sha256=hashlib.sha256((SOURCE_ROOT/'vincent2010_publisher.pdf').read_bytes()).hexdigest(),
                  original_jpeg_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),pdf_page_zero_based=3,pdf_image_xref=91,
                  graph_size_px=[799,708],x_time_anchors=[[59,0],[782,120]],temperature_anchors=[[53,1.2],[191,.6]],
                  discharge_anchors=[[655,0],[449,1.2]],foreground_threshold_gray=135,
                  point_method='At prespecified even-hour x columns, choose longest connected dark vertical run within panel a/c ROIs; median y. These are samples of published plotted lines, not original logger timestamps.',
                  sampled_regions_y=[[57,188],[430,640]],withheld_target_region_y=[214,414],
                  n_graph_samples=len(rows),temperature_pixel_error_C=2*.6/(191-53),discharge_pixel_error_m3_s=2*1.2/(655-449),time_pixel_error_h=2*120/(782-59),
                  error_status='Conservative graph-coordinate engineering bounds from +/-2 pixels, additional unknown interpolation/line-rendering/axis-anchor errors. Not statistical confidence intervals.',
                  source_temperature_accuracy_C=.1,source_discharge_accuracy_m3_s=.1,timezone='not established; retain source-naive clock',
                  preferred_window_hours=[12,108],raw_numerical_source_available=False,bed_target_values_extracted=False)
    (ROOT/'FORCING_PROVENANCE.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('Extracted forcing only',len(rows),'graph samples, preferred12–108h')
if __name__=='__main__':main()
