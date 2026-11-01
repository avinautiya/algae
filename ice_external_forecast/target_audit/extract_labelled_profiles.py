"""Digitize labelled Fig5 profiles as relative elevations, not invented hours."""
import argparse,csv,hashlib,json,statistics,datetime as dt
from pathlib import Path
from PIL import Image,ImageDraw

X_TICKS=[(-60,81),(-45,137),(-30,194),(-15,250),(0,308),(15,365),(30,421),(45,479),(60,536),(75,593),(90,650)]
Y_TICKS=[(40,6),(20,79),(0,152),(-20,225),(-40,298)]
SECTIONS={'A':319,'B':373,'C':433}
DATES=['2005-08-25','2005-08-26','2005-08-29']
SOURCE_SHA256='d5b1d127aa81cd10b14c543b27fa4cff16f51f94521541ce463433938b0ec389'

def fit(anchors):
    a=statistics.mean(x for x,y in anchors);b=statistics.mean(y for x,y in anchors)
    scale=sum((x-a)*(y-b) for x,y in anchors)/sum((x-a)**2 for x,y in anchors)
    return scale,b-scale*a

def color(date,rgb):
    r,g,b=rgb
    if date=='2005-08-25':return r>130 and r>1.5*g and r>1.5*b
    if date=='2005-08-26':return r>50 and b>70 and r>1.2*g and b>1.4*g
    return g>85 and b>85 and g>1.12*r and b>1.15*r

def main():
    a=argparse.ArgumentParser();a.add_argument('image',type=Path);a.add_argument('--out',type=Path,required=True);args=a.parse_args()
    if hashlib.sha256(args.image.read_bytes()).hexdigest()!=SOURCE_SHA256:raise ValueError('source figure hash mismatch')
    im=Image.open(args.image).convert('RGB')
    if im.size!=(709,367):raise ValueError('expected unscaled native embedded Fig5 raster')
    sx,ox=fit(X_TICKS);sy,oy=fit(Y_TICKS);rows=[];proof=im.copy();draw=ImageDraw.Draw(proof)
    for section,x in SECTIONS.items():
        for date in DATES:
            ys=[y for y in range(90,175) if color(date,im.getpixel((x,y)))]
            if not ys or max(ys)-min(ys)+1!=len(ys):raise ValueError(f'ambiguous/missing profile run: {section} {date}')
            cy=statistics.mean(ys);neighbors=[y for xx in range(x-2,x+3) for y in range(90,175) if color(date,im.getpixel((xx,y)))]
            bound=(max(abs(y-cy) for y in neighbors)+2)/abs(sy)
            rows.append(dict(section=section,date=date,time_of_day='',timezone='not established',
                source_pixel_x=x,source_pixel_y=cy,profile_color_run_ymin=min(ys),profile_color_run_ymax=max(ys),
                relative_elevation_to_central_00_tick_m=(cy-oy)/sy,
                graphical_elevation_halfwidth_m=bound,distance_from_plotted_breach_origin_m=(x-ox)/sx,
                graphical_distance_halfwidth_m=2/abs(sx),absolute_altitude_offset='not assigned; leading altitude digits clipped in published figure',
                target_kind='drawn measured longitudinal profile sampled at explicitly labelled section',
                not_raw_survey=True,measurement_error_separate='not separately established for individual channel-bed target'))
            draw.ellipse((x-3,cy-3,x+3,cy+3),outline='black',width=1)
    differences=[]
    for section in SECTIONS:
        selected={r['date']:r for r in rows if r['section']==section}
        for start,end in [('2005-08-25','2005-08-26'),('2005-08-26','2005-08-29'),('2005-08-25','2005-08-29')]:
            r0,r1=selected[start],selected[end];days=(dt.date.fromisoformat(end)-dt.date.fromisoformat(start)).days
            differences.append(dict(section=section,start_date=start,end_date=end,
                vertical_bed_lowering_m=r0['relative_elevation_to_central_00_tick_m']-r1['relative_elevation_to_central_00_tick_m'],
                endpoint_sum_graphical_bound_m=r0['graphical_elevation_halfwidth_m']+r1['graphical_elevation_halfwidth_m'],
                elapsed_time_lower_days=max(0,days-1),elapsed_time_upper_days=days+1,
                elapsed_time_kind='calendar-date-only envelope, not observed exact survey times',
                not_direct_thermal_volume=True,dependence='overlapping intervals and shared endpoints'))
    args.out.mkdir(parents=True,exist_ok=True)
    for name,data in [('FIG5_LABELLED_RELATIVE_ELEVATIONS.csv',rows),('FIG5_DATE_ONLY_BED_CHANGES.csv',differences)]:
        with (args.out/name).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
    proof.save(args.out/'FIG5_SELECTION_PROOF.png')
    manifest=dict(source_image_sha256=hashlib.sha256(args.image.read_bytes()).hexdigest(),source_figure=5,
        x_tick_calibration=X_TICKS,y_relative_tick_calibration=Y_TICKS,labelled_section_columns=SECTIONS,
        dates=DATES,points=len(rows),relative_differences=len(differences),
        no_absolute_altitude_inferred=True,no_exact_time_assigned=True,no_field_model_score=True,
        graphical_bounds='color-run/neighbor-column extent plus 2 vertical pixels; not confidence intervals or instrument error',
        parameter_freeze_reported_sha256='1b3734f3120eacaef7559130b89cde035d85b80c02594ed733c8e267188b347e')
    (args.out/'FIG5_EXTRACTION_MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(rows,indent=2))

if __name__=='__main__':main()
