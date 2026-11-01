"""Digitize saved corrected daily-temperature curves; no model fitting.

Usage: python extract.py --image SOURCE.png --out OUTPUT_DIRECTORY
Requires Pillow. Pixel coordinates refer to original 581 x 408 raster.
"""
import argparse, csv, datetime as dt, hashlib, json, statistics
from pathlib import Path
from PIL import Image, ImageDraw

EXPECTED = '3474fcd55ae569d9d3723feec37fab3a07e959fb07c7a2a69f74d694dac63e62'
X = [(0,86),(7,142),(14,197),(24,276)]
Y = [(0,326),(.1,267),(.2,208),(.3,148),(.4,89),(.5,30)]
AMBIGUOUS = {'P3': {17,18,20,23,26}, 'P5': {18,23}}

def fit(anchors):
    a = statistics.mean(x for x,y in anchors); b = statistics.mean(y for x,y in anchors)
    m = sum((x-a)*(y-b) for x,y in anchors)/sum((x-a)**2 for x,y in anchors)
    return m, b-m*a

def color_match(station, rgb):
    r,g,b = rgb
    return (g>75 and g>r*1.4 and g>b*1.25) if station=='P3' else (r>110 and r>g*3 and r>b*3)

def extract(image):
    if hashlib.sha256(image.read_bytes()).hexdigest()!=EXPECTED:
        raise ValueError('source image hash mismatch')
    im=Image.open(image).convert('RGB')
    if im.size!=(581,408): raise ValueError('wrong original raster dimensions')
    sx,ox=fit(X); sy,oy=fit(Y); rows=[]; proof=im.copy(); draw=ImageDraw.Draw(proof)
    for day in range(16,31):
        date=dt.datetime(2019,7,day,12); x=ox+sx*((date-dt.datetime(2019,7,8)).total_seconds()/86400)
        col=round(x)
        for station in ['P3','P5']:
            ys=[y for y in range(12,326) if color_match(station,im.getpixel((col,y)))]
            if not ys or max(ys)-min(ys)+1!=len(ys):
                raise ValueError(f'missing or multiple color run: {day}, {station}, {ys}')
            cy=statistics.mean(ys); temp=(cy-oy)/sy
            neighborhood=[y for xx in range(col-2,col+3) for y in range(12,326) if color_match(station,im.getpixel((xx,y)))]
            graphical=(max(abs(y-cy) for y in neighborhood)+2)/abs(sy)+.001
            identified=day not in AMBIGUOUS[station]
            rows.append(dict(station=station,date_time_noon_naive=date.isoformat(),timezone='not established',
                pixel_x=col,noon_subpixel_x=x,pixel_y=cy,color_run_ymin=min(ys),color_run_ymax=max(ys),
                temperature_C=temp,graphical_halfwidth_C=graphical,
                instrument_accuracy_bound_C=.05,unknown_bulk_minus_floor_offset_C='not identified',
                rounding_step_C=.001,visible_daily_vertex=identified,
                observation_kind='digitized corrected daily floor-temperature curve',
                near_collinear_curve_sample=not identified,raw_sensor_record=False,
                calibration_offsets_already_in_source=True))
            draw.ellipse((col-2,cy-2,col+2,cy+2),outline='black' if identified else 'purple')
    return rows,proof,dict(x_anchors=X,y_anchors=Y,x_scale=sx,x_intercept=ox,y_scale=sy,y_intercept=oy)

def main():
    p=argparse.ArgumentParser();p.add_argument('--image',type=Path,default=Path(__file__).resolve().parent.parent/'temperature_cell9_output0.png');p.add_argument('--out',type=Path,default=Path(__file__).resolve().parent);a=p.parse_args()
    rows,proof,axes=extract(a.image);a.out.mkdir(parents=True,exist_ok=True)
    for name,selected in [('daily_curve_samples_July16_30.csv',rows),('identifiable_vertices_July16_30.csv',[r for r in rows if r['visible_daily_vertex']])]:
        with (a.out/name).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(selected)
    proof.save(a.out/'selection_proof.png')
    (a.out/'manifest.json').write_text(json.dumps(dict(source_png_sha256=EXPECTED,
        source_repository='christopheogier/PlaineMorteLake_hydro_thermodynamics.jl',
        source_commit='9116f98dc7fef42be82dc4de50316f4a676ca8e9',source_notebook='Temperature_visualization.ipynb',
        source_cell_zero_based=9,source_output=0,axes=axes,curve_samples=len(rows),
        visible_vertices=sum(r['visible_daily_vertex'] for r in rows),
        time_convention='helper-functions.ipynb cell14: calendar-day average, labelled same-day noon; timezone unstated',
        sensor_source='https://tc.copernicus.org/articles/15/5133/2021/ section 3.1.2',
        bounds='graphical halfwidth and instrument accuracy are separate analyst/source bounds, not confidence intervals',
        no_bulk_representation_bound=True,model_fit=False),indent=2)+'\n')
    print(f'{len(rows)} curve samples; {sum(r["visible_daily_vertex"] for r in rows)} visible vertices')

if __name__=='__main__':main()
