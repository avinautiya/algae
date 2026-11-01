"""Audit observed, corrected Himalayan profile point clouds without melt labels.

Cross-profile medians are descriptive terrain summaries, not extracted ice
boundaries, debris thickness, measured flow, or attributed melt volumes.
"""
import argparse
from collections import defaultdict
import csv
import datetime
import hashlib
import json
import math
from pathlib import Path
from statistics import median

FIELDS=['glacier','cliff_id','profile_id','survey_date','point_x','point_y',
        'point_z','distance_along','distance_to']


def load_points(path):
    groups=defaultdict(list)
    with path.open(newline='') as stream:
        reader=csv.DictReader(stream)
        if reader.fieldnames!=FIELDS:raise ValueError('unexpected published profile schema')
        for row in reader:
            datetime.date.fromisoformat(row['survey_date'])
            nums={k:float(row[k]) for k in FIELDS[4:]}
            if not all(math.isfinite(v) for v in nums.values()):
                raise ValueError('nonfinite profile coordinates')
            if abs(nums['distance_to'])>2.501:raise ValueError('outside documented extraction strip')
            key=(row['glacier'],int(row['cliff_id']),int(row['profile_id']),row['survey_date'])
            groups[key].append((nums['distance_along'],nums['point_z'],nums['distance_to']))
    return groups


def bin_profile(points,bin_width_m=.5,max_perpendicular_m=.5,min_points=5):
    if not math.isfinite(bin_width_m) or bin_width_m<=0:raise ValueError('positive bin width required')
    if not math.isfinite(max_perpendicular_m) or max_perpendicular_m<=0:raise ValueError('positive strip width required')
    if isinstance(min_points,bool) or not isinstance(min_points,int) or min_points<1:raise ValueError('positive integer minimum required')
    buckets=defaultdict(list)
    for x,z,p in points:
        if not all(math.isfinite(v) for v in (x,z,p)):raise ValueError('nonfinite point')
        if abs(p)<=max_perpendicular_m:buckets[math.floor(x/bin_width_m)].append(z)
    return [dict(bin_index=k,distance_along_m=(k+.5)*bin_width_m,
                 median_elevation_m=median(v),point_count=len(v),
                 elevation_iqr_m=quantile(v,.75)-quantile(v,.25))
            for k,v in sorted(buckets.items()) if len(v)>=min_points]


def quantile(values,q):
    s=sorted(values);index=(len(s)-1)*q;lo=math.floor(index);hi=math.ceil(index)
    return s[lo]*(hi-index)+s[hi]*(index-lo) if hi!=lo else s[lo]


def summarize(root):
    path=root/'published_profiles'/'cliff_profile_point_clouds.csv'
    groups=load_points(path);rows=[];binned={}
    for key,points in sorted(groups.items()):
        bins=bin_profile(points);binned['|'.join(map(str,key))]=bins
        rows.append(dict(glacier=key[0],cliff_id=key[1],profile_id=key[2],survey_date=key[3],
                         points=len(points),retained_bins=len(bins),
                         raw_distance_min_m=min(p[0] for p in points),raw_distance_max_m=max(p[0] for p in points),
                         raw_elevation_min_m=min(p[1] for p in points),raw_elevation_max_m=max(p[1] for p in points)))
    byprofile=defaultdict(list)
    for key in groups:byprofile[key[:3]].append(key[3])
    result=dict(source_doi='10.5281/zenodo.13073856',source_file_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                evidence='published corrected point clouds; observed geometry only',
                point_count=sum(len(p) for p in groups.values()),
                profile_count=len(byprofile),profile_epoch_count=len(groups),
                processing=dict(bin_width_m=.5,max_perpendicular_distance_m=.5,min_points_per_bin=5),
                caveats=['No stream forcing or water temperature in this file.',
                         'No cliff-face, debris or water classification in this schema.',
                         'These are cross-section coordinates, not a 3D melt volume.',
                         'Within-bin IQR is spatial scatter, not survey error or a confidence interval.',
                         'No quantitative melt attribution without delineation, survey uncertainty and matching footprint.'],
                profiles=[dict(glacier=k[0],cliff_id=k[1],profile_id=k[2],dates=sorted(set(v))) for k,v in sorted(byprofile.items())],
                epochs=rows)
    (root/'observed_geometry_audit.json').write_text(json.dumps(result,indent=2)+'\n')
    (root/'observed_profile_bins.json').write_text(json.dumps(binned,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('point_count','profile_count','profile_epoch_count','processing')},indent=2))
    return result,binned


def plot(root,binned):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    profile_ids=sorted({tuple(k.split('|')[:3]) for k in binned})
    # All profile IDs, no outcome-based selection.
    ncol=2;nrow=math.ceil(len(profile_ids)/ncol)
    fig,axes=plt.subplots(nrow,ncol,figsize=(12,3.4*nrow),squeeze=False)
    for ax,profile in zip(axes.flat,profile_ids):
        for key,bins in binned.items():
            if tuple(key.split('|')[:3])!=profile or not bins:continue
            ax.plot([b['distance_along_m'] for b in bins],[b['median_elevation_m'] for b in bins],label=key.split('|')[3],lw=1.1)
        ax.set_title(f'{profile[0]} cliff {profile[1]}, profile {profile[2]}')
        ax.set_xlabel('Distance along profile (m)');ax.set_ylabel('Corrected elevation (m a.s.l.)')
        ax.grid(alpha=.2);ax.legend(fontsize=6,ncol=2)
    for ax in list(axes.flat)[len(profile_ids):]:ax.set_visible(False)
    fig.suptitle('Observed Himalayan terrain profiles — not inferred melt volumes',fontsize=14)
    fig.tight_layout(rect=(0,0,1,.98));fig.savefig(root/'himalayan_observed_profiles.png',dpi=170)
    fig.savefig(root/'himalayan_observed_profiles.pdf');plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--plot',action='store_true')
    args=p.parse_args();r,b=summarize(args.root)
    if args.plot:plot(args.root,b)
