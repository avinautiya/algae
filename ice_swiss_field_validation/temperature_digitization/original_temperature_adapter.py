"""Audit and normalize original-author corrected temperature products.

Never silently assigns a timezone or fills missing values. The raw-suffix
product is audited separately because its final additive offset is absent.
"""
import argparse, collections, csv, datetime as dt, hashlib, json, math, re
from pathlib import Path

OFFSETS={3:.035,5:.04}

def parse_time(text):
    values=re.fullmatch(r'(\d{4})-(\d{1,2})-(\d{1,2})T(\d{1,2}):(\d{1,2}):(\d{1,2})',text)
    if not values:raise ValueError('unsupported original naive timestamp: '+text)
    return dt.datetime(*map(int,values.groups()))

def source_rows(path):
    with path.open(newline='') as f:
        reader=csv.DictReader(f)
        if reader.fieldnames[0]!='Date_Time' or len(reader.fieldnames)!=2:raise ValueError('unexpected schema')
        value_name=reader.fieldnames[1]
        for i,row in enumerate(reader,start=2):
            yield i,row['Date_Time'],parse_time(row['Date_Time']),float(row[value_name])

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def verified(root,relative,manifest):
    entries=[e for e in manifest['members'] if e['path']==relative]
    if len(entries)!=1 or not entries[0]['verified_crc32']:raise ValueError('source member not CRC verified')
    path=root/relative
    if sha(path)!=entries[0]['sha256']:raise ValueError('source member SHA mismatch: '+relative)
    return path

def normalize(root,out):
    manifest=json.loads((root/'SOURCE_MANIFEST.json').read_text());out.mkdir(parents=True,exist_ok=True)
    all_rows=[];audit={'timezone':'not established by temperature headers; no UTC conversion performed',
        'source_manifest_sha256':sha(root/'SOURCE_MANIFEST.json'),'products':{},'raw_vs_daily':{},
        'sensor_accuracy_bound_C':.05,'instrument_source':'https://tc.copernicus.org/articles/15/5133/2021/ section 3.1.2',
        'temperature_representation':'point sensor at floor; unknown offset to bulk water',
        'no_missing_value_imputation':True,'final_offsets_already_in_daily_and_hourly':OFFSETS}
    for station,offset in OFFSETS.items():
        daily_values={}
        for period,suffix in [('daily','daily'),('hourly','hour')]:
            rel=f'products/temp{station}_{suffix}.csv';p=verified(root,rel,manifest);times=[];missing=[];finite=[]
            for line,original,time,value in source_rows(p):
                times.append(time);isvalid=math.isfinite(value)
                if not isvalid:missing.append(original)
                else:finite.append(value)
                if period=='daily':
                    if time.time()!=dt.time(12):raise ValueError('daily label is not noon')
                    start=time.replace(hour=0);end=start+dt.timedelta(days=1);daily_values[time.date()]=value
                else:
                    if time.minute or time.second:raise ValueError('hourly label is not hour start')
                    start=time;end=start+dt.timedelta(hours=1)
                all_rows.append(dict(station=f'P{station}',period=period,original_timestamp=original,
                    time_label_naive=time.isoformat(),window_start_naive=start.isoformat(),window_end_naive=end.isoformat(),
                    timezone='unresolved',temperature_C=value if isvalid else '',valid=isvalid,
                    source_path=rel,source_csv_line=line,source_sha256=sha(p),instrument_accuracy_bound_C=.05,
                    representation='floor sensor, not measured bulk mean',final_additive_offset_already_applied_C=offset))
            if len(set(times))!=len(times) or times!=sorted(times):raise ValueError('duplicate or nonmonotone time')
            audit['products'][rel]=dict(rows=len(times),valid=len(finite),missing=len(missing),
                first=times[0].isoformat(),last=times[-1].isoformat(),
                finite_min_C=min(finite),finite_max_C=max(finite),missing_timestamps=missing,
                gaps=[dict(start=a.isoformat(),end=b.isoformat()) for a,b in zip(times,times[1:])
                    if b-a!=(dt.timedelta(days=1) if period=='daily' else dt.timedelta(hours=1))])
        rawrel=f'products/temp{station}_corr_raw.csv';raw=verified(root,rawrel,manifest)
        grouped=collections.defaultdict(list);rawcount=collections.Counter();rawvalid=collections.Counter();previous=None;first_per_date={}
        for line,original,time,value in source_rows(raw):
            rawcount[time.date()]+=1
            if time.date() not in first_per_date:first_per_date[time.date()]=time
            if math.isfinite(value):grouped[time.date()].append(value);rawvalid[time.date()]+=1
        comparisons=[]
        for date,observed in sorted(daily_values.items()):
            vals=grouped[date]
            if not vals or not math.isfinite(observed):continue
            rawmean=math.fsum(vals)/len(vals)
            comparisons.append(dict(date=date.isoformat(),raw_rows=rawcount[date],valid_raw=rawvalid[date],
                mean_of_all_finite_exported_raw_C=rawmean,original_corrected_daily_C=observed,
                residual_if_raw_has_no_final_offset_C=observed-(rawmean+offset),
                residual_if_raw_already_final_corrected_C=observed-rawmean))
        audit['raw_vs_daily'][f'P{station}']=dict(final_offset_from_source_C=offset,
            max_abs_residual_C=max(abs(r['residual_if_raw_has_no_final_offset_C']) for r in comparisons),
            note='raw export has already been rounded; helper omissions/partial days affect means. Direct average is an audit, not a replacement product.',
            daily=comparisons)
    with (out/'original_corrected_temperature_normalized.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(all_rows[0]));w.writeheader();w.writerows(all_rows)
    (out/'ORIGINAL_TEMPERATURE_AUDIT.json').write_text(json.dumps(audit,indent=2)+'\n')
    return audit

def main():
    p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    audit=normalize(a.root,a.out)
    print(json.dumps({k:{x:y for x,y in v.items() if x!='missing_timestamps'} for k,v in audit['products'].items()},indent=2))

if __name__=='__main__':main()
