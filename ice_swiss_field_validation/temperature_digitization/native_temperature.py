"""Time-weighted integrals of author-native, quality-masked floor temperatures.

No timezone is inferred. Only adjacent finite endpoints separated by at most
max_gap_seconds support a piecewise-linear segment. Missing support is reported.
"""
import argparse,bisect,csv,datetime as dt,json,math,collections
from pathlib import Path
from original_temperature_adapter import source_rows,verified,OFFSETS

EPOCH=dt.datetime(2000,1,1)

def seconds(time):
    if time.tzinfo is not None:raise ValueError('expected explicit naive source-clock time; reconcile timezone upstream')
    return (time-EPOCH).total_seconds()

class TemperatureSeries:
    def __init__(self,records,max_gap_seconds=60):
        if not math.isfinite(max_gap_seconds) or max_gap_seconds<=0:raise ValueError('positive finite max gap required')
        self.max_gap_seconds=max_gap_seconds;self.segments=[];self.starts=[];self.ends=[];self._positive_prefix_cache={}
        self.prefix_integral=[0.];self.prefix_coverage=[0.];self.step_histogram=collections.Counter();self.invalid_pairs=0
        previous=None;self.first=None;self.last=None
        for time,value in records:
            t=seconds(time)
            if self.first is None:self.first=time
            self.last=time
            if previous:
                a,ta=previous;gap=t-a
                if gap<=0:raise ValueError('timestamps must be unique and increasing')
                self.step_histogram[gap]+=1
                if gap<=max_gap_seconds and math.isfinite(ta) and math.isfinite(value):
                    self.starts.append(a);self.ends.append(t);self.segments.append((a,t,ta,value))
                    self.prefix_integral.append(self.prefix_integral[-1]+gap*(ta+value)/2)
                    self.prefix_coverage.append(self.prefix_coverage[-1]+gap)
                else:self.invalid_pairs+=1
            previous=t,value
        if not self.segments:raise ValueError('no supported finite adjacent-temperature segments')

    @classmethod
    def from_archive(cls,root,station,max_gap_seconds=60):
        station=int(str(station).removeprefix('P'))
        if station not in OFFSETS:raise ValueError('only independently audited P3/P5 final offsets supported')
        root=Path(root);manifest=json.loads((root/'SOURCE_MANIFEST.json').read_text())
        relative=f'products/temp{station}_corr_raw.csv';path=verified(root,relative,manifest)
        # The misleading corr_raw file does not contain the final additive
        # correction used in author daily/hourly products. Apply exactly once.
        result=cls(((time,value+OFFSETS[station]) for _,_,time,value in source_rows(path)),max_gap_seconds)
        result.station=f'P{station}';result.source_relative_path=relative;result.final_offset_C=OFFSETS[station]
        return result

    def _primitive(self,t):
        j=bisect.bisect_right(self.starts,t)-1
        if j<0:return 0.,0.
        a,b,ta,tb=self.segments[j];u=min(max(t-a,0),b-a)
        return (self.prefix_integral[j]+ta*u+(tb-ta)*u*u/(2*(b-a)),
                self.prefix_coverage[j]+u)

    def integrate(self,start,end,extra_offset_C=0):
        a,b=seconds(start),seconds(end)
        if not math.isfinite(extra_offset_C) or b<=a:raise ValueError('positive interval and finite sensitivity offset required')
        ia,ca=self._primitive(a);ib,cb=self._primitive(b);covered=cb-ca
        cursor=a;largest=0.;j=max(0,bisect.bisect_right(self.ends,a))
        while j<len(self.segments) and self.starts[j]<b:
            left=max(a,self.starts[j]);right=min(b,self.ends[j])
            if right>left:
                largest=max(largest,left-cursor);cursor=right
            j+=1
        largest=max(largest,b-cursor)
        return dict(integral_C_s=ib-ia+extra_offset_C*covered,
            covered_seconds=covered,requested_seconds=b-a,coverage_fraction=covered/(b-a),
            max_uncovered_gap_seconds=largest,
            extra_offset_C=extra_offset_C,source_clock='naive/unresolved',
            representation='piecewise-linear floor-point temperature, not observed bulk temperature')

    @staticmethod
    def _positive_area(a,b,duration):
        if a>=0 and b>=0:return duration*(a+b)/2
        if a<=0 and b<=0:return 0.
        if a>0:return duration*a*a/(2*(a-b))
        return duration*b*b/(2*(b-a))

    def _positive_primitive(self,t,shift):
        if shift not in self._positive_prefix_cache:
            prefix=[0.]
            for a,b,ta,tb in self.segments:
                prefix.append(prefix[-1]+self._positive_area(ta+shift,tb+shift,b-a))
            self._positive_prefix_cache[shift]=prefix
        j=bisect.bisect_right(self.starts,t)-1
        if j<0:return 0.
        a,b,ta,tb=self.segments[j];u=min(max(t-a,0),b-a)
        edge=ta+(tb-ta)*u/(b-a)
        return self._positive_prefix_cache[shift][j]+self._positive_area(ta+shift,edge+shift,u)

    def integrate_positive(self,start,end,extra_offset_C=0):
        """Explicit melt-only max(T+offset,0) scenario, not raw-data correction."""
        result=self.integrate(start,end,extra_offset_C)
        result['unclipped_signed_integral_C_s']=result['integral_C_s']
        result['integral_C_s']=self._positive_primitive(seconds(end),extra_offset_C)-self._positive_primitive(seconds(start),extra_offset_C)
        result['melt_only_temperature_clamp']=True
        return result

def audit(root,out,max_gap_seconds):
    out.mkdir(parents=True,exist_ok=True);summary={};hourly=[];daily=[]
    for station in [3,5]:
        s=TemperatureSeries.from_archive(root,station,max_gap_seconds)
        start=s.first.replace(hour=0,minute=0,second=0);finish=s.last.replace(hour=0,minute=0,second=0)+dt.timedelta(days=1)
        time=start
        while time<finish:
            row=s.integrate(time,time+dt.timedelta(hours=1));row.update(station=f'P{station}',start=time.isoformat(),end=(time+dt.timedelta(hours=1)).isoformat());hourly.append(row);time+=dt.timedelta(hours=1)
        author={time.date():v for _,_,time,v in source_rows(root/f'products/temp{station}_daily.csv')}
        time=start
        while time<finish:
            row=s.integrate(time,time+dt.timedelta(days=1));mean=row['integral_C_s']/row['covered_seconds'] if row['covered_seconds'] else None
            observed=author.get(time.date(),float('nan'));row.update(station=f'P{station}',date=time.date().isoformat(),
                supported_time_weighted_mean_C=mean,author_sample_weighted_daily_C=observed if math.isfinite(observed) else '',
                author_minus_timeweighted_C=observed-mean if math.isfinite(observed) and mean is not None else '')
            daily.append(row);time+=dt.timedelta(days=1)
        summary[f'P{station}']=dict(final_offset_C=s.final_offset_C,segments=len(s.segments),
            rejected_pairs=s.invalid_pairs,max_gap_seconds=max_gap_seconds,
            source_time_step_histogram_seconds=dict(s.step_histogram),
            first=s.first.isoformat(),last=s.last.isoformat())
    for name,rows in [('NATIVE_HOURLY_INTEGRALS.csv',hourly),('NATIVE_DAILY_COMPARISON.csv',daily)]:
        with (out/name).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (out/'NATIVE_INTEGRAL_AUDIT.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('--out',type=Path,required=True);p.add_argument('--max-gap-seconds',type=float,default=60);a=p.parse_args();audit(a.root,a.out,a.max_gap_seconds)
