"""Explicit retrospective reconstruction alongside unchanged strict support.

Only finite valid endpoint gaps <=24 h may be reconstructed. Observed and
reconstructed heat integrals/durations remain separate. No extrapolation.
"""
import json,math
from pathlib import Path
from native_temperature import TemperatureSeries,seconds
from original_temperature_adapter import source_rows,verified,OFFSETS

METHODS=('linear_endpoint_bridge','constant_previous_endpoint','constant_next_endpoint')

class ExplicitSegments(TemperatureSeries):
    """Use declared disjoint segments, without inventing links between them."""
    def __init__(self,segments):
        self.segments=list(segments);self.starts=[];self.ends=[]
        self.prefix_integral=[0.];self.prefix_coverage=[0.];self._positive_prefix_cache={}
        previous_end=None
        for a,b,ta,tb in self.segments:
            if not all(math.isfinite(x) for x in [a,b,ta,tb]) or b<=a:raise ValueError('invalid segment')
            if previous_end is not None and a<previous_end:raise ValueError('overlapping segments')
            previous_end=b;self.starts.append(a);self.ends.append(b)
            self.prefix_integral.append(self.prefix_integral[-1]+(b-a)*(ta+tb)/2)
            self.prefix_coverage.append(self.prefix_coverage[-1]+(b-a))

class GapReconstructedTemperatureSeries:
    def __init__(self,records,max_observed_gap_seconds=60,max_reconstruction_gap_seconds=86400):
        if not math.isfinite(max_observed_gap_seconds) or not 0<max_observed_gap_seconds<=60:
            raise ValueError('observed adjacency limit must be positive and no larger than audited 60 seconds')
        if (not math.isfinite(max_reconstruction_gap_seconds) or
            max_reconstruction_gap_seconds<=0 or max_reconstruction_gap_seconds>86400):
            raise ValueError('reconstruction gap limit must be positive and no larger than 24 hours')
        records=list(records)
        self.observed=TemperatureSeries(records,max_observed_gap_seconds)
        self.max_observed_gap_seconds=max_observed_gap_seconds
        self.max_reconstruction_gap_seconds=max_reconstruction_gap_seconds
        candidate_gaps=[];last_finite=None;self.rejected_endpoint_gaps=[]
        for index,(time,value) in enumerate(records):
            if not math.isfinite(value):continue
            now=seconds(time)
            if last_finite:
                previous_index,a,ta=last_finite;duration=now-a
                # An original consecutive finite pair <=strict threshold is
                # already observed. Other finite endpoint pairs span missing
                # values or a recorded acquisition gap.
                strict_pair=(index==previous_index+1 and duration<=max_observed_gap_seconds)
                if not strict_pair:
                    if duration<=max_reconstruction_gap_seconds:
                        candidate_gaps.append((a,now,ta,value))
                    else:self.rejected_endpoint_gaps.append(dict(start_seconds=a,end_seconds=now,duration_seconds=duration))
            last_finite=index,now,value
        self.reconstructed={
            'linear_endpoint_bridge':ExplicitSegments(candidate_gaps),
            'constant_previous_endpoint':ExplicitSegments([(a,b,ta,ta) for a,b,ta,tb in candidate_gaps]),
            'constant_next_endpoint':ExplicitSegments([(a,b,tb,tb) for a,b,ta,tb in candidate_gaps])}
        self.all_support=ExplicitSegments(sorted([(a,b,0.,0.) for a,b,ta,tb in self.observed.segments+candidate_gaps]))
        self.reconstructible_endpoint_gaps=candidate_gaps

    @classmethod
    def from_archive(cls,root,station,max_observed_gap_seconds=60,max_reconstruction_gap_seconds=86400):
        station=int(str(station).removeprefix('P'))
        if station not in OFFSETS:raise ValueError('only P3/P5 audited offsets supported')
        root=Path(root);manifest=json.loads((root/'SOURCE_MANIFEST.json').read_text())
        relative=f'products/temp{station}_corr_raw.csv';path=verified(root,relative,manifest)
        result=cls(((t,v+OFFSETS[station]) for _,_,t,v in source_rows(path)),
                   max_observed_gap_seconds,max_reconstruction_gap_seconds)
        result.station=f'P{station}';result.final_offset_C=OFFSETS[station];result.source_path=relative
        return result

    def _integrate(self,start,end,extra_offset_C,method,positive):
        if method not in METHODS:raise ValueError('undeclared reconstruction method')
        operation='integrate_positive' if positive else 'integrate'
        observed=getattr(self.observed,operation)(start,end,extra_offset_C)
        filled=getattr(self.reconstructed[method],operation)(start,end,extra_offset_C)
        support=self.all_support.integrate(start,end)
        requested=support['requested_seconds']
        result=dict(integral_C_s=observed['integral_C_s']+filled['integral_C_s'],
            observed_integral_C_s=observed['integral_C_s'],reconstructed_integral_C_s=filled['integral_C_s'],
            observed_seconds=observed['covered_seconds'],reconstructed_seconds=filled['covered_seconds'],
            covered_seconds=observed['covered_seconds']+filled['covered_seconds'],requested_seconds=requested,
            observed_fraction=observed['covered_seconds']/requested,reconstructed_fraction=filled['covered_seconds']/requested,
            coverage_fraction=support['coverage_fraction'],max_uncovered_gap_seconds=support['max_uncovered_gap_seconds'],
            reconstruction_method=method,extra_offset_C=extra_offset_C,
            max_reconstruction_gap_seconds=self.max_reconstruction_gap_seconds,
            no_extrapolation=True,source_clock='naive/unresolved',
            representation='measured floor-point plus explicitly reconstructed gaps; not observed bulk temperature',
            melt_only_temperature_clamp=positive)
        if positive:result['unclipped_signed_integral_C_s']=observed['unclipped_signed_integral_C_s']+filled['unclipped_signed_integral_C_s']
        return result

    def integrate(self,start,end,extra_offset_C=0,method='linear_endpoint_bridge'):
        return self._integrate(start,end,extra_offset_C,method,False)

    def integrate_positive(self,start,end,extra_offset_C=0,method='linear_endpoint_bridge'):
        return self._integrate(start,end,extra_offset_C,method,True)
