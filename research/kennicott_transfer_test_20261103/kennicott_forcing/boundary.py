"""Source-verified thermistor-plane boundary, not radiometric surface temperature.
Only temperature/weather features are read. No melt target reader exists here.
"""
import bisect,csv,datetime as dt,hashlib,json,math
from pathlib import Path
DATA=Path(__file__).resolve().parent.parent/'data'/'kennicott'
PROFILE={'KDT_2020b':(.05,'Depth: 5 cm'),'KDT_2023a':(.02,'Depth: 2 cm'),'KDT_2023b':(.09,'Depth: 9.0 cm'),'KDT_2023c':(.02,'Depth: 2.0 cm')}
_CACHE={}
def _datetime(value):
 if isinstance(value,dt.datetime):out=value
 else:out=dt.datetime.fromisoformat(value.replace("Z","+00:00"))
 if out.tzinfo is not None:
  if out.utcoffset()!=dt.timedelta(0):raise ValueError('Use UTC timestamps only')
  out=out.replace(tzinfo=None)
 return out

def _seconds(value):return (_datetime(value)-dt.datetime(1970,1,1)).total_seconds()
def _label(t):return (dt.datetime(1970,1,1)+dt.timedelta(seconds=t)).isoformat()+'Z'
def _load(name,clock,col):
 manifest=json.loads((DATA/'FEATURE_EXPORT_MANIFEST.json').read_text());row=next(x for x in manifest['files'] if x['local_path']=='features/'+name+'.csv');raw=(DATA/row['local_path']).read_bytes()
 if hashlib.sha256(raw).hexdigest()!=row['sha256']:raise ValueError('Feature SHA mismatch')
 key=(name,col,row['sha256'])
 if key in _CACHE:return _CACHE[key]
 rows=list(csv.DictReader(raw.decode('utf-8-sig').splitlines()));times=[];values=[]
 for r in rows:
  try:v=float(r[col])
  except ValueError:continue
  if not math.isfinite(v):continue
  times.append(_seconds(r[clock]));values.append(v)
 if len(times)<2 or any(b<=a for a,b in zip(times,times[1:])):raise ValueError('Strict native finite clock needed')
 _CACHE[key]=(times,values,row['sha256']);return _CACHE[key]
def _at(times,values,t,max_gap_s):
 i=bisect.bisect_left(times,t)
 if i<len(times) and times[i]==t:return values[i]
 if i==0 or i==len(times):raise ValueError('Feature extrapolation forbidden')
 a,b=times[i-1],times[i]
 if b-a>max_gap_s:raise ValueError('Feature interpolation gap exceeds contract')
 w=(t-a)/(b-a);return (1-w)*values[i-1]+w*values[i]
def ols(x,y):
 if len(x)!=len(y) or len(x)<2 or not all(math.isfinite(v) for v in x+y):raise ValueError('finite paired regression data needed')
 mx=sum(x)/len(x);my=sum(y)/len(y);den=sum((v-mx)**2 for v in x)
 if den<=1e-20:raise ValueError('Air variance too small')
 slope=sum((a-mx)*(b-my) for a,b in zip(x,y))/den;return my-slope*mx,slope

def build_boundary(profile_id,start_utc,end_utc,*,air_station='Debris_T2',proxy_error_multiplier=0.,installation_days=5,fit_days=14):
 if profile_id not in PROFILE or air_station not in {'Debris_T2','Debris_T3'}:raise ValueError('Profile/air source not allowlisted')
 if not math.isfinite(proxy_error_multiplier):raise ValueError('finite error multiplier needed')
 if installation_days!=5 or fit_days!=14:raise ValueError('Fixed pre-outcome protocol requires5/14days')
 start,end=_seconds(start_utc),_seconds(end_utc)
 if end<=start or (end-start)%3600:raise ValueError('Positive integer-hour horizon required')
 plane,col=PROFILE[profile_id];tt,tv,th=_load(profile_id,'Datetime',col);at,av,ah=_load(air_station,'Date_Time','T_1_75m')
 reset_end=tt[0]+5*86400;fit_end=reset_end+14*86400;fit_times=[reset_end+i*3600 for i in range(337)]
 x=[_at(at,av,t,7200) for t in fit_times];y=[_at(tt,tv,t,600) for t in fit_times]
 a7,b7=ols(x[:168],y[:168]);err=[a7+b7*v-w for v,w in zip(x[168:],y[168:])];rmse=math.sqrt(sum(v*v for v in err)/len(err));a,b=ols(x,y)
 times=[start+i*3600 for i in range(int((end-start)/3600)+1)];temps=[];flags=[];air_values=[]
 for t in times:
  use_proxy=t<reset_end
  air_value=_at(at,av,t,7200);air_values.append(air_value)
  temp=a+b*air_value+proxy_error_multiplier*2*rmse if use_proxy else _at(tt,tv,t,600)
  temps.append(temp);flags.append(use_proxy)
 return {'times_utc':[_label(t) for t in times],'elapsed_s':[t-start for t in times],'boundary_C':temps,'proxy_boundary_C':[a+b*v+proxy_error_multiplier*2*rmse for v in air_values],'air_C':air_values,'is_proxy':flags,'metadata':{'profile_id':profile_id,'air_station':air_station,'upper_boundary_plane_m':plane,'boundary_quantity':'temperature at first measured buried thermistor plane, not actual surface','native_clock':'UTC from source metadata','resampling':'hourly sample labels; linear interpolation, no averaging or clipping','proxy_error_multiplier':proxy_error_multiplier,'proxy_offset_C':proxy_error_multiplier*2*rmse,'conditional_scenario_not_confidence_interval':True,'fit_intercept_C':a,'fit_slope':b,'fit_rows':337,'fit_start_UTC':_label(reset_end),'fit_end_UTC':_label(fit_end),'blocked_fit_rows':168,'blocked_test_rows':169,'blocked_RMSE_C':rmse,'blocked_bias_C':sum(err)/len(err),'blocked_max_abs_C':max(abs(e) for e in err),'initial_profile_UTC':_label(tt[0]),'proxy_until_UTC_exclusive':_label(reset_end),'proxy_hours':sum(flags[:-1]),'total_hours':len(times)-1,'source_profile_SHA256':th,'source_air_SHA256':ah,'uses_melt_targets':False,'warning':'Temperature-conditioned retrospective reconstruction; future thermistors used. Not an advance weather-only forecast. Source missing initial period retained.'}}
