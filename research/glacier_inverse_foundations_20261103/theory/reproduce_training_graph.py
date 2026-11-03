"""Development-only diagnostics. Refuses overwrite; opens no new melt targets."""
import argparse,csv,datetime as dt,json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from paths import require_piramide
from weak_graph import fit,lowest_edge_enclosure
def main():
    a=argparse.ArgumentParser();a.add_argument('output',type=Path);x=a.parse_args()
    if x.output.exists():a.error('new output file required')
    old=require_piramide();src=old/'data/decoded/debtemp_Tower1.csv'
    with src.open()as f:rows=[r for r in csv.DictReader(f)if '2022-12-02T03:00:00'<=r['SOURCE_LOCAL_LABEL']<'2022-12-16T03:00:00']
    labels=[dt.datetime.fromisoformat(r['SOURCE_LOCAL_LABEL'])for r in rows]
    t=[(z-labels[0]).total_seconds()for z in labels];T=[[float(r['DT'+str(j)])for j in range(1,6)]for r in rows]
    C=json.loads((old/'data/TRAINING_CAPACITY_AUDIT.json').read_text())['finite_sensible_capacity_scenes']['training_valid_probe_proxy_not_measured_depth_mean']['C_J_m3K']
    z=[.02,.125,.25,.375,.5]
    result=dict(scope='Tower1 TRAIN development only',new_validation=False,Tower2_camera_values_opened=False,
        primary_window_hours=24,source_capacity_is_conditional=True,
        window_diagnostics=[fit(z,t,T,capacity_J_m3K=C,window_s=h*3600)for h in[6,12,24]],
        necessary_lowest_edge_enclosure=lowest_edge_enclosure(z,t,T,capacity_J_m3K=C,window_s=86400))
    with x.output.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
if __name__=='__main__':main()
