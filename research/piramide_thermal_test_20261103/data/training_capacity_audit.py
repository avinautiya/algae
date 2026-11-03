"""Feature-only TRAIN interval heat-capacity sensitivity, never camera outcomes.

Training is the predeclared half-open interval [Dec2 03, Dec16 03).
Native moisture samples outside [0,phi] are flagged/excluded from proxy means,
never clipped. Proxy means are NOT measured whole-column moisture.
"""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
from read_verified_process import read_process
ROOT=Path(__file__).parent
START='2022-12-02T03:00';END='2022-12-16T03:00';PHI=.44

def compute():
    sm=read_process(ROOT/'decoded','sm','Tower1');temps=read_process(ROOT/'decoded','debtemp','Tower1')
    sm=sm.loc[(sm.index>=START)&(sm.index<END)];temps=temps.loc[(temps.index>=START)&(temps.index<END)]
    assert len(sm)==len(temps)==672
    info={}
    for col in sm:
        good=sm[col].between(0,PHI)&sm[col].notna()
        info[col]={'rows':len(sm),'valid':int(good.sum()),'negative':int((sm[col]<0).sum()),'above_assumed_porosity':int((sm[col]>PHI).sum()),
            'mean_valid_samples':float(sm.loc[good,col].mean()) if good.any() else None,'unclipped_min':float(sm[col].min()),'unclipped_max':float(sm[col].max())}
    allgood=((sm>=0)&(sm<=PHI)&sm.notna()).all(axis=1)
    theta_proxy=float(sm.loc[allgood].mean(axis=1).mean()) if allgood.any() else None
    rho=2684.9333333333334;cp=796.8666666666668;rho_w_cw=1000*4180.;rho_air_cp=1.2*1005.
    rock=(1-PHI)*rho*cp
    def C(theta):return rock+rho_w_cw*theta+rho_air_cp*(PHI-theta)
    scenes={'dry_assumed_pore_air':{'theta':0.,'C_J_m3K':C(0.)},'training_valid_probe_proxy_not_measured_depth_mean':{'theta':theta_proxy,'C_J_m3K':C(theta_proxy) if theta_proxy is not None else None},'fully_water_saturated_conditional':{'theta':PHI,'C_J_m3K':C(PHI)}}
    return {'scope':'Training-only environmental features; no camera values, melt predictions, fitting or scores',
       'training_interval_half_open_source_local':[START,END],'rows':672,'porosity_assumed_not_measured':PHI,
       'source_files':{n:hashlib.sha256((ROOT/'decoded'/n).read_bytes()).hexdigest() for n in ['sm_Tower1.csv','debtemp_Tower1.csv']},
       'moisture_qc':info,'all_three_valid_rows':int(allgood.sum()),'proxy_definition':'Unweighted mean across three probes per complete physically valid row, then mean across TRAIN rows; not a volume-weighted field measurement',
       'temperature_negative_samples':{n:int((temps[n]<0).sum()) for n in temps},
       'sensible_C_model':'(1-phi)*rho_rock*cp_rock + theta*rho_water*cp_water + (phi-theta)*rho_air*cp_air',
       'constants':{'rho_rock_kg_m3':rho,'cp_rock_J_kg_K':cp,'water_rho_cp_J_m3K':rho_w_cw,'air_rho_cp_J_m3K':rho_air_cp},
       'finite_sensible_capacity_scenes':scenes,
       'maximum_pore_phase_latent_J_m3_assuming_water_volume_phi':PHI*1000*334000,
       'phase_warning':'Finite sensible C bounds do not bound apparent heat capacity at freezing. Pore-ice latent storage requires enthalpy/phase state. Cold-glacier-ice admittance alone does not resolve wet-debris freezing.',
       'melt_validation_ready':False}
if __name__=='__main__':
    (ROOT/'TRAINING_CAPACITY_AUDIT.json').write_text(json.dumps(compute(),indent=2)+'\n')
    print('Training-only capacity audit saved; no camera outcomes accessed.')
