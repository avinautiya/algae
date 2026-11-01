# V5 local-hydraulics development reanalysis

Frozen after inspecting V3/V4 scores and before V5 scores. This is an explicitly target-exposed **development reanalysis**, not new held-out validation. Its motivation is the diagnosed assumption that one P3 Reynolds series and one transfer length impose the same heat-transfer coefficient at P3 and P5 despite different measured local depths. No local gain is fitted.

## Sole physical feature change

For each day and station s, impose original author daily Q from `discharge_channel_full_daily.csv`, depth h_s from `stage3_full_daily.csv` or `stage5_daily.csv`, and a shared rectangular width w. Derive

`S_s=w*h_s; v_s=Q/S_s; D_s=4*w*h_s/(w+2*h_s); Re_s=v_s*D_s/nu; h_transfer_s=k*Nu(Re_s,Pr)/D_s`.

Use nu=1.8e−6m²/s, Pr=13.5,k=0.56W/(mK),ice density917kg/m³,latent334000J/kg. These prescribed0°C engineering constants are unchanged or explicitly added as a fixed property, not fit to test outcomes. Width scenarios1.5/2.0/2.5m come from the study's2±0.5m segment-width assumption; the same w applies at both stations. No station-specific width calibration, temperature adjustment, gain or parameter is allowed.

Imposing downstream Q at upstream P5 assumes negligible local contributions/losses and meltwater-feedback changes along this segment. Daily quantities are author processed products; Q and P3 stage include reconstructed periods, not independent daily direct hydraulic observations. Pressure and discharge source clocks inherit the existing unresolved processed-file convention; manual-window shifts −2/0/+2h remain explicit scenarios. Daily Q and stage are held constant over calendar days. No hydraulic prognostic accuracy is claimed.

All source files must match SOURCE_MANIFEST.json hashes. Missing/NaN/nonpositive Q or depth rejects the interval across all methods. All standard Re/Pr gates remain. Do not infer a missing local depth from the melt target or replace missing pressure with a convenient constant.

**Pre-execution quality audit:** Temperature_and_Stage notebook cells91 and104 mark P5 periods with a sensor above the bed or in air; the code deletes pressure intervals separately from temperature. Cells100–105 compare processed pressure with manual depth checks but do not demonstrate a continuous manual-depth correction. Consequently finite daily pressure-derived stage is only a **conditional local depth proxy**, potentially averaging partial valid support. Do not call it continuously observed full water depth or assert all sensor-height error is removed. Preserve source quality flags and disclose that the nonlinear hydraulic closure is evaluated using daily stage means, not the mean of instantaneous h-transfer. This is a stage-proxy geometry experiment; its sensitivity cannot recover unmeasured sensor elevation above the floor.

## Everything else retained

Retain the same six thermal families, common per-scenario case membership, pooled one-scalar training-only fits, original per-station25July split, original endpoints/errors, nonnegative temperature clamp and all reporting rules. Native-time temperature reconstruction uses the independently selected linear endpoint bridge primary with left/right persistence sensitivity, at most24h finite-endpoint gaps and no extrapolation. All27 clock/temperature/width scenarios are executed for each imputer,81total. No physical closure or feature family will be selected from these81runs by minimum test error.

Report every nominal per-station and pooled result and all scenarios. Compare against both fair trained simple baselines and every unfitted comparator. Also compare V4 saved predictions on **identical accepted target intervals** to isolate local-hydraulics effects from membership changes. Record daily D_s,Re_s,v_s and the input data used. Local-hydraulics geometry is an existing physical construction; an error reduction is not itself a new theorem or proof of full coupled cliff-model field accuracy. Any superiority found here is developmental evidence requiring independent prospective data.
