# Hydraulic input quality and interpretation

These flags apply to V5 output cells; they are not fitted parameters or an assertion of recovered missing observations.

| Input | Meaning and quality flag |
|---|---|
| `discharge_channel_full_daily.csv` | Author-derived P3 channel discharge. Rating-curve/inferred forcing; reconstructed during the13–24July pressure-data gap using a lake-balance assumption. Not a direct continuous salt-dilution record. |
| `stage3_full_daily.csv` | Author-processed P3 depth proxy, including reconstructed hydraulic support. Not a continuous independent manual floor-depth observation. |
| `stage5_daily.csv` | Author-processed pressure-derived depth proxy; NaN days refused. Finite daily means may use partial native support. Sensor pressure measures submergence, which equals full local depth only while the sensor sits on the floor. |
| Conserved Q at P5 | Assumes no substantial intervening sources, losses or meltwater-feedback correction between P5 and P3. This is a model assumption, not a colocated P5 discharge measurement. |
| Shared width1.5/2.0/2.5m | Segment-width engineering scenarios based on2±0.5m study assumption; not continuously measured local width at both stations. |

`Temperature_and_Stage.ipynb` cell91 marks P5 sensor “Not on bottom” and “in the air”; cell104 removes specified pressure intervals separately from temperature. Cells100–105 compare manual depths and processed pressure without demonstrating a continuous manual-height correction. The V5 reconstruction therefore remains conditional on pressure-derived daily stage as a geometry proxy. No vertical sensor-position bias is silently declared fixed.

Source uncertainty labels must remain intact: original manual bed files call their third column standard deviation; sums of two epoch errors used in a robustness scenario are not confidence intervals. None of these proxy quantities establishes prognostic hydraulic accuracy. V5 tests only how physically localising the imposed heat-transfer geometry changes an already-exposed thermal-component experiment.
