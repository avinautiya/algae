# Observation operator and environmental interpretation

Let eta be the water-surface elevation, z_b the channel bed and z_s the pressure-sensor elevation, all in the same reference frame. Hydrostatic pressure head is

    p/(rho*g) = eta - z_s.

The actual full flow depth is

    h = eta - z_b = p/(rho*g) + (z_s - z_b).

Consequently pressure head is full depth only when the sensor is on the bed (or its clearance is independently known). An unknown positive sensor clearance admits multiple geometries with the same pressure record. Since area, hydraulic diameter and flow resistance depend on actual depth, a deterministic geometry/velocity prediction from pressure alone is then underdetermined. This is an elementary observation identity, not a new mathematical theorem or a glacier discovery.

This matters in the recovered archive. Temperature_and_Stage notebook cell79 identifies periods when the P3 sensor was not on the channel floor; cell80 removes those pressure records while retaining otherwise valid temperatures. The P5 plots also label an off-bottom period and subsequent air exposure. Quality-approved temperature, quality-approved pressure and independently measured bed elevation therefore have distinct support. A missing temperature or pressure must not be interpreted as zero water contact or zero flow.

Author stage3_full_daily and discharge_channel_full_daily contain reconstructions. Treating those as imposed forcing can support a conditional component hindcast; it does not validate prediction of the underlying geometry or discharge. Likewise a floor-temperature measurement with a documented accuracy bound is not an independently observed bulk-water temperature.

The hydraulic source-expression audit is separate: positive depths and floor elevations imply the piezometric drop (z_up+h_up)-(z_down+h_down). Correcting its downstream sign does not determine missing sensor clearance, kinetic-head changes or time-varying friction. The saved negative velocity-transfer result demonstrates why the algebraic correction alone cannot be advertised as a superior field model.

For a complete cliff model, atmospheric heat, water-contact evolution, debris transmission, floor motion and face retreat must remain distinct states and fluxes. A successful channel thermal-component hindcast would not, by itself, validate a vertical cliff's full melt volume or lifetime. The present package preserves that boundary explicitly.
