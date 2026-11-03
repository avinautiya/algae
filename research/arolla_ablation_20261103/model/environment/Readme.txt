Codes:
1) Functions for thermal conductivity estimations:
	a) RunThermalConductivity_allSites (to calculate CRh or CRi) -> in here choose the parameters to consider (sensor uncertainty, soil moisture, lithological properties) and the method to use (CRh or CRi)
	b) ThermalConductivity_CRh
	c) ThermalConductivity_CRi
	d) Conductivity_Brock (to calculate NYB)
	e) Calibrate_deb_par_PIR_k_z0 (calculates optimised value of k (and z0))

2) Functions for aerodynamic surface roughness lenght:
	a) Runz0 (main) -> in here choose the parameters (d or no d, surface temperature, sensor uncertainty, period of measurements)
	b) z0Calc
	c) Calibrate_deb_par_PIR_k_z0 -> it calibrates z0 if chosen


Data:
1) Ablation_Tower1.mat, Ablation_Tower2.mat and Ablation_Tower3.mat: Timetables containing the automatic distances measured from the ablation stake pictures in m and the calculated ablation in mm w.e.

2) DATA_Piramide.mat: 
	a) debtemp: structure containing timetables of thermistor data (°C) for the three sites, where DT1 is the closest to the surface and DT5 the closest to the ice.
	b) sm: structure containing timetables of moisture content (m3/m3) data for the three sites, where SM1 is the closest to the surface and SM3 the closest to the ice.
	c) trh: structure containing timetables of temperature (°C) and relative humidity (%) data for the three sites, where XX_1, XX_2, XX_3 and XX_4 correspond to the variables at 0.5, 1, 2 and 2.7 m above the surface. TA is the air temperature, RH is the relative humidity and DTA is the dew temperature.
	d) wind: structure containing timetables of wind speed (m/s) and direction (°) data for the three sites, where L1, L2, L3 and L4 correspond to 0.5, 1, 2 and 2.7 m above the surface and FF is the wind speed, FFmx is the maximum wind speed, FFstd is the standard deviation, and DIR is the direction.

3) LithologyData: structure containing the density (kg/m3) and specific heat capacity (J/kg/K) for each site.

4) SM_data: structure containing the moisture content (m3/m3) for the selected period, as well as the mean value for each site

5) Calibration: data required to run the k and z0 calibrations with T&C


