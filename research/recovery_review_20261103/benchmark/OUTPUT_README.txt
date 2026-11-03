Model outputs from the Debris-Covered Glacier Melt Model Intercomparison Experiment.

The DCG-MIP_Model_outputs.zip file contains all model outputs, sorted in folders by the model acronym, and the site acronym, as well as the experiments (EXP1, EXP2, EXP3, EXP4) as per the manuscript. 

The model output files are organised in matrices, where rows are time steps (either hourly or daily),
and columns are the debris thickness of the model run. 

Column one contains always the actual debris at the site. The other 12 columns are the model
runs with variable debris thickness, according to experiment 2:
The sites thicknesses of column 1 are: [6 10 61 30 20 18 11 80 30] in centimetres. 
The thickness for columns 2 to 13 are: [1 2 4 6 8 10 12 15 20 30 50 100] in centimetres. 

File naming: {MODELNAME}_{SITENAME}_{EXPX}_{U}.csv

Where models names are: {'GRO17B','GRO17A','DEBCF','ROU15','AMelt','THRED','d2EB', ...
    'DEBPG','MCC19','DETImC','DETImU','KP1','KM1','KO2', ...
    'HyperfitU','HyperfitC','DDFdebrisC'};

Site names are: {'ARO','CN','DJA','LIR','MIA','PIR','SDF','TAP','TAS'};

And there are four experiments: EXP1, EXP2, EXP3, EXP4. 

For EXP3 and EXP4, there are 100 files, from 001 to 100, corresponding to the uncertainty simulations.

The folder Dates_Rows contains an array with hourly and daily time steps for each site. 

See DCG-MIP manuscript for more details. 

