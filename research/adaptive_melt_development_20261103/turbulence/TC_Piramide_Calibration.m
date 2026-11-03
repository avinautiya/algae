% Run Tethys-Chloris at the point scale

% T&C scripts can be found on https://hyd.ifu.ethz.ch/research-data-models/t-c.html

% Define site name and simulation number
catchment_name = 'Piramide';
% Tower = 'Tower1';
SIMnm =   [catchment_name '_Calibration_k_' num2str(kOpt) '_z0_' num2str(z0ms(zcal))];

% % Define path to T&C folder
% path_tc = 'C:\Users\jmelovel\Documents\T&C\'; % Path to git repository, to update

% path_TC = [path_tc 'TC']; % HIMAL master branch
path_TC = [path_tc 'TC_model']; % Mike's branch
% addpath(path_TC) 
cd(path_TC)

% % Add path to functions
% path_func = [path_tc '\TC_setups\Functions\'];
% addpath(path_func)

% Add path to inputs
path_in = [root 'Data\Calibration'];% path_in = [path_tc 'TC_setups\Piramide\RUNS\INPUTS\' SIMnm];
addpath(path_in)

% Add path to store output file
% path_out = [path_tc 'TC_outputs\' catchment_name '\PointScale\' datestr(datetime("today"))];
path_out = [path_tc 'TC_outputs\' SIMnm];

% Add path to save figures
path_fig = [path_out '\Figures\'];

% Create output folder if it doesn't exist
if ~exist(path_out, 'dir')
    mkdir(path_out); addpath(path_out)
    mkdir(path_fig)
end 

% Specify parameters
dt = 3600; % Time step (s)
dth = 1; % Time step (h)
ms = 8; % Number of soil layers 
cc = 1; % Crown area

% Load meteorological and other input data 
load([path_in '\run_data_' catchment_name '_' Tower '_26092023_01-Nov-2022-24-Apr-2023.mat']) % same as 'NANs_AllInter1301_Rswcorr_lr_PiramideT2_09_point.mat'

cld = N;
N = Latm; 

% Subset meteo data if necessary
dateTimes = datetime(Date,'ConvertFrom','datenum');
x1 = find(dateTimes==AWS.DATE(1));
x2 = find(dateTimes==AWS.DATE(end));
NN = x2-x1+1;
Date = Date(x1:x2);
Pr = Pr(x1:x2);
Pre = Pre(x1:x2);
Ta = Ta(x1:x2);
Ws = Ws(x1:x2); 
ea = ea(x1:x2);  
SAD1 = SAD1(x1:x2);
SAD2 = SAD2(x1:x2); 
SAB1 = SAB1(x1:x2);
SAB2 = SAB2(x1:x2); 
N = N(x1:x2); 
cld = cld(x1:x2); 
Tdew = Tdew(x1:x2);
esat = esat(x1:x2);
PARB = PARB(x1:x2);
PARD = PARD(x1:x2);

% Further meteo processing
Ds = esat-ea; % Vapor pressure deficit (Pa)
Ds(Ds < 0) = 0;
Ws(Ws <= 0) = 0.01;

% Specify integration intervals for solar variables if not already loaded
% with meteo data
% t_bef = 0.75; 
% t_aft = 0.25;

% Load atmospheric CO2 data
load([path_in '\Ca_Data.mat']);
Ca_dateTime = datetime(Date_CO2,'ConvertFrom','datenum');
dateTime = datetime(Date,'ConvertFrom','datenum');
Ca = interp1(Ca_dateTime,Ca,dateTime,'nearest','extrap');
clear d1 d2 Date_CO2 Ca_dateTime
Oa = 210000; % Intercellular partial pressure of oxygen (umolO2/mol)

% Format date variable
[YE,MO,DA,HO,~,~] = datevec(Date);
Datam(:,1) = YE; Datam(:,2)= MO; Datam(:,3)= DA; Datam(:,4)= HO;
clear YE MO DA HO MI SE

% Filename for model parameters
PARAM_IC = [path_tc 'TC_setups\MOD_PARAM_calibration.m'];

% Run model
MAIN_FRAME;

% Save model outputs
save([path_out '\' SIMnm '.mat']);

