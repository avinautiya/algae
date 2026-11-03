%% Calculate z0 with Ts and including d

% This script calculates the roughness length (z0) using the temperature
% structure parameter method, incorporating debris temperature data. It
% reads data from various tower levels, fits empirical models to the data,
% and filters the results based on specific conditions.

%% Read file
cd(path_data)

if tower == 1
    clear Results Param stats Stats WIND TRH DT
end

dataname = [SITE '_d'];

% Load data
load([path_data '\DATA_' SITE '.mat']) % structure wih structures:
% 1) wind: timetables of towers data different elevations with: L1_FF (windspeed, m/s), L1_FFmx (max windspeed, m/s), L1_FFstd (windspeed std, wm/s), L1_Pulse, L2_FF, L2_FFmx, L2_FFstd, L2_Pulse, BATT_1, L3_FF, L3_FFmx, L3_FFstd, DIR_L3 (wind direction, °), L4_FF, L4_FFmx, L4_FFstd, DIR_L4, BATT_2
% 2) trh: timetables of towers with data at different elevations: TA_1 (air temperature), RH_1 (relative humidity), DTA_1 (dew temperature), TA_2, RH_2, DTA_2, TA_3, RH_3, DTA_3, TA_4, RH_4, DTA_4
% 3) dt: timetables of towers with data at different depths: DT1, DT2, DT3, DT4, DT5

% Choose tower
Tower = char(horzcat('Tower',num2str(tower)));

% Period of Tower with fewer measurements
date_len = min([length(wind.Tower1.DATE),length(wind.Tower2.DATE),length(wind.Tower3.DATE)]);
idx = find([length(wind.Tower1.DATE),length(wind.Tower2.DATE),length(wind.Tower3.DATE)]==date_len);
Tower_min = ['Tower' char(num2str(idx))];
date_ini = wind.(Tower_min).DATE(1); date_end = wind.(Tower_min).DATE(end);

% Get wind, temperature and DT data
WIND.(Tower).data = wind.(Tower); % wind.(Tower) is a time table that contains the wind speed the max, the std and the pulse for all the levels of (Tower) and the wind direction for levels 3 and 4
TRH.(Tower).data = trh.(Tower); % trh.(Tower) is a time table that contains the air temperature, the relative humidity and dew temperature for all levels
DT.(Tower).data = debtemp.(Tower); % debtemp.(Tower) is a time table that contains the debristemperature for all levels

% Retime WIND and TRH to match DT time steps 
endDate = WIND.(Tower).data.DATE(end);
endDate = endDate-minutes(minute(endDate));
newTimes = DT.(Tower).data.DATE(1):minutes(30):endDate;
WIND.(Tower).data = retime(WIND.(Tower).data,newTimes,'linear');
TRH.(Tower).data = retime(TRH.(Tower).data,newTimes,'linear');
DT.(Tower).data = retime(DT.(Tower).data,newTimes);

% Extract time step 
time = minutes(WIND.(Tower).data.DATE(2)-WIND.(Tower).data.DATE(1));

if t3 == 1
    % Cut to period of Tower with fewer measurements
    ST3 = timerange(date_ini,date_end+minutes(time));
    WIND.(Tower).data = WIND.(Tower).data(ST3,:);
    TRH.(Tower).data = TRH.(Tower).data(ST3,:);
    DT.(Tower).data = DT.(Tower).data(ST3,:);
end

cd(path_codes)

%% General variables

% Define variables
alpha_m = 5; % Sicart et al. 2014, Miles et al. 2017 -> 0.5
alpha_h = alpha_m;
Pr = 1; % Sicart et al. 2014 and Miles et al., 2017
k = 0.40; % Sicart et al. 2014 and Miles et al. 2017, k = 0.41 Rounce et al. 2015
g = 9.81; % m s^-2

% Indexes for the wind and temperature data
if surfT == 0
    L_FF = ["L2_FF","L3_FF","L4_FF"];
    TA = ["TA_2","TA_3","TA_4"];
    % Sensor elevations and uncertainties
    z = [1, 2, 2.7];
    sigu = [0.45, 0.89, 0.89]; % [m/s] 0.45 m/s for the lower ones and 2mph (0.89 m/s) for the upper ones
    sigT = 0.35; % [°C] between 0 and 50 °C
elseif surfT == 1
    L_FF = ["L1_FF","L2_FF","L3_FF","L4_FF"];
    TA = ["TA_1","TA_2","TA_3","TA_4"];
    % Sensor elevations and uncertainties
    z = [0.5, 1, 2, 2.7];
    sigu = [0.45, 0.45, 0.89, 0.89]; % [m/s] 0.45 m/s for the lower ones and 2mph (0.89 m/s) for the upper ones
    sigT = 0.35; % [°C] between 0 and 50 °C
end

% Remove level 3 for Kyzlysu because temperature did not work
if strcmp(SITE,'Kyzylsu')
    if tower == 3
        L_FF = ["L1_FF","L2_FF","L4_FF"];
        TA = ["TA_1","TA_2","TA_4"];
        z = [0.5, 1, 2.7];
        sigu = [0.45, 0.45, 0.89];
    end
end

z4L = erase(string(z),".");
L_names = strcat("L", z4L); L_names = convertStringsToChars(L_names);

wind_var = string(WIND.(Tower).data.Properties.VariableNames);
trh_var = string(TRH.(Tower).data.Properties.VariableNames);
w_ind = find(ismember(wind_var,L_FF));
t_ind = find(ismember(trh_var,TA));

%% Stability with Ts from debris thermistor

if surfT == 0
    % Convert to K for Rib calculation
    for i = t_ind
    TRH.(Tower).data.(i) = TRH.(Tower).data.(i)+273.15;
    end
    TRH.(Tower).data.TA_1 = TRH.(Tower).data.TA_1+273.15;
    
    len = length(TRH.(Tower).data.DATE);
    % Richardson number: Rib = g*(Tz-Ts)/(T0*u^2)*(z-zs)
    for i = 1:len
        Tz = [];
        u = [];
        for jj = 1:length(t_ind)
            Tz = [Tz TRH.(Tower).data.(t_ind(jj))(i)]; % Temperature at height z [K]
            u = [u WIND.(Tower).data.(w_ind(jj))(i)]; % wind vel at height z 
        end
        Ts = TRH.(Tower).data.TA_1(i); % Temperature at the surface [K]
        T0 = mean(Tz); % Mean temperature of the air layer [K]
        zs = 0.5;
        TRH.(Tower).data.Rib(i) = mean(g.*(Tz-Ts)./(T0.*u.^2).*(z-zs)); 
    end
elseif surfT == 1
        % Convert to K for Rib calculation
    for i = t_ind
    TRH.(Tower).data.(i) = TRH.(Tower).data.(i)+273.15;
    end
    DT.(Tower).data.DT1 = DT.(Tower).data.DT1+273.15;
    
    len = length(TRH.(Tower).data.DATE);
    % Richardson number: Rib = g*(Tz-Ts)/(T0*u^2)*(z-zs)
    for i = 1:len
        Tz = [];
        u = [];
        for jj = 1:length(t_ind)
            Tz = [Tz TRH.(Tower).data.(t_ind(jj))(i)]; % Temperature at height z [K]
            u = [u WIND.(Tower).data.(w_ind(jj))(i)]; % wind vel at height z 
        end
        Ts = DT.(Tower).data.DT1(i); % Temperature at the surface [K]
        T0 = mean(Tz); % Mean temperature of the air layer [K]
        zs = 0;
        TRH.(Tower).data.Rib(i) = mean(g.*(Tz-Ts)./(T0.*u.^2).*(z-zs)); 
    end
end

% Filter data with -0.03<Rib<0.03 - Miles et al. 2017
idx_rib = TRH.(Tower).data.Rib<0.03 & TRH.(Tower).data.Rib>-0.03;
TRH.(Tower).Rib = TRH.(Tower).data(idx_rib,:);
WIND.(Tower).Rib = WIND.(Tower).data(idx_rib,:);
DT.(Tower).Rib = DT.(Tower).data(idx_rib,:);

% Filter data with u2m>=1.5 - Miles et al. 2017
idx_u2m = WIND.(Tower).data.L3_FF>=1.5;
TRH.(Tower).u2m = TRH.(Tower).data(idx_u2m,:);
WIND.(Tower).u2m = WIND.(Tower).data(idx_u2m,:);
DT.(Tower).u2m = DT.(Tower).data(idx_u2m,:);

% Filter both - Miles et al. 2017
idx_ribu2m = TRH.(Tower).data.Rib<0.03 & TRH.(Tower).data.Rib>-0.03 & WIND.(Tower).data.L3_FF>=1.5;
TRH.(Tower).Ribu2m = TRH.(Tower).data(idx_ribu2m,:);
WIND.(Tower).Ribu2m = WIND.(Tower).data(idx_ribu2m,:);
DT.(Tower).Ribu2m = DT.(Tower).data(idx_ribu2m,:);

%% z0 fmincon with costwapper for uncertainties

warning('off', 'MATLAB:table:RowsAddedExistingVars')

% Set solver options
maxIter = 1000;
stepSize = 0.001;
tolX = 0.001;
options = optimoptions('fmincon','MaxIter',maxIter,...
    'UseParallel',false,'FiniteDifferenceStepSize',stepSize,'Display',...
    'off','StepTolerance',tolX,'MaxFunctionEvaluations',inf);

% Define number of iterations 
iter = 10; % 

fn = fieldnames(WIND.(Tower));
f1 = waitbar(0, 'Starting');
if ~exist('Param.(fn{jj})','var')
    Results.(Tower).(fn{1}) = table();
    Param.(Tower).(fn{1}) = table();
end
len = length(WIND.(Tower).data.DATE);
% CM = parula(len);
for i = 1:len % Iterate thorugh time
    waitbar(i/len, f1, sprintf('Calculating fits: %d %%', floor(i/len*100)));
    t = [];
    u = [];
    for II = 1:length(t_ind) % Creates a vector with the temperatures and one for the wind speeds for the current timestep
    t = [t TRH.(Tower).data.(t_ind(II))(i)]; % Temperature at height z [K]
    u = [u WIND.(Tower).data.(w_ind(II))(i)]; % wind vel at height z 
    end
    if surfT == 0
        Ts = TRH.(Tower).data.TA_1(i); % Temperature at the surface [K]
    elseif surfT == 1
        Ts = DT.(Tower).data.DT1(i);
    end
    
    if dval == 0
        % Initial conditions and bounds
        x0 = [1, 0.01, 1, 0.01];
        lb =  [0 0 -1 0];
        ub = [10 0.1 10 0.1];
        L = 10^8.*ones(1,length(w_ind));
    elseif dval == 1
        % Initial conditions and bounds
        x0 = [1, 0.01, 1, 0.01, 0.01];
        lb =  [0 0 -1 0 -0.1];
        ub = [10 0.1 10 0.1 0.1];
        L = 10^8.*ones(1,length(w_ind));
    end

    for ii = 1:iter % Iterate 10 times and then stop
        Results.(Tower).data.Iteration(ii) = ii;
        for l = 1:length(L_names)
            Results.(Tower).data.(L_names{l})(ii) = L(l);
        end

        % Minimise F
            F = @(x) costwrapper_all(x,z,u,t,sigu,sigT,k,Pr,alpha_m,alpha_h,L,Ts,dval,uncert);
            problem = createOptimProblem("fmincon","objective",F,"x0",...
                x0,"lb",lb,"ub",ub,'options',options);
            [x,fval,exitflag,output] = fmincon(problem);
        
        if dval == 0
            
            % Fitted values
            uf_fit = x(1);
            z0_fit = x(2);
            tf_fit = x(3);
            zt_fit = x(4);
            
            % Store values
            j = costwrapper_all(x,z,u,t,sigu,sigT,k,Pr,alpha_m,alpha_h,L,Ts,dval,uncert);
            Results.(Tower).data.uf(ii) = uf_fit;
            Results.(Tower).data.z0(ii) = z0_fit;
            Results.(Tower).data.tf(ii) = tf_fit;
            Results.(Tower).data.zt(ii) = zt_fit;
            Results.(Tower).data.j(ii) = j;
            
            % Calculate new L values
            L = (t.*uf_fit^2)/(k*g*tf_fit);
        
            % Define new initial conditions for the next iteration
            x0 = [uf_fit z0_fit tf_fit zt_fit];
        
            if ii<iter
                for l = 1:length(L_names)
                Results.data.(L_names{l})(ii+1) = L(l);
                end
            end
             
            % Fitted values
            u_fitted = (uf_fit/k).*(log((z)./z0_fit)+alpha_m.*((z)./L));
            t_fitted = ((tf_fit./k).*(Pr.*log((z)./zt_fit)+alpha_h.*((z)./L)))+Ts;

        elseif dval == 1
        
            % Fitted values
            uf_fit = x(1);
            z0_fit = x(2);
            tf_fit = x(3);
            zt_fit = x(4);
            d_fit = x(5);
            
            % Store values
            j = costwrapper_all(x,z,u,t,sigu,sigT,k,Pr,alpha_m,alpha_h,L,Ts,dval,uncert);
            Results.(Tower).data.uf(ii) = uf_fit;
            Results.(Tower).data.z0(ii) = z0_fit;
            Results.(Tower).data.tf(ii) = tf_fit;
            Results.(Tower).data.zt(ii) = zt_fit;
            Results.(Tower).data.d(ii) = d_fit;
            Results.(Tower).data.j(ii) = j;
            
            % Calculate new L values
            L = (t.*uf_fit^2)/(k*g*tf_fit);
        
            % Define new initial conditions for the next iteration
            x0 = [uf_fit z0_fit tf_fit zt_fit d_fit];
        
            if ii<iter
                for l = 1:length(L_names)
                Results.data.(L_names{l})(ii+1) = L(l);
                end
            end
             
            % Fitted values
            u_fitted = (uf_fit/k).*(log((z+d_fit)./z0_fit)+alpha_m.*((z+d_fit)./L));
            t_fitted = ((tf_fit./k).*(Pr.*log((z+d_fit)./zt_fit)+alpha_h.*((z+d_fit)./L)))+Ts;
        end

        % Calculate r2
        [r2_u, rmse_u] = rsquare(u,u_fitted);
        [r2_t, rmse_t] = rsquare(t,t_fitted);
    end

    % Store parameters of the last iteration
    Param.(Tower).data.DATE(i)=WIND.(Tower).data.DATE(i);
    for l = 1:length(L_names)
        Param.(Tower).data.(L_names{l})(i) = L(l);
    end

    Param.(Tower).data.uf(i)=uf_fit;
    Param.(Tower).data.z0(i)=z0_fit;
    Param.(Tower).data.tf(i)=tf_fit;
    Param.(Tower).data.zt(i)=zt_fit;
    if dval == 1
        Param.(Tower).data.d(i)=d_fit;
    end
    Param.(Tower).data.j(i) = j;
    Param.(Tower).data.r2_u(i) = r2_u;
    Param.(Tower).data.r2_t(i) = r2_t;
    Param.(Tower).data.rmse_u(i) = rmse_u;
    Param.(Tower).data.rmse_t(i) = rmse_t;

    Ldiff = [];
    for l = 1:length(L_names)
        Ldiff = [Ldiff,Results.(Tower).data.(L_names{l})(ii)-Results.(Tower).data.(L_names{l})(ii-1)];
    end
    Ldiff = mean(Ldiff);
    Param.(Tower).data.Converges(i) = Ldiff;
end
close(f1)

%% Filter Rib, u2m and Ribu2m
Param.(Tower).Rib = Param.(Tower).data(idx_rib,:);
Param.(Tower).u2m = Param.(Tower).data(idx_u2m,:);
Param.(Tower).Ribu2m = Param.(Tower).data(idx_ribu2m,:);

sets = length(fieldnames(Param.(Tower)));

for jj = 1:length(fn)
    % Calculate z0 and zt stats
    mz0 = mean(Param.(Tower).(fn{jj}).z0);
    gmz0 = geomean(Param.(Tower).(fn{jj}).z0);
    stdz0 = std(Param.(Tower).(fn{jj}).z0);
    mez0 = median(Param.(Tower).(fn{jj}).z0);
    MADz0 = mad(Param.(Tower).(fn{jj}).z0);
    c=-1/(sqrt(2)*erfcinv(3/2)); 
    NMADz0 = c*MADz0;
    minz0 = min(Param.(Tower).(fn{jj}).z0);
    maxz0 = max(Param.(Tower).(fn{jj}).z0);
    mzt = mean(Param.(Tower).(fn{jj}).zt);
    stdzt = std(Param.(Tower).(fn{jj}).zt);
    mezt = median(Param.(Tower).(fn{jj}).zt);
    MADzt = mad(Param.(Tower).(fn{jj}).zt);
    NMADzt = c*MADzt;
    minzt = min(Param.(Tower).(fn{jj}).zt);
    maxzt = max(Param.(Tower).(fn{jj}).zt);
    if dval == 1
        md = mean(Param.(Tower).(fn{jj}).d);
        med = median(Param.(Tower).(fn{jj}).d);
        MADd = mad(Param.(Tower).(fn{jj}).d);
        NMADd = c*MADd;
        mind = min(Param.(Tower).(fn{jj}).d);
        maxd = max(Param.(Tower).(fn{jj}).d);
        stats.(Tower).(fn{jj})(1,:) = [mz0,gmz0,stdz0,mez0,MADz0,NMADz0,minz0,maxz0,...
        mzt,stdzt,mezt,MADzt,NMADzt,minzt,maxzt,...
        md,med,MADd,NMADd,mind,maxd];
    elseif dval == 0
        stats.(Tower).(fn{jj})(1,:) = [mz0,gmz0,stdz0,mez0,MADz0,NMADz0,minz0,maxz0,mzt,stdzt,mezt,MADzt,NMADzt,minzt,maxzt];
    end
end

%% Filter profiles with L condition and j/R^2 condition

f1 = waitbar(0, 'Starting');
for jj = 1:sets
    waitbar(jj/sets, f1, sprintf('Dataset filtering: %d', jj));
    filter1 = strcat(fn{jj},'_j');
    filter2 = strcat(fn{jj},'_r2');
    % % Filter out results with j>1 and that do not converge (Ldiff=/=0))
    idx1 = ((Param.(Tower).(fn{jj}).("j")<1) & (Param.(Tower).(fn{jj}).("Converges")==0));
    Param.(Tower).(filter1) = Param.(Tower).(fn{jj})(idx1,:);
    WIND.(Tower).(filter1) = WIND.(Tower).(fn{jj})(idx1,:);
    TRH.(Tower).(filter1) = TRH.(Tower).(fn{jj})(idx1,:);
    DT.(Tower).(filter1) = DT.(Tower).(fn{jj})(idx1,:);
    % Filter out results with r2_u<0.75 & r2_t<0.75 and that do not converge (Ldiff=/=0))
    idx2 = ((Param.(Tower).(fn{jj}).("r2_u")>0.75 & Param.(Tower).(fn{jj}).("r2_t")>0.75) & Param.(Tower).(fn{jj}).("Converges")==0);
    Param.(Tower).(filter2) = Param.(Tower).(fn{jj})(idx2,:);
    WIND.(Tower).(filter2) = WIND.(Tower).(fn{jj})(idx2,:);
    TRH.(Tower).(filter2) = TRH.(Tower).(fn{jj})(idx2,:); 
    DT.(Tower).(filter2) = DT.(Tower).(fn{jj})(idx2,:);
    Param.(Tower).(fn{jj}).idx1 = idx1;
    Param.(Tower).(fn{jj}).idx2 = idx2;
end
close(f1)

% Plot measurements and estimations of the profiles that meet the criteria
fn_f = fieldnames(Param.(Tower));
fn_f = fn_f(length(sets)+1:end);
ti = ["j<1","r2>0.75","j<1","r2>0.75","j<1","r2>0.75","j<1","r2>0.75"];

% Plot measurements and calculate stats of the profiles that meet the criteria
f2 = waitbar(0, 'Starting');
for kk = 4:length(fn_f) % 10-11 = Ribu2m_j and Rib2um_r2
    if isempty(Param.(Tower).(fn_f{kk}))==1
        stats.(Tower).(fn_f{kk}) = nan(size(stats.(Tower).(fn{jj})));
        continue
    end
    if strcmp(fn_f{kk},'Ribu2m_j')
        if plt == 1
        fig=figure;
        set(gcf,'position',[2000,1000,1000,500])
        end
        c_line = parula(12);
        c = phasemap(24);
        shift = 13;
        CM = [c(shift+1:end,:);c(1:shift,:)];
%         CM = vertcat(c_line,flipud(c_line));
        for i = 1:length(Param.(Tower).(fn_f{kk}).DATE)
            waitbar(i/length(Param.(Tower).(fn_f{kk}).DATE), f2, sprintf('Plotting: %d %%', floor(i/length(Param.(Tower).(fn_f{kk}).DATE)*100)));
            if hour(Param.(Tower).(fn_f{kk}).DATE(i))==0
                col = CM(24,:);
            else 
                col = CM(hour(Param.(Tower).(fn_f{kk}).DATE(i)),:);
            end
            t = [];
            u = [];
            for II = 1:length(t_ind)
            t = [t TRH.(Tower).(fn_f{kk}).(t_ind(II))(i)]; % Temperature at height z [K]
            u = [u WIND.(Tower).(fn_f{kk}).(w_ind(II))(i)]; % wind vel at height z 
            end
            Ts = DT.(Tower).(fn_f{kk}).DT1(i);
            uf_fit = Param.(Tower).(fn_f{kk}).uf(i);
            z0_fit = Param.(Tower).(fn_f{kk}).z0(i);
            tf_fit = Param.(Tower).(fn_f{kk}).tf(i);
            zt_fit = Param.(Tower).(fn_f{kk}).zt(i);
            L = [];
            for l = 1:length(L_names)
                L = [L, Param.(Tower).(fn_f{kk}).(L_names{l})(i)];
            end
            if dval == 0
                u_fitted = (uf_fit/k).*(log((z)./z0_fit)+alpha_m.*((z)./L));
                t_fitted = ((tf_fit./k).*(Pr.*log((z)./zt_fit)+alpha_h.*((z)./L)))+Ts;
            elseif dval == 1
                u_fitted = (uf_fit/k).*(log((z+d_fit)./z0_fit)+alpha_m.*((z+d_fit)./L));
                t_fitted = ((tf_fit./k).*(Pr.*log((z+d_fit)./zt_fit)+alpha_h.*((z+d_fit)./L)))+Ts;
            end
            
            if plt == 1
            subplot(1,2,1)
            hold on; grid on; box on;
            errorbar(log(z), u, sigu,'*', 'color', col)
            plot(log(z), u_fitted, 'color', col)
            title(['Tower ' num2str(tower) ' filtered wind profiles (' num2str(length(Param.(Tower).(fn_f{kk}).j)) '/' num2str(length(WIND.(Tower).data.DATE)) ')'])
            xlabel('ln(z)')
            ylabel('W_s [m/s]')
            legend('Data','Fit')
            subplot(1,2,2)
            hold on; grid on; box on
            errorbar(log(z), t, repelem(sigT,length(z)),'*', 'color', col)
            plot(log(z), t_fitted, 'color', col)
            colormap(CM)
            colbar = colorbar; colbar.Ticks = linspace(0,0.95,12); colbar.TickLabels = 0:2:23; colbar.Label.String = 'Hour';  
            title(['Tower ' num2str(tower) ' filtered temperature profiles (' num2str(length(Param.(Tower).(fn_f{kk}).j)) '/' num2str(length(WIND.(Tower).data.DATE)) ')'])
            xlabel('ln(z)')
            ylabel('T_a [K]')
            legend('Data','Fit')
            end
        end
        
        if plt == 1
        % Save profile plots
        cd(path_figures)
        if exist('ST3','var')
            figname_png = ['Tower_' num2str(tower) '_Filtered_Profiles-' (fn_f{kk}) '_'...
                dataname '_' char(string(time)) '_Tower3period.png'];
            figname_pdf = ['Tower_' num2str(tower) '_Filtered_Profiles-' (fn_f{kk}) '_'...
                dataname '_' char(string(time)) '_Tower3period.pdf'];
        %     saveas(fig,figname)  
            exportgraphics(fig, figname_png, 'Resolution', 300);
            exportgraphics(fig, figname_pdf, 'Resolution', 300);
            close(fig)
        else
            figname_png = ['Tower_' num2str(tower) '_Filtered_Profiles-' (fn_f{kk}) '_'...
            dataname '_' char(string(time)) '.png'];
            figname_pdf = ['Tower_' num2str(tower) '_Filtered_Profiles-' (fn_f{kk}) '_'...
            dataname '_' char(string(time)) '.pdf'];
        %     saveas(fig,figname)  
            exportgraphics(fig, figname_png, 'Resolution', 300);
            exportgraphics(fig, figname_pdf, 'Resolution', 300);
            close(fig)
        end
        end
    end

    % Calculate new z0 and zt stats (filtered)
    avgz0_f = mean(Param.(Tower).(fn_f{kk}).z0);
    gmeanz0_f = geomean(Param.(Tower).(fn_f{kk}).z0);
    stdz0_f = std(Param.(Tower).(fn_f{kk}).z0);
    medz0_f = median(round(Param.(Tower).(fn_f{kk}).z0,3));
    MADz0_f = mad(Param.(Tower).(fn_f{kk}).z0);
    c=-1/(sqrt(2)*erfcinv(3/2)); 
    NMADz0_f = c*MADz0_f;
    minz0_f = min(Param.(Tower).(fn_f{kk}).z0);
    maxz0_f = max(Param.(Tower).(fn_f{kk}).z0);
    avgzt_f = mean(Param.(Tower).(fn_f{kk}).zt);
    stdzt_f = std(Param.(Tower).(fn_f{kk}).zt);
    medzt_f = median(round(Param.(Tower).(fn_f{kk}).zt,3));
    MADzt_f = mad(Param.(Tower).(fn_f{kk}).zt);
    NMADzt_f = c*MADzt_f;
    minzt_f = min(Param.(Tower).(fn_f{kk}).zt);
    maxzt_f = max(Param.(Tower).(fn_f{kk}).zt);
    if dval == 1 
        avgd_f = mean(Param.(Tower).(fn_f{kk}).d);
        medd_f = median(round(Param.(Tower).(fn_f{kk}).d,3));
        MADd_f = mad(Param.(Tower).(fn_f{kk}).d);
        NMADd_f = c*MADd_f;
        mind_f = min(Param.(Tower).(fn_f{kk}).d);
        maxd_f = max(Param.(Tower).(fn_f{kk}).d);
    
        stats.(Tower).(fn_f{kk}) = [avgz0_f,gmeanz0_f,stdz0_f,medz0_f,MADz0_f,NMADz0_f,minz0_f,maxz0_f,...
            avgzt_f,stdzt_f,medzt_f,MADzt_f,NMADzt_f,minzt_f,maxzt_f,...
            avgd_f,medd_f,MADd_f,NMADd_f,mind_f,maxd_f];
    elseif dval == 0
        stats.(Tower).(fn_f{kk}) = [avgz0_f,gmeanz0_f,stdz0_f,medz0_f,MADz0_f,NMADz0_f,minz0_f,maxz0_f,...
            avgzt_f,stdzt_f,medzt_f,MADzt_f,NMADzt_f,minzt_f,maxzt_f];
    end
    
end
% close(f1)
close(f2)

fn_s = fieldnames(stats.(Tower));
if dval == 0
    Names = ["Subset","z0 mean","z0 Geomean","z0 std","z0 median","z0 MAD","z0 NMAD","z0 min",...
        "z0 max","zt mean","zt std","zt median","zt MAD","zt NMAD","zt min","zt max"];
    varTypes = ["string","double","double","double","double","double","double","double",...
        "double","double","double","double","double","double","double","double"];
elseif dval == 1
    Names = ["Subset","z0 mean","z0 Geomean","z0 std","z0 median","z0 MAD","z0 NMAD","z0 min",...
        "z0 max","zt mean","zt std","zt median","zt MAD","zt NMAD","zt min","zt max","d mean",...
        "d median","d MAD","d NMAD","d min","d max"];
    varTypes = ["string","double","double","double","double","double","double","double",...
        "double","double","double","double","double","double","double","double",...
        "double","double","double","double","double","double"];
end
sz = [length(fn_s) length(Names)];
Stats.(Tower) =  table('Size',sz,'VariableNames',...
    Names,'VariableTypes',varTypes);
Stats.(Tower).Subset = string(fn_s);
for ll = 1:length(fn_s)
    Stats.(Tower)(ll,2:end) = num2cell(stats.(Tower).(fn_s{ll}));
end


%% Save
cd(path_results)

save(['Parameters_z0_' SITE '_' char(string(time)) '_' d_text '_' sT_text '_' unc_text '_' t3_text '.mat'],"Param");
save(['Statistics_z0_' SITE '_' char(string(time)) '_' d_text '_' sT_text '_' unc_text '_' t3_text '.mat'],"Stats");





