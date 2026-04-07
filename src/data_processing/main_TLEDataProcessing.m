% =========================================================================
% TLE-Based Orbital Decay Reconstruction and Data Preparation
%
% Description:
% This script processes TLE data for a selected satellite constellation,
% removes maneuver phases, computes orbital decay and equivalent tangential
% acceleration, filters outliers, converts orbital elements to LLA, and
% generates structured datasets for further regression and modeling.
% =========================================================================

clear;
close all;

%% ------------------------------------------------------------------------
% 1. CONSTANTS
% ------------------------------------------------------------------------

Re = 6371;        % Earth radius [km]
mu = 398600e9;    % Gravitational parameter [m^3/s^2]

%% ------------------------------------------------------------------------
% 2. SATELLITE SELECTION
% ------------------------------------------------------------------------

SatNum = 314;

switch SatNum
    case 9
        SatString = 'YG-9';
        start_date = datetime('2011-01-01');
        stop_date  = datetime('2017-07-10');
    case 16
        SatString = 'YG-16';
        start_date = datetime('2013-03-01');
        stop_date  = datetime('2019-12-01');
    case 17
        SatString = 'YG-17';
        start_date = datetime('2013-10-01');
        stop_date  = datetime('2022-01-01');
    case 20
        SatString = 'YG-20';
        start_date = datetime('2015-01-01');
        stop_date  = datetime('2022-01-01');
    case 25
        SatString = 'YG-25';
        start_date = datetime('2015-03-01');
        stop_date  = datetime('2024-01-01');
    case 312
        SatString = 'YG-31-2';
        start_date = datetime('2021-03-01');
        stop_date  = datetime('2024-12-01');
    case 313 
        SatString = 'YG-31-3';
        start_date = datetime('2021-05-01');
        stop_date  = datetime('2024-12-01');
    case 314
        SatString = 'YG-31-4';
        start_date = datetime('2021-07-01');
        stop_date  = datetime('2024-12-01');
end

%% ------------------------------------------------------------------------
% 3. TLE DATA LOADING AND ORBITAL ELEMENT EXTRACTION
% ------------------------------------------------------------------------

start_jday = juliandate(start_date);
stop_jday  = juliandate(stop_date);

filepath_A = sprintf('data/TLE_raw/%s 01A.txt', SatString);
filepath_B = sprintf('data/TLE_raw/%s 01B.txt', SatString);
filepath_C = sprintf('data/TLE_raw/%s 01C.txt', SatString);

[time_A,julian_time_A,sma_A,ecc_A,inc_A,OMG_A,aop_A,M_A,r_A] = ...
    TLEDataHandle(filepath_A,start_date,stop_date);

[time_B,julian_time_B,sma_B,ecc_B,inc_B,OMG_B,aop_B,M_B,r_B] = ...
    TLEDataHandle(filepath_B,start_date,stop_date);

[time_C,julian_time_C,sma_C,ecc_C,inc_C,OMG_C,aop_C,M_C,r_C] = ...
    TLEDataHandle(filepath_C,start_date,stop_date);

%% ------------------------------------------------------------------------
% 4. MANEUVER PHASE REMOVAL
% ------------------------------------------------------------------------

mutation_time = DeletePhaseMutations(SatNum);

delete_day = 5;
delete_period_jday = zeros(length(mutation_time),2);

for idx = 1:length(mutation_time)
    delete_period_jday(idx,1) = mutation_time(idx) - delete_day;
    delete_period_jday(idx,2) = mutation_time(idx) + delete_day;
end

[Time_day_A,time1_A,sma1_A,ecc1_A,inc1_A,OMG1_A,aop1_A,M1_A,r1_A] = ...
    TLEelments_Mutation(delete_period_jday,time_A,julian_time_A,sma_A,ecc_A,inc_A,OMG_A,aop_A,M_A,r_A);

[Time_day_B,time1_B,sma1_B,ecc1_B,inc1_B,OMG1_B,aop1_B,M1_B,r1_B] = ...
    TLEelments_Mutation(delete_period_jday,time_B,julian_time_B,sma_B,ecc_B,inc_B,OMG_B,aop_B,M_B,r_B);

[Time_day_C,time1_C,sma1_C,ecc1_C,inc1_C,OMG1_C,aop1_C,M1_C,r1_C] = ...
    TLEelments_Mutation(delete_period_jday,time_C,julian_time_C,sma_C,ecc_C,inc_C,OMG_C,aop_C,M_C,r_C);

%% ------------------------------------------------------------------------
% 5. ORBITAL DECAY AND TANGENTIAL ACCELERATION
% ------------------------------------------------------------------------

[Time_dadt_A,dadt_A,Delta_A] = DadtGenerate(sma1_A,Time_day_A);
[Time_dadt_B,dadt_B,Delta_B] = DadtGenerate(sma1_B,Time_day_B);
[Time_dadt_C,dadt_C,Delta_C] = DadtGenerate(sma1_C,Time_day_C);

% Velocity magnitude and true anomaly
for idx = 1:length(sma1_A)
    [~,~,v_norm_A(idx),tht1_A(idx)] = ...
        element2eci(sma1_A(idx),ecc1_A(idx),inc1_A(idx),OMG1_A(idx),aop1_A(idx),M1_A(idx),r1_A(idx));
end

for idx = 1:length(sma1_B)
    [~,~,v_norm_B(idx),tht1_B(idx)] = ...
        element2eci(sma1_B(idx),ecc1_B(idx),inc1_B(idx),OMG1_B(idx),aop1_B(idx),M1_B(idx),r1_B(idx));
end

for idx = 1:length(sma1_C)
    [~,~,v_norm_C(idx),tht1_C(idx)] = ...
        element2eci(sma1_C(idx),ecc1_C(idx),inc1_C(idx),OMG1_C(idx),aop1_C(idx),M1_C(idx),r1_C(idx));
end

% Tangential acceleration (elliptical orbit formulation)
f_TLE_A = -dadt_A .* mu ./ (2 .* (sma1_A*1e3).^2 .* v_norm_A);
f_TLE_B = -dadt_B .* mu ./ (2 .* (sma1_B*1e3).^2 .* v_norm_B);
f_TLE_C = -dadt_C .* mu ./ (2 .* (sma1_C*1e3).^2 .* v_norm_C);

[Time_dadt_A,f_TLE_A,sma1_A,ecc1_A,inc1_A,OMG1_A,aop1_A,M1_A,r1_A,tht1_A,Delta_A] = ...
    DeleteNan(Time_dadt_A,f_TLE_A,sma1_A,ecc1_A,inc1_A,OMG1_A,aop1_A,M1_A,r1_A,tht1_A,Delta_A);

[Time_dadt_B,f_TLE_B,sma1_B,ecc1_B,inc1_B,OMG1_B,aop1_B,M1_B,r1_B,tht1_B,Delta_B] = ...
    DeleteNan(Time_dadt_B,f_TLE_B,sma1_B,ecc1_B,inc1_B,OMG1_B,aop1_B,M1_B,r1_B,tht1_B,Delta_B);

[Time_dadt_C,f_TLE_C,sma1_C,ecc1_C,inc1_C,OMG1_C,aop1_C,M1_C,r1_C,tht1_C,Delta_C] = ...
    DeleteNan(Time_dadt_C,f_TLE_C,sma1_C,ecc1_C,inc1_C,OMG1_C,aop1_C,M1_C,r1_C,tht1_C,Delta_C);

%% ------------------------------------------------------------------------
% 6. OUTLIER REMOVAL
% ------------------------------------------------------------------------

[Time_dadt_A_clean,Delta_A_clean,f_TLE_A_clean,sma_A_clean,ecc_A_clean,inc_A_clean,OMG_A_clean,aop_A_clean,M_A_clean,r_A_clean,tht_A_clean] = ...
    RemoveOutliersBydSMA(Time_dadt_A,f_TLE_A,sma1_A,ecc1_A,inc1_A,OMG1_A,aop1_A,M1_A,r1_A,tht1_A,Delta_A);

[Time_dadt_B_clean,Delta_B_clean,f_TLE_B_clean,sma_B_clean,ecc_B_clean,inc_B_clean,OMG_B_clean,aop_B_clean,M_B_clean,r_B_clean,tht_B_clean] = ...
    RemoveOutliersBydSMA(Time_dadt_B,f_TLE_B,sma1_B,ecc1_B,inc1_B,OMG1_B,aop1_B,M1_B,r1_B,tht1_B,Delta_B);

[Time_dadt_C_clean,Delta_C_clean,f_TLE_C_clean,sma_C_clean,ecc_C_clean,inc_C_clean,OMG_C_clean,aop_C_clean,M_C_clean,r_C_clean,tht_C_clean] = ...
    RemoveOutliersBydSMA(Time_dadt_C,f_TLE_C,sma1_C,ecc1_C,inc1_C,OMG1_C,aop1_C,M1_C,r1_C,tht1_C,Delta_C);

%% ------------------------------------------------------------------------
% 7. ORBITAL ELEMENT → LLA CONVERSION
% ------------------------------------------------------------------------

Time_jday_A = julian_time_A(1) + Time_dadt_A_clean;
Time_jday_B = julian_time_B(1) + Time_dadt_B_clean;
Time_jday_C = julian_time_C(1) + Time_dadt_C_clean;

[Lat_A,Lon_A,Alt_A] = LLAGenerate(sma_A_clean,ecc_A_clean,inc_A_clean,OMG_A_clean,aop_A_clean,M_A_clean,r_A_clean,Time_jday_A);
[Lat_B,Lon_B,Alt_B] = LLAGenerate(sma_B_clean,ecc_B_clean,inc_B_clean,OMG_B_clean,aop_B_clean,M_B_clean,r_B_clean,Time_jday_B);
[Lat_C,Lon_C,Alt_C] = LLAGenerate(sma_C_clean,ecc_C_clean,inc_C_clean,OMG_C_clean,aop_C_clean,M_C_clean,r_C_clean,Time_jday_C);

%% ------------------------------------------------------------------------
% 8. FINAL DATA STRUCTURE GENERATION
% ------------------------------------------------------------------------

Data_A = AllDataCollection_ElemCalc_Final(Time_jday_A,sma_A_clean,ecc_A_clean,inc_A_clean,OMG_A_clean,aop_A_clean,M_A_clean,Lat_A,Lon_A,Alt_A,Delta_A_clean,f_TLE_A_clean);
Data_B = AllDataCollection_ElemCalc_Final(Time_jday_B,sma_B_clean,ecc_B_clean,inc_B_clean,OMG_B_clean,aop_B_clean,M_B_clean,Lat_B,Lon_B,Alt_B,Delta_B_clean,f_TLE_B_clean);
Data_C = AllDataCollection_ElemCalc_Final(Time_jday_C,sma_C_clean,ecc_C_clean,inc_C_clean,OMG_C_clean,aop_C_clean,M_C_clean,Lat_C,Lon_C,Alt_C,Delta_C_clean,f_TLE_C_clean);

%% ------------------------------------------------------------------------
% 9. VISUALIZATION
% ------------------------------------------------------------------------

figure

subplot(2,1,1)
plot(Data_A.t_jday,Data_A.sma,'.'); hold on; grid on;
plot(Data_B.t_jday,Data_B.sma,'.');
plot(Data_C.t_jday,Data_C.sma,'.');
xlabel('Time [day]'); ylabel('\ita \rm[km]');

subplot(2,1,2)
plot(Data_A.t_jday,-Data_A.Delt_a,'.'); hold on; grid on;
plot(Data_B.t_jday,-Data_B.Delt_a,'.');
plot(Data_C.t_jday,-Data_C.Delt_a,'.');
xlabel('Time [day]'); ylabel('\Delta\ita \rm[m]');

%% save data

save('data/TLE_IRI_dataset/YG31-4_RegressionData.mat','Data_A','Data_B','Data_C');

%% ========================================================================
% SUBFUNCTIONS
% ========================================================================

function [Lat,Lon,Alt] = LLAGenerate(sma,ecc,inc,OMG,aop,M,r,jd)
% Convert orbital elements to latitude, longitude, altitude
Re = 6371;
Len = length(sma);

Lat = zeros(Len,1);
Lon = zeros(Len,1);

for idx = 1:Len
    [r_eci,v_eci,~,~] = element2eci(sma(idx),ecc(idx),inc(idx),OMG(idx),aop(idx),M(idx),r(idx));
    [r_ecef,~] = eci2ecef_dj(r_eci,v_eci,jd(idx));
    [Lat(idx),Lon(idx),~] = ecef2geod(r_ecef(1),r_ecef(2),r_ecef(3));
end

Alt = r - Re;
end

function [Time_dadt,dadt,Delt_a] = DadtGenerate(sma1,Time_day)
% Compute semi-major axis derivative and variation
len = length(sma1);
dadt = zeros(len-1,1);
Delt_a = zeros(len-1,1);

for idx = 1:len-1
    if (Time_day(idx+1)-Time_day(idx)) >= 5
        dadt(idx) = NaN;
        Delt_a(idx) = NaN;
    else
        dadt(idx) = (sma1(idx+1)-sma1(idx))/(Time_day(idx+1)-Time_day(idx))/86400*1e3;
        Delt_a(idx) = (sma1(idx+1)-sma1(idx))*1e3;
    end
end

dadt(len) = dadt(end);
Delt_a(len) = Delt_a(end);

dadt(dadt>0) = NaN;
Delt_a(Delt_a>0) = NaN;

Time_dadt = Time_day;
end

function [Time1,f1,sma1,ecc1,inc1,OMG1,aop1,M1,r1,tht1,dadt1] = DeleteNan(Time,f,sma,ecc,inc,OMG,aop,M,r,tht,dadt)
% Remove NaN values
mask = ~isnan(f);

Time1 = Time(mask);
f1 = f(mask);
sma1 = sma(mask);
ecc1 = ecc(mask);
inc1 = inc(mask);
OMG1 = OMG(mask);
aop1 = aop(mask);
M1 = M(mask);
r1 = r(mask);
tht1 = tht(mask);
dadt1 = dadt(mask);
end

function [Time1,dadt1,f1,sma1,ecc1,inc1,OMG1,aop1,M1,r1,tht1] = RemoveOutliersBydSMA(Time,f,sma,ecc,inc,OMG,aop,M,r,tht,dadt)
% Remove outliers based on tangential acceleration threshold
threshold = 1e-6;
mask = ~(f > threshold);

Time1 = Time(mask);
f1 = f(mask);
sma1 = sma(mask);
ecc1 = ecc(mask);
inc1 = inc(mask);
OMG1 = OMG(mask);
aop1 = aop(mask);
M1 = M(mask);
r1 = r(mask);
tht1 = tht(mask);
dadt1 = dadt(mask);
end