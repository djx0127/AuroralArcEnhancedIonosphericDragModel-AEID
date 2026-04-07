% =========================================================================
% Averaged Orbital Decay Analysis and Correlation with F10.7
%
% Description:
% This script averages the orbital decay (tangential acceleration) of
% three satellites within each constellation and analyzes its correlation
% with solar activity (F10.7 index).
% =========================================================================

clear;
close all;

%% ------------------------------------------------------------------------
% 1. SURFACE-TO-MASS RATIO NORMALIZATION
% ------------------------------------------------------------------------

% YG9–YG25 satellites
A1 = 1.4*1.2 + 2*1.4*2;
m1 = 640;

% YG31 satellites
A2 = 1.5*1.3 + 2*3.5*1.5;
m21 = 780; % Satellite A
m22 = 720; % Satellites B/C

% Surface-to-mass ratios and scaling factors
SvM_1  = A1/m1;
SvM_21 = A2/m21; c_SvM21 = SvM_1/SvM_21;
SvM_22 = A2/m22; c_SvM22 = SvM_1/SvM_22;

%% ------------------------------------------------------------------------
% 2. DATA LOADING AND NORMALIZATION
% ------------------------------------------------------------------------

% Load Dst data
load('data/TLEdata_collection/DstData_2011to2024.mat');

% Dataset list
TLE_data_files = {
    'data/TLE_IRI_dataset/YG9_IRIData1.mat', 'YG-9';
    'data/TLE_IRI_dataset/YG16_IRIData1.mat', 'YG-16';
    'data/TLE_IRI_dataset/YG17_IRIData1.mat', 'YG-17';
    'data/TLE_IRI_dataset/YG20_IRIData1.mat', 'YG-20';
    'data/TLE_IRI_dataset/YG25_IRIData1.mat', 'YG-25';
    'data/TLE_IRI_dataset/YG31-2_IRIData1.mat', 'YG-31-2';
    'data/TLE_IRI_dataset/YG31-3_IRIData1.mat', 'YG-31-3';
    'data/TLE_IRI_dataset/YG31-4_IRIData1.mat', 'YG-31-4';
};

% Preallocation
all_TLE_data = cell(length(TLE_data_files), 3);
time_str = struct();
ft_str   = struct();

for i = 1:length(TLE_data_files)
    
    filename = TLE_data_files{i,1};
    satellite_name = TLE_data_files{i,2};
    
    % Load A/B/C data
    load(filename, 'Data_A');
    load(filename, 'Data_B');
    load(filename, 'Data_C');
    
    all_TLE_data{i,1} = Data_A;
    all_TLE_data{i,2} = Data_B;
    all_TLE_data{i,3} = Data_C;
    
    % Raw tangential acceleration
    ft_A_raw = Data_A.ft;
    ft_B_raw = Data_B.ft;
    ft_C_raw = Data_C.ft;
    
    % Normalize YG31 series
    if contains(satellite_name, 'YG-31')
        ft_A_norm = ft_A_raw * c_SvM21;
        ft_B_norm = ft_B_raw * c_SvM22;
        ft_C_norm = ft_C_raw * c_SvM22;
        fprintf('Normalized ft for %s.\n', satellite_name);
    else
        ft_A_norm = ft_A_raw;
        ft_B_norm = ft_B_raw;
        ft_C_norm = ft_C_raw;
        fprintf('Reference group: %s (no normalization).\n', satellite_name);
    end
    
    % Store time
    time_str.(['timeA' num2str(i)]) = Data_A.t_jday;
    time_str.(['timeB' num2str(i)]) = Data_B.t_jday;
    time_str.(['timeC' num2str(i)]) = Data_C.t_jday;
    
    % Store normalized force
    ft_str.(['ftA' num2str(i)]) = ft_A_norm;
    ft_str.(['ftB' num2str(i)]) = ft_B_norm;
    ft_str.(['ftC' num2str(i)]) = ft_C_norm;
    
    % Save into structure
    all_TLE_data{i,1}.ft_norm = ft_A_norm;
    all_TLE_data{i,2}.ft_norm = ft_B_norm;
    all_TLE_data{i,3}.ft_norm = ft_C_norm;
end

%% ------------------------------------------------------------------------
% 3. BASIC VISUALIZATION (SMA, ft, F107)
% ------------------------------------------------------------------------

tjd_start = min([all_TLE_data{1,1}.t_jday(1), all_TLE_data{1,2}.t_jday(1), all_TLE_data{1,3}.t_jday(1)]);
tjd_stop  = max([all_TLE_data{8,1}.t_jday(end),all_TLE_data{8,2}.t_jday(end),all_TLE_data{8,3}.t_jday(end)]);

figure
% Semi-major axis
subplot(3,1,1); hold on; grid on;
for i = 1:size(all_TLE_data,1)
    for j = 1:3
        plot(all_TLE_data{i,j}.t_jday - tjd_start, ...
             all_TLE_data{i,j}.sma,'.');
    end
end
ylabel('\ita \rm[km]');

% Tangential acceleration
subplot(3,1,2); hold on; grid on;
for i = 1:size(all_TLE_data,1)
    for j = 1:3
        semilogy(all_TLE_data{i,j}.t_jday - tjd_start, ...
                 all_TLE_data{i,j}.ft_norm,'.');
    end
end
ylabel('\itf_t \rm[N/kg]');

% F107
subplot(3,1,3); hold on; grid on;
for i = 1:size(all_TLE_data,1)
    for j = 1:3
        plot(all_TLE_data{i,j}.t_jday - tjd_start, ...
             all_TLE_data{i,j}.F107(:,1),'.');
    end
end
ylabel('F107 [sfu]');
xlabel('Time [day]');

%% ------------------------------------------------------------------------
% 4. ION DENSITY (Ni, Ne)
% ------------------------------------------------------------------------

figure
subplot(2,1,1); hold on; grid on;
for i = 1:size(all_TLE_data,1)
    for j = 1:3
        semilogy(all_TLE_data{i,j}.t_jday - tjd_start, ...
                 all_TLE_data{i,j}.Ni(:,1),'.');
    end
end
ylabel('\itN_i \rm[m^{-3}]');

subplot(2,1,2); hold on; grid on;
for i = 1:size(all_TLE_data,1)
    for j = 1:3
        semilogy(all_TLE_data{i,j}.t_jday - tjd_start, ...
                 all_TLE_data{i,j}.Ne(:,1),'.');
    end
end
ylabel('\itN_e \rm[m^{-3}]');

%% ------------------------------------------------------------------------
% 5. TIME-BINNED STATISTICAL ANALYSIS
% ------------------------------------------------------------------------

% Merge all data
all_times_jd  = [];
all_ft_values = [];

ft_fields = fieldnames(ft_str);

for k = 1:length(ft_fields)
    tf = ft_fields{k};
    tt = strrep(tf,'ft','time');
    all_times_jd  = [all_times_jd; time_str.(tt)];
    all_ft_values = [all_ft_values; ft_str.(tf)];
end

% Remove invalid data
valid = isfinite(all_times_jd) & isfinite(all_ft_values) & (all_ft_values >= 1e-12);
all_times_jd  = all_times_jd(valid);
all_ft_values = all_ft_values(valid);

% Time binning
bin_width_days = 50;
time_bins = floor(min(all_times_jd)) : bin_width_days : ceil(max(all_times_jd));

mean_ft     = nan(length(time_bins)-1,1);
upper_bound = nan(length(time_bins)-1,1);
lower_bound = nan(length(time_bins)-1,1);
time_centers = nan(length(time_bins)-1,1);

for k = 1:length(time_bins)-1
    
    idx = (all_times_jd >= time_bins(k)) & (all_times_jd < time_bins(k+1));
    
    if any(idx)
        log_ft = log10(all_ft_values(idx));
        mu = mean(log_ft);
        sigma = std(log_ft);
        
        mean_ft(k)     = 10^mu;
        upper_bound(k) = 10^(mu + sigma);
        lower_bound(k) = 10^(mu - sigma);
        
        time_centers(k) = (time_bins(k) + time_bins(k+1))/2;
    end
end

%% ------------------------------------------------------------------------
% 6. FINAL VISUALIZATION (ft vs F107)
% ------------------------------------------------------------------------

figure
yyaxis left; hold on; grid on;

fill_x = [time_centers; flipud(time_centers)];
fill_y = [upper_bound; flipud(lower_bound)];

valid = ~(isnan(fill_x) | isnan(fill_y));
fill(fill_x(valid)-tjd_start, fill_y(valid), ...
    [0 0.4470 0.7410],'FaceAlpha',0.3,'EdgeColor','none');

plot(time_centers - tjd_start, mean_ft,'-','LineWidth',2);
set(gca,'YScale','log');
ylabel('\itf_t \rm[m/s^2]');

% F107 overlay
yyaxis right
all_times_F107 = [];
all_F107_values = [];

for i = 1:size(all_TLE_data,1)
    for j = 1:3
        all_times_F107  = [all_times_F107; all_TLE_data{i,j}.t_jday];
        all_F107_values = [all_F107_values; all_TLE_data{i,j}.F107(:,1)];
    end
end

plot(all_times_F107 - tjd_start, all_F107_values,'.');
ylabel('F107 [sfu]');
xlabel('Time [day]');