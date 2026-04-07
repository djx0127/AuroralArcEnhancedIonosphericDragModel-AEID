%% ===============================================================
%  Auroral Arc Charging Parameter Optimization
%  ---------------------------------------------------------------
%  Case: YG-9 – First-Level Optimization
%
%  Workflow:
%    1) Define geometric area-to-mass ratios and normalization factors
%    2) Load TLE/IRI data and Dst data
%    3) Extract and preprocess A/B/C satellite data
%    4) Run initial model (with activation function + theoretical model)
%    5) Merge and average multi-satellite time series
%    6) Optimize model parameter using fminsearch
%    7) Re-evaluate model with optimized parameter (DTW & RMSE)
%    8) Visualize solar/geophysical indices and model performance
%  ===============================================================
clear;
close all;

%% ---------------------------------------------------------------
% 1. Area-to-mass ratio normalization parameters
% ---------------------------------------------------------------
% YG-9 / YG-25
A1 = 1.4*1.2 + 2*1.4*2;
m1 = 640;
% YG-31 (2, 3, 4)
A2 = 1.5*1.3 + 2*3.5*1.5;
m21 = 780;   % YG-31 A-satellite mass
m22 = 720;   % YG-31 B/C-satellite mass
% Three types of area-to-mass ratios and scaling factors
SvM_1  = A1 / m1;
SvM_21 = A2 / m21;  c_SvM21 = SvM_1 / SvM_21;
SvM_22 = A2 / m22;  c_SvM22 = SvM_1 / SvM_22;

%% ---------------------------------------------------------------
% 2. Load TLE / IRI data (managed via cell arrays)
% ---------------------------------------------------------------
load('data/TLE_data_collection/DstData_2011to2024.mat');
TLE_data_files = {
    'data/TLE_IRI_dataset/YG9_IRIData1.mat', 'YG-9';
};
% Pre-allocate cell arrays for all formations and their A/B/C satellites
% each row: one formation (file)
% each column: A / B / C satellite data
all_TLE_data  = cell(length(TLE_data_files), 3);
legend_labels = cell(length(TLE_data_files), 3);
idx = 1;
filename        = TLE_data_files{idx,1};
satellite_name  = TLE_data_files{idx,2};
% Load A/B/C satellite data from the same file
load(filename, 'Data_A');
load(filename, 'Data_B');
load(filename, 'Data_C');
all_TLE_data{idx,1} = Data_A;
all_TLE_data{idx,2} = Data_B;
all_TLE_data{idx,3} = Data_C;
idx = 1;   % keep index explicit for later calls

%% ---------------------------------------------------------------
% 3. Extract and organize data for A/B/C satellites
% ---------------------------------------------------------------
[elem_A, track_A] = IRIDataLoading(all_TLE_data{idx,1}, allData);
[elem_B, track_B] = IRIDataLoading(all_TLE_data{idx,2}, allData);
[elem_C, track_C] = IRIDataLoading(all_TLE_data{idx,3}, allData);

%% ---------------------------------------------------------------
% 4. Initial model evaluation (activation function + theoretical model)
% ---------------------------------------------------------------
c_0 = 0;
[yj_init_A, y_star_A] = ModelOptimizationComputation_1(c_0, elem_A, track_A);
[yj_init_B, y_star_B] = ModelOptimizationComputation_1(c_0, elem_B, track_B);
[yj_init_C, y_star_C] = ModelOptimizationComputation_1(c_0, elem_C, track_C);
% Initial evaluation metrics
Dst_1  = dtw(yj_init_A, y_star_A);                  % Dynamic time warping distance
RMSE_1 = sqrt(mean((yj_init_A - y_star_A).^2));      % Root mean square error

%% ---------------------------------------------------------------
% 5. Merge A/B/C series and compute mean & bounds
% ---------------------------------------------------------------
[time_all, y_all, upper_bound, lower_bound] = mergeAndSortTimeSeries( ...
    track_A.time_jday, y_star_A, ...
    track_B.time_jday, y_star_B, ...
    track_C.time_jday, y_star_C);
% Remove NaN entries for later global comparison
y_all_0nan = y_all;
y_all_0nan(isnan(y_all_0nan)) = [];

%% ---------------------------------------------------------------
% 6. Parameter optimization by fminsearch
% ---------------------------------------------------------------
c_initial = 15;
params_A = fminsearch(@(p) objectiveFunction(p, elem_A, track_A, y_all), c_initial);

%% ---------------------------------------------------------------
% 7. Re-run model with optimized parameter & evaluate
% ---------------------------------------------------------------
[y_model_A, ~] = ModelOptimizationComputation_1(params_A, elem_A, track_A);
[y_model_B, ~] = ModelOptimizationComputation_1(params_A, elem_B, track_B);
[y_model_C, ~] = ModelOptimizationComputation_1(params_A, elem_C, track_C);
% DTW-based performance metrics (vs merged series)
Dst_A = dtw(y_model_A, y_all_0nan);
Dst_B = dtw(y_model_B, y_all_0nan);
Dst_C = dtw(y_model_C, y_all_0nan);
% RMSE-based performance metrics (vs each satellite's reference y_star_*)
RMSE_A = sqrt(mean((y_model_A - y_star_A).^2));
RMSE_B = sqrt(mean((y_model_B - y_star_B).^2));
RMSE_C = sqrt(mean((y_model_C - y_star_C).^2));

%% ---------------------------------------------------------------
% 8. Visualization of solar/geophysical indices and density
% ---------------------------------------------------------------
time_start = all_TLE_data{1,1}.t_jday(1);
time_stop  = all_TLE_data{1,3}.t_jday(end);
plot_color = [0, 0.4470, 0.7410];   % MATLAB default blue
fill_color = plot_color;

window_width = 50;
F107_smooth = movmean(track_B.F107, window_width);
Dst_smooth  = movmean(track_B.Dst,  window_width);
ni_smooth   = movmean(track_A.ni,   window_width);
% 8.1 Solar and geomagnetic activity indices
figure(2)
subplot(2,1,1)
hold on; grid on; box on;
plot(track_B.time_jday - time_start, track_B.F107, 'b.', 'LineWidth', 1);
plot(track_B.time_jday - time_start, F107_smooth, '-', 'LineWidth', 1);
xlim([0, time_stop - time_start]);
xlabel('Time'); ylabel('F107');
subplot(2,1,2)
hold on; grid on; box on;
plot(track_B.time_jday - time_start, track_B.Dst, 'b.', 'LineWidth', 1);
plot(track_B.time_jday - time_start, Dst_smooth, '-', 'LineWidth', 1);
xlim([0, time_stop - time_start]);
xlabel('Time'); ylabel('Dst');
% 8.2 Plasma density (ni) in logarithmic scale
figure(3)
semilogy(track_A.time_jday - time_start, track_A.ni, 'b.', 'LineWidth', 1);
hold on; grid on; box on;
semilogy(track_A.time_jday - time_start, ni_smooth, '-', 'LineWidth', 1);
xlim([0, time_stop - time_start]);
xlabel('Time'); ylabel('ni');

%% ---------------------------------------------------------------
% 9. Visualization of Residual analysis (log-scale)
% ---------------------------------------------------------------
Res_log_A = log10(y_star_A) - log10(y_model_A);  
Res_log_ref = log10(y_star_A) - log10(yj_init_A); 
figure(11)
set(gcf, 'Position', [100, 100, 630, 450]);
subplot(2,1,1)
h_act = semilogy(track_A.time_jday-time_start,yj_init_A,'k.','MarkerSize',6);
hold on; grid on; box on;
h_opt = semilogy(track_A.time_jday-time_start,y_model_A,'b.','MarkerSize',6);
fill_x = [time_all; flipud(time_all)];
fill_y = [upper_bound; flipud(lower_bound)];
nan_indices = isnan(fill_x) | isnan(fill_y);
fill_x(nan_indices) = [];
fill_y(nan_indices) = [];
h_fill = fill(fill_x - time_start, fill_y, fill_color);
set(h_fill, 'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
h_TLE = semilogy(time_all-time_start,y_all,'-', 'Color', plot_color, 'LineWidth', 2.0);
xlim([0,time_stop-time_start]); yticks([1e-2,1e-1,1]);
xlabel('Time [day]'); ylabel('da/dt [km/day]');
legend([h_TLE,h_act,h_opt],{'Actual Data','Theoretical Data','Data of AEID Model'},...
    'NumColumns',3,'Location', 'northoutside');
set(gca,'FontSize',14); set(gca,'FontName','Times New Roman');

subplot(2,1,2)
yyaxis left
h_resm = plot(track_A.time_jday-time_start,Res_log_A,'b-','LineWidth',1.5);
hold on; grid on; box on;
h_reso = plot(track_A.time_jday-time_start,Res_log_ref,'k-','LineWidth',1.5);
xlim([0,time_stop-time_start]); ylim([-1,1.5]);
xlabel('Time [day]'); ylabel('RE (log_{10})');
set(gca,'FontSize',14); set(gca,'FontName','Times New Roman');
yyaxis right
plot(track_A.time_jday-time_start,track_A.ni,'LineWidth',1);
hold on; grid on; box on;
xlim([0,time_stop-time_start]); ylim([0,1e5]);
ylabel('Ni [cm^{-3}]');
legend([h_reso,h_resm],{'RE of Theoretical Model','RE of AEID Model'},...
    'NumColumns',2,'Location', 'northoutside');
set(gca,'FontSize',14); set(gca,'FontName','Times New Roman');

%% ========================================================================
% SUBFUNCTIONS
% =========================================================================

function [elem,track] = IRIDataLoading(IRIData,allData)
% Read the data in the struct IRIData 
% Restore in track and elem at TLE-intervals and sub-intervals
track_time_jday = IRIData.t_jday;
Len_time = length(track_time_jday);
track_time_day = track_time_jday-track_time_jday(1);
time_end = track_time_day(end);
track_Dst = zeros(Len_time,1);
for idx = 1:Len_time
    [yr,mn,dy,hr,~,~] = invjday(track_time_jday(idx));
    track_Dst(idx) = GetDstValue(allData,yr,mn,dy,hr);
end

track_ni = IRIData.Ni(:,1)/1e6; % from m-3 to cm-3
track_m = Len_time;

% sub-interval variables i=1~n
elem.Te = IRIData.Te;
elem.Ti = IRIData.Ti;
elem.ne = IRIData.Ne/1e6; % from m-3 to cm-3
elem.ni = IRIData.Ni/1e6; % from m-3 to cm-3
elem.sma = IRIData.a1_seris;
elem.Lat = IRIData.Lat_seris;
elem.time_jday = IRIData.time_seris; 
elem.F107 = IRIData.F107;
elem.MLAT = IRIData.MLAT;
elem.MLT = IRIData.MLT;

% TLE-interval variables j=1~m
track.F107 = IRIData.F107(:,1);
track.Dst  = track_Dst;
track.Kp = IRIData.Kp(:,1);
track.sma  = IRIData.sma; 
track.delt_a = IRIData.Delt_a;
track.time_jday = track_time_jday;
track.time_day = track_time_day;
track.ni   = track_ni;
track.m    = track_m;
track.ft   = IRIData.ft;
end

function [time_centers, y_mean, upper_bound, lower_bound] = mergeAndSortTimeSeries(time_A, y_A, time_B, y_B, time_C, y_C)
% Combine three sets of time series data and sort them by time
combined_time = [time_A(:); time_B(:); time_C(:)];
combined_y    = [y_A(:);    y_B(:);    y_C(:)];
[time_all, sort_order] = sort(combined_time);
y_all = combined_y(sort_order);

bin_width_days = 20;
min_time = floor(min(time_all));
max_time = ceil(max(time_all));
time_bins = min_time : bin_width_days : max_time;
y_mean = nan(length(time_bins) - 1, 1);
upper_bound = nan(length(time_bins) - 1, 1);
lower_bound = nan(length(time_bins) - 1, 1);
time_centers = nan(length(time_bins) - 1, 1);
for k = 1:length(time_bins) - 1
    indices_in_bin = (time_all >= time_bins(k)) & (time_all < time_bins(k+1));
    if any(indices_in_bin)
        y_in_bin = y_all(indices_in_bin);
        log_y = log10(y_in_bin);
        mu = mean(log_y);
        sigma = std(log_y);
        
        y_mean(k) = 10^mu;
        upper_bound(k) = 10^(mu + sigma);
        lower_bound(k) = 10^(mu - sigma);
        time_centers(k) = (time_bins(k) + time_bins(k+1)) / 2;
    end
end
end

function distance = objectiveFunction(c_params, input_elem, input_track, y_all)
% object function for the model optimization
    [y_model,~] = ModelOptimizationComputation_1(c_params,input_elem,input_track);
    y_true = y_all;
    distance = dtw(y_model, y_true);
end
