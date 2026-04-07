%% ===============================================================
%  Auroral Arc Charging Parameter Optimization
%  ---------------------------------------------------------------
%  Case: YG-31-2,3,4 – Third-Level Optimization
%  ===============================================================
clear;
close all;

%% ---------------------------------------------------------------
% 1. Area-to-mass ratio normalization parameters
% ---------------------------------------------------------------
%% Area-to-mass ratio normalization parameters
% YG-9 to YG-25
A1 = 1.4*1.2+2*1.4*2;
m1 = 640;
% YG-31
A2 = 1.5*1.3 + 2*3.5*1.5;
m21 = 780; % YG-31 A-satellite
m22 = 720; % YG-31 B/C-satellites
% Three types of area-to-mass ratios and scaling factors
SvM_1  = A1/m1;
SvM_21 = A2/m21; c_SvM21 = SvM_1/SvM_21;
SvM_22 = A2/m22; c_SvM22 = SvM_1/SvM_22;

%% ---------------------------------------------------------------
% 2. Load TLE/IRI data (managed with cell arrays)
% ---------------------------------------------------------------
load('DstData\DstData_2011to2024.mat');
TLE_data_files = {
    'RegressionData2026/YG31-2_IRIData1.mat', 'YG-31-2';
    'RegressionData2026/YG31-3_IRIData1.mat', 'YG-31-3';
    'RegressionData2026/YG31-4_IRIData1.mat', 'YG-31-4';
};
% Use cell arrays to store processed data for all satellites
num_formations = size(TLE_data_files, 1);
num_sats_per_formation = 3;
total_sats = num_formations * num_sats_per_formation;

processed_data  = cell(num_formations, 1); % Processed data for optimization
satellite_info  = cell(num_formations, 1); % Satellite labels for plotting
all_TLE_data    = cell(length(TLE_data_files), 3); 

for i = 1:length(TLE_data_files)
    filename  = TLE_data_files{i,1};
    base_name = TLE_data_files{i,2};
    
    % Load A/B/C satellite data
    load(filename, 'Data_A', 'Data_B', 'Data_C');
    dadt_rawA = Data_A.Delt_a;
    dadt_rawB = Data_B.Delt_a;
    dadt_rawC = Data_C.Delt_a;
    
    % Normalize Delta a by area-to-mass factors
    Data_A.Delt_a = dadt_rawA * c_SvM21; % A-satellite uses c_SvM21
    Data_B.Delt_a = dadt_rawB * c_SvM22; % B-satellite uses c_SvM22
    Data_C.Delt_a = dadt_rawC * c_SvM22; % C-satellite uses c_SvM22
    
    % Use A-satellite data for global optimization
    raw_data = Data_A;
    [elem,track] = IRIDataLoading(raw_data,allData);
    c_0 = [0,1,1,1,0];
    [y_init,y_star] = ModelOptimizationComputation_3(c_0,elem,track);
    
    processed_data{i} = struct('elem', elem, 'track', track, 'y_star', y_star);
    satellite_info{i} = sprintf('%s-%s', base_name, 'A');
    
    % ========== Prepare plotting data for A/B/C of this formation ==========
    all_TLE_data{i,1} = Data_A;
    all_TLE_data{i,2} = Data_B;
    all_TLE_data{i,3} = Data_C;
    [elem_A,track_A] = IRIDataLoading(all_TLE_data{i,1},allData);
    [elem_B,track_B] = IRIDataLoading(all_TLE_data{i,2},allData);
    [elem_C,track_C] = IRIDataLoading(all_TLE_data{i,3},allData);
    
    c_0 = [0,1,1,1,0];
    [all_TLE_data{i,1}.y_init_A,all_TLE_data{i,1}.y_star_A] ...
        = ModelOptimizationComputation_3(c_0,elem_A,track_A);
    [y_init_B,y_star_B] ...
        = ModelOptimizationComputation_3(c_0,elem_B,track_B);
    [y_init_C,y_star_C] ...
        = ModelOptimizationComputation_3(c_0,elem_C,track_C);
    
    % Merge, sort and compute mean & bounds for visualization
    [all_TLE_data{i,1}.time_all, all_TLE_data{i,1}.y_all, ...
     all_TLE_data{i,1}.upper_bound, all_TLE_data{i,1}.lower_bound]  ...
        = mergeAndSortTimeSeries( ...
            track_A.time_jday, all_TLE_data{i,1}.y_star_A, ...
            track_B.time_jday, y_star_B, ...
            track_C.time_jday, y_star_C);
    
    all_TLE_data{i,1}.y_all_0nan = all_TLE_data{i,1}.y_all;
    all_TLE_data{i,1}.y_all_0nan(isnan(all_TLE_data{i,1}.y_all_0nan)) = [];
    
    % Initial evaluation metrics
    all_TLE_data{i,1}.Dst_1  = dtw(all_TLE_data{i,1}.y_init_A, all_TLE_data{i,1}.y_all_0nan);
    all_TLE_data{i,1}.RMSE_1 = sqrt(mean((all_TLE_data{i,1}.y_init_A - all_TLE_data{i,1}.y_star_A).^2));
end
fprintf('All data of %d satellites have been processed.\n\n', total_sats);
%% ---------------------------------------------------------------
% 3. Global parameter optimization by ga
% ---------------------------------------------------------------
% Parameter bounds
lb = [  1, 0.01,  50, 0.001,  1];  % Lower bounds
rb = [400,    1, 200,   0.5, 30];  % Upper bounds
opts = optimoptions('ga', 'Display', 'iter', ...       
    'MaxGenerations', 100, 'PopulationSize', 200, 'UseParallel', true);     
num_vars = 5; % Number of parameters to be optimized
[optimal_params, fval] = ga(@(p) combinedObjectiveFunction(p, processed_data), ...
                            num_vars, [], [], [], [], lb, rb, [], opts);
disp(optimal_params);
%% ---------------------------------------------------------------
% 4. Re-evaluate model with optimized parameters & visualization
% ---------------------------------------------------------------
for i = 1:length(TLE_data_files)
    elem_A  = processed_data{i}.elem;
    track_A = processed_data{i}.track;

    % Model output with globally optimized parameters
    [all_TLE_data{i,1}.y_model_A,~] = ...
        ModelOptimizationComputation_3(optimal_params,elem_A,track_A);
    
    % DTW evaluation after optimization
    all_TLE_data{i,1}.Dst_A = dtw(all_TLE_data{i,1}.y_model_A, ...
                                  all_TLE_data{i,1}.y_all_0nan);
    % RMSE evaluation after optimization
    all_TLE_data{i,1}.RMSE_A = sqrt(mean((all_TLE_data{i,1}.y_model_A - ...
                                          all_TLE_data{i,1}.y_star_A).^2));
   
    % Key parameters for plotting
    time_start = all_TLE_data{i,1}.t_jday(1);
    time_stop  = all_TLE_data{i,1}.t_jday(end);
    plot_color = [0, 0.4470, 0.7410]; % MATLAB default blue
    fill_color = plot_color;
   
    figure(i)
    set(gcf, 'Position', [100, 100, 630, 450]);
    
    % 4.1 Decay-rate comparison (initial vs optimized vs merged)
    subplot(2,1,1)
    h_act = semilogy(track_A.time_jday-time_start, ...
                     all_TLE_data{i,1}.y_init_A,'k.','MarkerSize',6);
    hold on; grid on; box on;
    h_opt = semilogy(track_A.time_jday-time_start, ...
                     all_TLE_data{i,1}.y_model_A,'b.','MarkerSize',6);
    
    fill_x = [all_TLE_data{i,1}.time_all; flipud(all_TLE_data{i,1}.time_all)];
    fill_y = [all_TLE_data{i,1}.upper_bound; flipud(all_TLE_data{i,1}.lower_bound)];
    nan_indices = isnan(fill_x) | isnan(fill_y);
    fill_x(nan_indices) = [];
    fill_y(nan_indices) = [];
    h_fill = fill(fill_x - time_start, fill_y, fill_color);
    set(h_fill, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    
    h_TLE = semilogy(all_TLE_data{i,1}.time_all-time_start, ...
                     all_TLE_data{i,1}.y_all,'-', ...
                     'Color', plot_color, 'LineWidth', 2.0);
    xlim([0,time_stop-time_start]); yticks([1e-2,1e-1,1]);
    xlabel('Time [day]'); ylabel('da/dt [km/day]');
    legend([h_TLE,h_act,h_opt], ...
        {'Actual Data','Theoretical Data','Data of AEID Model'}, ...
        'NumColumns',3,'Location', 'northoutside');
    set(gca,'FontSize',14); set(gca,'FontName','Times New Roman');
    
    % 4.2 Residuals (log-space) + density
    res_model = log10(all_TLE_data{i,1}.y_star_A) - ...
                log10(all_TLE_data{i,1}.y_model_A);
    res_orign = log10(all_TLE_data{i,1}.y_star_A) - ...
                log10(all_TLE_data{i,1}.y_init_A);
    
    subplot(2,1,2)
    yyaxis left
    h_resm = plot(track_A.time_jday-time_start,res_model,'b-','LineWidth',1.5);
    hold on; grid on; box on;
    h_reso = plot(track_A.time_jday-time_start,res_orign,'k-','LineWidth',1.5);
    xlim([0,time_stop-time_start]); ylim([-1,2]);    
    xlabel('Time [day]'); ylabel('RE (log_{10})');
    set(gca,'FontSize',14); set(gca,'FontName','Times New Roman');
    
    yyaxis right
    plot(track_A.time_jday-time_start,track_A.ni,'LineWidth',1);
    hold on; grid on; box on;
    xlim([0,time_stop-time_start]); ylim([0,2e5]);
    ylabel('Ni [cm^{-3}]');
    legend([h_reso,h_resm], ...
        {'RE of Theoretical Model','RE of AEID Model'}, ...
        'NumColumns',2,'Location', 'northoutside');
    set(gca,'FontSize',14); set(gca,'FontName','Times New Roman');
end

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
        time_centers(k) = (time_bins(k) + time_bins(k+1)) / 2; % 使用箱子中心作为时间点
    end
end

end

function total_error = combinedObjectiveFunction(c_params, Data_All)
% object function for the model optimization
    total_error = 0;
    for i = 1:length(Data_All)
        elem   = Data_All{i}.elem;
        track  = Data_All{i}.track;
        y_star = Data_All{i}.y_star;
        [y_model,~] = ModelOptimizationComputation_3(c_params,elem,track);
        temp_error = sqrt(mean((y_model-y_star).^2));
        total_error = total_error + temp_error;
    end
end
