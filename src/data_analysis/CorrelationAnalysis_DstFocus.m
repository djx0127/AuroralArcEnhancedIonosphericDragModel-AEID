% =========================================================================
% Correlation Analysis under Strong Dst Disturbance Conditions
%
% Description:
% This script analyzes the correlation between tangential acceleration
% (ionospheric drag proxy) and geomagnetic activity (Dst index) 
%
% The analysis is based on multi-satellite constellation data (A/B/C),
% including:
%   - TLE-derived orbital parameters
%   - IRI-coupled ionospheric parameters
%   - Geomagnetic indices (Dst)
% =========================================================================

clear;
close all;

%% ------------------------------------------------------------------------
% 1. TLE DATA LOADING AND ORGANIZATION
% ------------------------------------------------------------------------

% List of TLE dataset files (MATLAB format)
TLE_data_files = {
    'data/TLE_data_collection/YG-9_MatlabData_collection.mat', 'YG-9';
    'data/TLE_data_collection/YG-16_MatlabData_collection.mat', 'YG-16';
    'data/TLE_data_collection/YG-17_MatlabData_collection.mat', 'YG-17';
    'data/TLE_data_collection/YG-20_MatlabData_collection.mat', 'YG-20';
    'data/TLE_data_collection/YG-25_MatlabData_collection.mat', 'YG-25';
    'data/TLE_data_collection/YG-31-2_MatlabData_collection.mat', 'YG-31-2';
    'data/TLE_data_collection/YG-31-3_MatlabData_collection.mat', 'YG-31-3';
    'data/TLE_data_collection/YG-31-4_MatlabData_collection.mat', 'YG-31-4';
};

% Preallocate storage: each row = one constellation, columns = A/B/C satellites
all_TLE_data = cell(length(TLE_data_files), 3);
legend_labels = cell(length(TLE_data_files), 3);

% Structures for merged time and force sequences
time_str = struct();
ft_str   = struct();

for i = 1:length(TLE_data_files)
    filename = TLE_data_files{i,1};
    satellite_name = TLE_data_files{i,2};
    
    % Load A/B/C satellite data
    temp_A = load(filename, 'data_A');
    temp_B = load(filename, 'data_B');
    temp_C = load(filename, 'data_C');
    
    all_TLE_data{i,1} = temp_A.data_A;
    all_TLE_data{i,2} = temp_B.data_B;
    all_TLE_data{i,3} = temp_C.data_C;
    
    % Merge time series (Julian day)
    time_str.(['timeA' num2str(i)]) = all_TLE_data{i,1}.Time_jday;
    time_str.(['timeB' num2str(i)]) = all_TLE_data{i,2}.Time_jday;
    time_str.(['timeC' num2str(i)]) = all_TLE_data{i,3}.Time_jday;
    
    % Merge tangential acceleration (ft)
    ft_str.(['ftA' num2str(i)]) = all_TLE_data{i,1}.ft;
    ft_str.(['ftB' num2str(i)]) = all_TLE_data{i,2}.ft;
    ft_str.(['ftC' num2str(i)]) = all_TLE_data{i,3}.ft;
end

% Load geomagnetic Dst data
load('data/TLE_data_collection/DstData_2011to2024.mat');

%% ------------------------------------------------------------------------
% 2. ORBITAL AND FORCE VISUALIZATION
% ------------------------------------------------------------------------

figure
% --- Semi-major axis ---
subplot(2,1,1)
hold on; grid on; box on;
colors = lines(size(all_TLE_data, 1));
lg = zeros(1,size(all_TLE_data, 1));

for i = 1:size(all_TLE_data, 1)
    for j = 1:3
        data_struct = all_TLE_data{i,j};
        lg(i) = plot(data_struct.Time_jday, data_struct.sma,'.', ...
            'MarkerSize', 8, 'Color', colors(i,:));
    end
end

legend(lg,{'YG-9','YG-16','YG-17','YG-20','YG-25','YG-31-2','YG-31-3','YG-31-4'});
xlabel('Time [day]');
ylabel('\ita \rm[km]');
xlim([all_TLE_data{1,1}.Time_jday(1),all_TLE_data{5,1}.Time_jday(end)])
set(gca,'FontSize',14,'FontName','Times New Roman');

% --- Tangential acceleration ---
subplot(2,1,2)
hold on; grid on; box on;

for i = 1:size(all_TLE_data, 1)
    for j = 1:3
        data_struct = all_TLE_data{i,j};
        lg(i) = semilogy(data_struct.Time_jday, data_struct.ft,'.', ...
            'MarkerSize', 8, 'Color', colors(i,:));
    end
end

legend(lg,{'YG-9','YG-16','YG-17','YG-20','YG-25','YG-31-2','YG-31-3','YG-31-4'});
xlabel('Time [day]');
ylabel('\itf_t \rm[N/kg]');
xlim([all_TLE_data{1,1}.Time_jday(1),all_TLE_data{5,1}.Time_jday(end)])
set(gca,'FontSize',14,'FontName','Times New Roman');

%% ------------------------------------------------------------------------
% 3. IRI DATA LOADING
% ------------------------------------------------------------------------

IRI_data_YG9  = IRI_Data_Load('YG-9');
IRI_data_YG16 = IRI_Data_Load('YG-16');
IRI_data_YG17 = IRI_Data_Load('YG-17');
IRI_data_YG20 = IRI_Data_Load('YG-20');
IRI_data_YG25 = IRI_Data_Load('YG-25');
IRI_data_YG312 = IRI_Data_Load('YG-31-2');
IRI_data_YG313 = IRI_Data_Load('YG-31-3');
IRI_data_YG314 = IRI_Data_Load('YG-31-4');

IRI_data_all = {IRI_data_YG9, IRI_data_YG16, IRI_data_YG17, ...
                IRI_data_YG20, IRI_data_YG25, ...
                IRI_data_YG312, IRI_data_YG313, IRI_data_YG314};

%% ------------------------------------------------------------------------
% 4. SPACE ENVIRONMENT PARAMETERS VISUALIZATION
% ------------------------------------------------------------------------

figure
colors = lines(8);
lg = zeros(1,8);

% --- ap index ---
subplot(2,1,1)
hold on; grid on; box on;

for i = 1:8
    data = IRI_data_all{i};
    lg(i) = plot(data.satA.time_jday, data.satA.ap,'.','Color',colors(i,:));
    plot(data.satB.time_jday, data.satB.ap,'.','Color',colors(i,:));
    plot(data.satC.time_jday, data.satC.ap,'.','Color',colors(i,:));
end

xlabel('Time'); ylabel('ap index');

% --- F10.7 index ---
subplot(2,1,2)
hold on; grid on; box on;

for i = 1:8
    data = IRI_data_all{i};
    lg(i) = plot(data.satA.time_jday, data.satA.F107,'.','Color',colors(i,:));
    plot(data.satB.time_jday, data.satB.F107,'.','Color',colors(i,:));
    plot(data.satC.time_jday, data.satC.F107,'.','Color',colors(i,:));
end

xlabel('Time'); ylabel('F10.7 (sfu)');
ylim([50,300]);

%% ------------------------------------------------------------------------
% 5. ION DENSITY VISUALIZATION
% ------------------------------------------------------------------------

figure
hold on; grid on; box on;

for i = 1:8
    data = IRI_data_all{i};
    lg(i) = semilogy(data.satA.time_jday, data.satA.ni_Op,'.','Color',colors(i,:));
    semilogy(data.satB.time_jday, data.satB.ni_Op,'.','Color',colors(i,:));
    semilogy(data.satC.time_jday, data.satC.ni_Op,'.','Color',colors(i,:));
end

xlabel('Time');
ylabel('n_i [cm^{-3}]');

%% ------------------------------------------------------------------------
% 6. TIME-ALIGNED ANALYSIS (Dst vs Drag)
% ------------------------------------------------------------------------

% Re-rank time and force sequences
[time_rank, ft_rank] = TimeNForce_Reranking(time_str, ft_str);

time_flag_start = time_rank(1);
time_flag_stop  = time_rank(end);

time_mask = (time_rank > time_flag_start) & (time_rank < time_flag_stop);
time_new = time_rank(time_mask);
ft_new   = abs(ft_rank(time_mask));

% Time discretization (12-hour interval)
time_interval = 12/24;
time_nodes = (time_flag_start:time_interval:time_flag_stop)';

len_new = length(time_nodes);
Dst_6h_avg = zeros(len_new,1);
ft_6h_avg  = zeros(len_new,1);
ft_6h_std  = zeros(len_new,1);

time_display = time_nodes - time_rank(1);

for idx = 1:len_new
    time_current = time_nodes(idx);
    window_start = time_current - time_interval;
    window_stop  = time_current;
    
    window_mask = (time_new>=window_start) & (time_new<window_stop);
    
    ft_in_window = ft_new(window_mask);
    
    if ~isempty(ft_in_window)
        ft_6h_avg(idx) = mean(ft_in_window);
        ft_6h_std(idx) = std(ft_in_window);
        
        if ft_6h_avg(idx) > 1e-6
            ft_6h_avg(idx) = NaN;
            ft_6h_std(idx) = NaN;
        end
    else
        ft_6h_avg(idx) = NaN;
        ft_6h_std(idx) = NaN;
    end
    
    % Compute averaged Dst
    Dst_6h_avg(idx) = GetEachDstData(allData,time_current,24,12);
end

%% ------------------------------------------------------------------------
% 7. SMOOTHED CORRELATION VISUALIZATION
% ------------------------------------------------------------------------

ft_sm  = smoothdata(ft_6h_avg,'rlowess',20);
Dst_sm = smoothdata(Dst_6h_avg,'rlowess',20);

figure

yyaxis left
semilogy(time_display, ft_6h_avg,'b-'); hold on;
semilogy(time_display, ft_sm,'r-');
ylabel('\itf_t \rm[N/kg]');

yyaxis right
plot(time_display, -Dst_6h_avg,'k-'); hold on;
plot(time_display, -Dst_sm,'g-');
ylabel('-Dst [nT]');

xlabel('Time');

%% ------------------------------------------------------------------------
% 8. FUNCTIONS
% ------------------------------------------------------------------------

function Dst = GetEachDstData(allData,time_day,start_Dt,stop_Dt)
% Compute averaged Dst over a time window prior to a given epoch

Dt = start_Dt - stop_Dt;
time_nodes = linspace(time_day-start_Dt/24, ...
                      time_day-(stop_Dt-1)/24, Dt);

len_hour = length(time_nodes);
temp = zeros(1,len_hour);

for jdx = 1:len_hour
    [yr,mn,dy,hr,~,~] = invjday(time_nodes(jdx));
    temp(jdx) = GetDstValue(allData,yr,mn,dy,hr);
end

Dst = mean(temp);

end

function data_str = IRI_Data_Load(MissionName)
% Load IRI dataset for a given satellite constellation

folder = ['CrossScaleData/',MissionName];

fileA = fullfile(folder, [MissionName,'A_pythonData_collection.mat']);
fileB = fullfile(folder, [MissionName,'B_pythonData_collection.mat']);
fileC = fullfile(folder, [MissionName,'C_pythonData_collection.mat']);

dataA = load(fileA);
dataB = load(fileB);
dataC = load(fileC);

data_str = struct();
data_str.satA = dataA;
data_str.satB = dataB;
data_str.satC = dataC;

end