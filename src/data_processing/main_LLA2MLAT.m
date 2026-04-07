%% ===============================================================
%  LLA-to-MLAT/MLT Batch Conversion Script
%  ---------------------------------------------------------------
%  Description:
%    - Load original data files containing structures (Data_A/B/C)
%    - Convert geographic coordinates (Lat/Lon/Alt) to
%      magnetic latitude (MLAT) and magnetic local time (MLT)
%    - Append MLAT/MLT to the structures
%    - Save to new .mat files
%  Dependencies:
%    - invjday:    Convert Julian day to calendar date
%    - aacgm_v2_convert: from AACGM-v2
%    - magneticLocalTime: from AACGM-v2
%  ===============================================================

clear; clc;

%% ---------------------------------------------------------------
% 1. Define input and output data file paths
% ---------------------------------------------------------------
TLE_data_files = {
    'RegressionData2026/YG9_IRIData.mat';
    'RegressionData2026/YG16_IRIData.mat';
    'RegressionData2026/YG17_IRIData.mat';
    'RegressionData2026/YG20_IRIData.mat';
    'RegressionData2026/YG25_IRIData.mat';
    'RegressionData2026/YG31-2_IRIData.mat';
    'RegressionData2026/YG31-3_IRIData.mat';
    'RegressionData2026/YG31-4_IRIData.mat';
};

TLE_data_files_New = {
    'RegressionData2026/YG9_IRIData1.mat';
    'RegressionData2026/YG16_IRIData1.mat';
    'RegressionData2026/YG17_IRIData1.mat';
    'RegressionData2026/YG20_IRIData1.mat';
    'RegressionData2026/YG25_IRIData1.mat';
    'RegressionData2026/YG31-2_IRIData1.mat';
    'RegressionData2026/YG31-3_IRIData1.mat';
    'RegressionData2026/YG31-4_IRIData1.mat';
};

% Names of structures to be processed in each .mat file
structNames = {'Data_A', 'Data_B', 'Data_C'};

%% ---------------------------------------------------------------
% 2. Initialize parallel pool
% ---------------------------------------------------------------
pool = gcp('nocreate');           % Get current pool, but do not create a new one
if isempty(pool)
    parpool(8);                   % Start a pool with 8 workers if none exists
end

%% ---------------------------------------------------------------
% 3. Loop over all files and process each structure
% ---------------------------------------------------------------
for i = 1:length(TLE_data_files)
    filePath = TLE_data_files{i};
    savePath = TLE_data_files_New{i};
    load(filePath);               % Load Data_A, Data_B, Data_C, etc.
    
    % Loop over each structure name
    for j = 1:length(structNames)
        sName = structNames{j};
        tempStruct = eval(sName);
        
        % Extract original data
        % (Assume Lat_seris, Lon_seris, Alt_seris, time_seris are
        %  arrays of the same size)
        lat  = tempStruct.Lat_seris;
        lon  = tempStruct.Lon_seris;
        alt  = tempStruct.Alt_seris * 1e-3;   % Convert altitude from m to km
        jday = tempStruct.time_seris;         % Julian day
        num_day   = size(jday, 1);
        num_seris = size(jday, 2);
        
        MLAT = zeros(num_day, num_seris);
        MLT  = zeros(num_day, num_seris);
        
        % -------------------------------------------------------
        % 3.1 Parallel loop over time rows
        % -------------------------------------------------------
        parfor idx = 1:num_day
            fprintf('Computing file %d, struct %s, row %d\n', i, sName, idx);
            for jdx = 1:num_seris
                % Convert Julian day to calendar date
                [year, mon, day, ~, ~, ~] = invjday(jday(idx, jdx));
                matlab_date = datetime(year, mon, day);
                
                % Convert to AACGM magnetic latitude and longitude
                [MLAT(idx, jdx), c_lon, ~] = aacgm_v2_convert( ...
                    lat(idx, jdx), ...
                    lon(idx, jdx), ...
                    alt(idx, jdx), ...
                    matlab_date, ...
                    0, 0);
                
                % Compute magnetic local time (MLT)
                [MLT(idx, jdx), ~, ~] = magneticLocalTime(matlab_date, c_lon);
            end
        end
        
        % Append new fields to the structure
        tempStruct.MLAT = MLAT;
        tempStruct.MLT  = MLT;
        
        % Write updated structure back to workspace variable
        assignin('base', sName, tempStruct);
    end
    
    %% -----------------------------------------------------------
    % 4. Save updated structures and clean up
    % -----------------------------------------------------------
    save(savePath, 'Data_A', 'Data_B', 'Data_C');
    clear('Data_A', 'Data_B', 'Data_C');
end

fprintf('\nAll files have been processed.\n');

%% ---------------------------------------------------------------
% 5. Shut down parallel pool
% ---------------------------------------------------------------
delete(gcp('nocreate'));