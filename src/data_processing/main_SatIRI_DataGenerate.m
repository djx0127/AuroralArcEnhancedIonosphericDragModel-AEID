% Associate IRI data with TLE data from main_TLEDataProcessing
% load: 'data/TLE_IRI_dataset/YG9_RegressionData.mat'
% save: 'data/TLE_IRI_dataset/YG9_IRIData.mat'

clear; clc;
close all; 
% Complete the data processing for 9, 16, 17, 20, 25, 31-2, 31-3, 31-4 in sequence.
load('data/TLE_IRI_dataset/YG9_RegressionData.mat');

Data_A = IRI_Generate(Data_A);
Data_B = IRI_Generate(Data_B);
Data_C = IRI_Generate(Data_C);

save('data/TLE_IRI_dataset/YG9_IRIData.mat','Data_A','Data_B','Data_C');

function Data = IRI_Generate(Data)
% handle single satellite
t_jday   = datetime(Data.t_jday, 'convertfrom', 'juliandate');
jd_minor = datetime(Data.time_seris, 'convertfrom', 'juliandate');
Lat = Data.Lat_seris;
Lon = Data.Lon_seris;
Alt = Data.Alt_seris/1e3;

N = length(t_jday);
Kp   = zeros(N,100);
F107 = zeros(N,100); 
Ni   = zeros(N,100);
Ne   = zeros(N,100); 
Ti   = zeros(N,100); 
Te   = zeros(N,100); 
for idx = 1:N
    time_vec = jd_minor(idx,:);
    lat_vec = Lat(idx,:);
    lon_vec = Lon(idx,:);
    alt_vec = Alt(idx,:);
    Out = IRI2020(time_vec, lat_vec, lon_vec, alt_vec, 'sat');
    Kp(idx,:) = Out.Kp';
    F107(idx,:) = Out.F107';
    Ni(idx,:) = Out.O_p';
    Ne(idx,:) = Out.dens';
    Ti(idx,:) = Out.Ti';
    Te(idx,:) = Out.Te';
end

Data.Kp = Kp;
Data.F107 = F107;
Data.Ni = Ni;
Data.Ne = Ne;
Data.Ti = Ti;
Data.Te = Te;
end
