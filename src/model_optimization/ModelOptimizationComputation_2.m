function [y_j,y_star] = ModelOptimizationComputation_2(c_current,elem,track)
% Compute actual and theoretical semi-major axis decay rates
% Method: mixed large TLE-interval + small sub-interval computation
%
% Input:
%   c_current : current model parameter
%   elem      : small-interval struct, storing Te Ti ne ni sma MLAT MLT time_jday
%   track     : large-interval struct, storing F107 Dst sma delt_a ni time_jday, m
%
% Output:
%   y_j   : actual decay rate
%   y_star: theoretical decay rate

mu = 398601;        % (km3/s2)  

% Small-interval parameters i = 1..n
elem_Te = elem.Te;        % [K]
elem_Ti = elem.Ti;        % [K]
elem_ne = elem.ne;        % [cm-3]
elem_ni = elem.ni;        % [cm-3]
elem_sma = elem.sma;      % [km]
elem_MLAT = elem.MLAT;    % [deg]
elem_MLT  = elem.MLT;     % [h]
elem_time_jday = elem.time_jday; % [jd]

% Large-interval parameters j = 1..m
track_F107 = abs(track.F107); % [sfu]
track_Dst  = abs(track.Dst);  % [nT] Dst has been converted to positive values
track_sma  = track.sma;       % [km]
track_delt_a = track.delt_a;  % [km]
track_time_jday = track.time_jday; % [jd]
track_ni = track.ni;

m = track.m;

%% Smooth Dst and F107
window_width = 50;
track_Dst = movmean(track_Dst, window_width);
track_F107 = movmean(track_F107, window_width);
track_ni = movmean(track_ni, window_width);

%% Compute target function Theta
Veq = zeros(m,100);
fion = zeros(m,100);
elem_delt_a = zeros(m,1);
for idx = 1:m  % Loop over large intervals
    % Equilibrium potential Veq
    Veq(idx,:) = Veq_calc(elem_Te(idx,:),elem_Ti(idx,:),elem_MLAT(idx,:),...
        elem_MLT(idx,:),track_F107(idx),track_ni(idx),c_current);
    % Ion drag acceleration fion
    fion(idx,:) = IonAcc_calc(elem_Te(idx,:),elem_ne(idx,:),elem_ni(idx,:),...
        elem_sma(idx,:),Veq(idx,:));
    % Numerical integration for Delta_a
    ft_temp = fion(idx,:);
    a0_temp = track_sma(idx);
    t_temp = (elem_time_jday(idx,:) - elem_time_jday(idx,1))*86400;
    
    % Integral Int1 for target function
    int_temp1 = 2 .* ft_temp .* sqrt(a0_temp.^3 ./ mu);
    elem_delt_a(idx) = trapz(t_temp, int_temp1);
end

% Unify time step
temp = track_time_jday(2:end);
Delt_t = temp - track_time_jday(1:end-1);
Delt_t(m) = Delt_t(m-1);

% Actual decay rate (from model integration)
elem_dadt = elem_delt_a./Delt_t;
elem_dadt_fill = filloutliers(elem_dadt,'linear','movmedian',20); 
elem_dadt_fill = smoothdata(elem_dadt_fill,'rlowess',20);

% Theoretical decay rate (from orbit-fit Delta a)
track_dadt = abs(track_delt_a)./Delt_t;
track_dadt_fill = filloutliers(track_dadt,'linear','movmedian',30);
track_dadt_fill = smoothdata(track_dadt_fill,'rlowess',30);

y_j = elem_dadt_fill;
y_star = track_dadt_fill;

end

function g = G_MLAT_sigmoid(mlat)
% MLAT activation function (sigmoid in |MLAT|)
    k = 1; lambda = 60;
    g = 1 ./ (1 + exp(-k .* (abs(mlat) - lambda)));
end

function h = H_mlt_gaussian(mlt)
% MLT activation function (wrapped Gaussian in MLT)
    mu = 22;    % Peak position (mean)
    sigma = 2;  % Width (standard deviation)
    dist = abs(mlt - mu);
    dist = min(dist, 24 - dist); % Key step: take shortest periodic distance
    % Apply Gaussian function
    h = exp(-dist.^2 / (2 * sigma^2));
end

function a = Amplitude_F107NDst_LowSolar(F107,ni,c)
% Amplitude modulation function in low solar activity
A1  = c(1);
k_F = c(2);
F0  = c(3);
p_n = c(4);
if A1 == 0
    A0 = 0;
else
    A0 = 15.1387;
end
x_F107 = F107;
a = A0 + ( A1 ./ ((ni).^p_n) ).*1./(1+exp(k_F.*(x_F107-F0)));
end

function DVeq = Delta_Veq_calc(MLAT,MLT,F107,ni,c)
% Compute extreme charging voltage increment from parameters
H_MLT = H_mlt_gaussian(MLT);
G_MLAT = G_MLAT_sigmoid(MLAT);
A_SunMag = Amplitude_F107NDst_LowSolar(F107,ni,c);
DVeq = A_SunMag.*G_MLAT.*H_MLT;
end

function Veq = Veq_calc(Te,Ti,MLAT,MLT,F107,ni,c)
% Compute equilibrium potential from parameters and solar/geomagnetic indices
k  = 1.380649e-23;  % Boltzmann constant (J/K) 
e  = 1.6021766e-19; % Elementary charge (C)
k_eV  = k/e;        % (eV/K)
mi = 15.999;
me = 5.48579909e-4; % Electron mass (u)

% Equilibrium potential at middle/low latitude
Veq0 = -k_eV .* Te ./ 2 .* log(mi .* Te ./ me ./ Ti);

% Increment of extreme charging equilibrium potential at polar region
DVeq = Delta_Veq_calc(MLAT,MLT,F107,ni,c);
Veq = abs(Veq0) + abs(DVeq);
end

function fion = IonAcc_calc(Te,Ne,Ni,sma,Veq)
% Input : Te [K], Ne [cm^-3], Ni [cm^-3], sma [km], Veq [V]
% Output: fion - quantity used in the integration to compute y_j

e  = 1.6021766e-19; % Elementary charge (C)
mi = 15.999;
u2kg  = @(m) 1.66053906660e-27 .* m;  % Convert atomic mass unit to kg
mi_kg = u2kg(mi);                     % Ion mass (kg)

% Satellite parameters
mu = 398601;        % (km3/s2)  
As_all = 1.4*1.2*6+2*1.4*2*2;
p0 = 1e-5;
ms = 640;
vs  = sqrt(mu ./ sma) * 1e3;  % Orbital speed (m/s)

lambda_D = 6.9 .* sqrt(Te ./ Ne);  % Debye length (cm)
C_Di = 4 .* (e .* Veq ./ (mi_kg .* vs.^2)).^2 .* log(lambda_D ./ p0);

% Quantity used for target function: k * N * Veq^2
fion = C_Di .* 0.5 .* As_all/ms .* mi_kg .* (Ni .* 1e6) .* vs.^2;

end

