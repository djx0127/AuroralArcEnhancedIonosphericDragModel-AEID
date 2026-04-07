function [r_eci, v_eci, v_norm, tht] = element2eci(sma,ecc,inc,OMG,aop,M,r)
% Convert Keplerian orbital elements to ECI coordinates
% Input:
%   sma : semi-major axis        [km]
%   ecc : eccentricity           [-]
%   inc : inclination            [deg]
%   OMG : right ascension of ascending node (RAAN) [deg]
%   aop : argument of perigee    [deg]
%   M   : mean anomaly           [deg]
%   r   : orbital radius (optional) [km]
% Output:
%   r_eci : position vector in ECI frame [m]
%   v_eci : velocity vector in ECI frame [m/s]
%   v_norm: velocity magnitude           [m/s]
%   tht   : true anomaly                 [rad]

% Earth's gravitational parameter (km^3/s^2)
mu = 398600.4418;

% Convert angles from degrees to radians
inc = deg2rad(inc);
OMG = deg2rad(OMG);
aop = deg2rad(aop);
M   = deg2rad(M);

% 1. Solve Kepler's equation for eccentric anomaly E
E = solveKeplerEquation(M, ecc);

% 2. Compute true anomaly tht
tht = 2 * atan2(sqrt(1+ecc) * sin(E/2), sqrt(1-ecc) * cos(E/2));
% tht = 2 * atan(sqrt((1+ecc)/(1-ecc))*tan(E/2));

% If radius not provided, compute it from SMA and eccentric anomaly
if nargin < 7 || isempty(r)
    r = sma * (1 - ecc * cos(E)); % [km]
end

% 3. Position in the orbital plane coordinates
r_orb = [r ; 0; 0];

% Compute velocity components in the orbital plane
p = sma * (1 - ecc^2); % semi-latus rectum
h = sqrt(mu * p);      % specific angular momentum
vr = sqrt(mu/p)*ecc*sin(tht);       % radial component
vu = sqrt(mu/p)*(1+ecc*cos(tht));   % transverse component
v_norm = sqrt(vr^2+vu^2);
v_orb = [vr; vu; 0];

% 5. Rotation matrix (from orbital plane to ECI)
u = aop + tht;
Rz_OMG = [cos(OMG), -sin(OMG),  0;
          sin(OMG),  cos(OMG),  0;
                 0,         0,  1];
Rx_inc = [1, 0,          0;
          0, cos(inc), -sin(inc);
          0, sin(inc),  cos(inc)];
Rz_aop = [cos(u), -sin(u),  0;
          sin(u),  cos(u),  0;
               0,       0,  1];

R = Rz_OMG * Rx_inc * Rz_aop;

% 6. Transform to ECI frame
r_eci = R * r_orb * 1e3; % m
v_eci = R * v_orb * 1e3; % m/s
v_norm = v_norm * 1e3;   % m/s

end

function E = solveKeplerEquation(M, ecc)
% Solve Kepler's equation M = E - ecc*sin(E) using Newton-Raphson iteration
    max_iter = 50;
    tol = 1e-12;
    
    % Initial guess
    if M < pi
        E0 = M + ecc/2;
    else
        E0 = M - ecc/2;
    end
    
    for i = 1:max_iter
        f       = E0 - ecc*sin(E0) - M;
        f_prime = 1 - ecc*cos(E0);
        
        E1 = E0 - f/f_prime;
        
        if abs(E1 - E0) < tol
            E = E1;
            return;
        end
        
        E0 = E1;
    end
    
    % If not converged, use the last iteration value
    E = E1;
end