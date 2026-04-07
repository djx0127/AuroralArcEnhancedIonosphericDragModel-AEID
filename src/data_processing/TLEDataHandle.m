function [time,julian_time, sma, ecc, inc, raan, argp, M, r] = TLEDataHandle(filename, start_date, end_date)
%==========================================================================
% Description:
%   This function reads a Two-Line Element (TLE) file and extracts the
%   corresponding orbital elements within a specified time range.
%   It parses the TLE format, converts epochs to datetime, computes the
%   semi-major axis, and outputs the filtered orbital data.
%
% Inputs:
%   filename    - Path to the TLE file
%   start_date  - Start date (datetime object)
%   end_date    - End date (datetime object)
%
% Outputs:
%   time         - Time array (datetime)
%   julian_time  - Julian date
%   sma          - Semi-major axis [km]
%   ecc          - Eccentricity [-]
%   inc          - Inclination [deg]
%   raan         - Right Ascension of Ascending Node [deg]
%   argp         - Argument of Perigee [deg]
%   M            - Mean anomaly [deg]
%   r            - Orbital radius [km]
%==========================================================================

%% -------------------- File Reading --------------------
% Read TLE file content and split into lines
data = fileread(filename);
lines = strsplit(data, '\n');

%% -------------------- Constants --------------------
mu = 398600.4418; % Earth's gravitational parameter [km^3/s^2]

%% -------------------- Initialization --------------------
time = [];
sma  = [];
ecc  = [];
inc  = [];
raan = [];
argp = [];
M    = [];

%% -------------------- TLE Parsing Loop --------------------
% Process TLE in pairs (Line 1 and Line 2)
for i = 1:2:length(lines)
    
    if i+1 > length(lines)
        break;
    end
    
    line1 = lines{i};
    line2 = lines{i+1};
    
    % Validate TLE format identifiers
    if ~strncmp(line1, '1 ', 2) || ~strncmp(line2, '2 ', 2)
        continue; % Skip invalid entries
    end
    
    %% -------- Epoch Extraction --------
    if length(line1) >= 32
        epoch_str = line1(19:32);
        
        % Parse epoch format: YYDDD.DDDDDDDD
        year_str = epoch_str(1:2);
        day_str  = epoch_str(3:end);
        
        year = str2double(year_str);
        day_of_year = str2double(day_str);
        
        % Convert to full year
        if year < 57
            year = year + 2000;
        else
            year = year + 1900;
        end
        
        % Convert day-of-year to datetime
        base_date = datetime(year, 1, 0);
        time_val = base_date + days(day_of_year);
        
        %% -------- Orbital Elements Extraction --------
        if length(line2) >= 63
            
            % Inclination [deg]
            inc_str = line2(9:16);
            inclination = str2double(inc_str);
            
            % RAAN [deg]
            raan_str = line2(18:25);
            right_ascension = str2double(raan_str);
            
            % Eccentricity (decimal point assumed)
            ecc_str = line2(27:33);
            eccentricity = str2double(['0.' ecc_str]);
            
            % Argument of perigee [deg]
            argp_str = line2(35:42);
            argument_of_perigee = str2double(argp_str);
            
            % Mean anomaly [deg]
            M_str = line2(44:51);
            mean_anomaly = str2double(M_str);
            
            % Mean motion [rev/day]
            mean_motion_str = line2(53:63);
            mean_motion = str2double(mean_motion_str);
            
            %% -------- Semi-major Axis Computation --------
            % Convert mean motion to rad/s
            n_rad_per_sec = mean_motion * 2 * pi / 86400;
            
            % Compute semi-major axis using Kepler's third law
            a = (mu / (n_rad_per_sec^2))^(1/3);
            
            %% -------- Store Data --------
            time = [time; time_val];
            sma  = [sma; a];
            ecc  = [ecc; eccentricity];
            inc  = [inc; inclination];
            raan = [raan; right_ascension];
            argp = [argp; argument_of_perigee];
            M    = [M; mean_anomaly];
        end
    end
end

%% -------------------- Time Filtering --------------------
% Select data within the specified time window
time_mask = (time >= start_date) & (time <= end_date);

time = time(time_mask);
sma  = sma(time_mask);
ecc  = ecc(time_mask);
inc  = inc(time_mask);
raan = raan(time_mask);
argp = argp(time_mask);
M    = M(time_mask);

%% -------------------- Time Conversion --------------------
% Convert datetime to Julian date
julian_time = juliandate(time);

%% -------------------- Radius Computation --------------------
% Compute orbital radius (approximate from orbital elements)
r = sma .* (1 - ecc.^2) ./ (1 + ecc .* cosd(argp));

end