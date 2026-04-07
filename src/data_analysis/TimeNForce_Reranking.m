function [time_rank, ft_rank] = TimeNForce_Reranking(time_str, ft_str)
% -------------------------------------------------------------------------
% TimeNForce_Reranking
%
% Reorganize and merge time and force sequences from multiple fields,
% then sort them in ascending order of time.
%
% Inputs:
%   time_str - structure containing time series (Julian day)
%   ft_str   - structure containing corresponding force values
%
% Outputs:
%   time_rank - merged and time-sorted time sequence
%   ft_rank   - corresponding force sequence after sorting
% -------------------------------------------------------------------------

% Get all field names from input structures
time_fields = fieldnames(time_str);
ft_fields   = fieldnames(ft_str);

% Check consistency between structures
if length(time_fields) ~= length(ft_fields)
    error('Mismatch between number of fields in time and force structures');
end

% Initialize merged arrays
all_times  = [];
all_forces = [];
field_indices = []; % Record source field index for each data point

% Merge all time and force data
for i = 1:length(time_fields)
    
    % Extract current time and force sequences
    current_times  = time_str.(time_fields{i});
    current_forces = ft_str.(ft_fields{i});
    
    % Check length consistency
    if length(current_times) ~= length(current_forces)
        error(['Time and force length mismatch in field: ', time_fields{i}]);
    end
    
    % Append to merged arrays
    all_times  = [all_times; current_times];
    all_forces = [all_forces; current_forces];
    
    % Track origin of data points (optional, for debugging or analysis)
    field_indices = [field_indices, repmat(i, 1, length(current_times))];
end

% Sort by time
[sorted_times, sort_idx] = sort(all_times);
sorted_forces = all_forces(sort_idx);

% Output sorted sequences
time_rank = sorted_times;
ft_rank   = sorted_forces;

end