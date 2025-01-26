function [centroid_positions, structure_labels, structure_lifetimes] = dimple_tracker(filtered_dimples, max_distance)
% DIMPLE_TRACKER - Track filtered regions using a line-based distance approach
%
%   Input:
%       filtered_dimples - 3D binary array (x, y, t) representing filtered regions across timesteps.
%       max_distance     - Maximum distance to associate centroids between frames.
%
%   Output:
%       centroid_positions - Cell array where each cell contains centroids of regions at each timestep.
%       structure_labels   - Cell array where each cell contains labels for regions at each timestep.
%       structure_lifetimes - Array where each element represents the lifetime of a structure.
%
%   Description:
%       This function tracks regions (e.g., dimples) across multiple timesteps by computing the
%       distance between region centroids in consecutive frames. Regions in the current frame
%       are associated with regions in the previous frame if the distance between their centroids
%       is below the specified threshold. It also calculates the lifetime of each structure.

% Number of timesteps
num_timesteps = size(filtered_dimples, 3);

% Preallocate outputs
centroid_positions = cell(num_timesteps, 1);  % Store centroids for each timestep
structure_labels = cell(num_timesteps, 1);   % Store region labels for tracking

% Loop through each timestep to extract region centroids and labels
for t = 1:num_timesteps
    disp(['Processing timestep: ', num2str(t)])
    
    % Get the binary mask for the current timestep
    binary_mask = filtered_dimples(:, :, t) > 0;

    % Label connected components in the binary mask
    connected_components = bwconncomp(binary_mask);

    % Measure region properties (centroids)
    region_props = regionprops(connected_components, 'Centroid');

    % Store centroids and labels for the current timestep
    if ~isempty(region_props)
        centroids = cat(1, region_props.Centroid); % Extract centroids
        centroid_positions{t} = centroids;

        % Assign labels to structures based on line-based distance to previous frame
        if t > 1 && ~isempty(centroid_positions{t-1})
            prev_centroids = centroid_positions{t-1};
            structure_labels{t} = zeros(size(centroids, 1), 1);

            % Check each centroid in the current frame
            for i = 1:size(centroids, 1)
                matched = false;
                for j = 1:size(prev_centroids, 1)
                    % Calculate distance between current and previous centroids
                    distance = norm(centroids(i, :) - prev_centroids(j, :));

                    % If distance is less than max_distance, match the structure
                    if distance < max_distance
                        structure_labels{t}(i) = structure_labels{t-1}(j);  % Same structure
                        matched = true;
                        break;
                    end
                end

                % Assign a new label if no match is found
                if ~matched
                    if t > 1
                        max_label = max(structure_labels{t-1}, [], 'omitnan');
                    else
                        max_label = 0;
                    end
                    structure_labels{t}(i) = max_label + 1; % New label
                end
            end
        else
            % Assign new labels for the first frame or if no previous centroids exist
            structure_labels{t} = (1:size(centroids, 1))';
        end
    else
        % No regions in the current timestep
        centroid_positions{t} = [];
        structure_labels{t} = [];
    end
end

%% Calculate lifetime of each structure
num_structures = max(cellfun(@max, structure_labels));  % Total number of structures
structure_lifetimes = zeros(num_structures, 1);  % Preallocate lifetimes array

% Loop through all timesteps to count occurrences of each structure
for t = 1:num_timesteps
    if ~isempty(structure_labels{t})
        unique_labels = unique(structure_labels{t});  % Get unique structure labels
        for label = unique_labels'
            if label > 0  % Ignore unmatched structures (label = 0)
                structure_lifetimes(label) = structure_lifetimes(label) + 1;
            end
        end
    end
end
end
