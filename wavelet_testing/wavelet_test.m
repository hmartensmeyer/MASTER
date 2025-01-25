%% Kun development. Må gjøre om denne til en funksjon tenker jeg.

eta = load ('..\data\ceiling_vid_mobile.mat');
eta = eta.grayscaleVideo_short;
%%

t_index = 145;
snapshot = eta(:, :, t_index);

% Perform 2D continuous wavelet transform with the Mexican hat wavelet
scales = 1:10;  % Adjust scale range based on feature size
cwt_result = cwtft2(snapshot, 'Wavelet', 'mexh', 'Scales', scales);

% Extract wavelet coefficients at a specific scale
selected_scale = 7;  % Example scale index
wavelet_coefficients = cwt_result.cfs(:,:,selected_scale);

% Define the threshold
W_thr = 0.15;

% Create a mask for regions where W > W_thr
mask = wavelet_coefficients > W_thr;

% Apply the mask to the wavelet coefficients
filtered_coefficients = wavelet_coefficients .* mask;

% Label connected regions in the binary mask
connected_components = bwconncomp(mask);

% Measure properties of connected regions
region_props = regionprops(connected_components, 'Eccentricity','Area', 'Solidity');

% Create a new mask for regions with eccentricity < 0.85 or circularity >
% 0.85
eccentricity_threshold = 0.85;
solidity_threshold = 0.6;
eccentric_regions = ismember(labelmatrix(connected_components), ...
    find([region_props.Eccentricity] < eccentricity_threshold  & [region_props.Solidity] > solidity_threshold));

% Apply the new mask to the wavelet coefficients
filtered_by_eccentricity = wavelet_coefficients .* eccentric_regions;
%filtered_by_eccentricity = 1 - eccentric_regions; %binary version

% Plot original surface elevation
figure('Name', 'Wavelet Analysis', 'Position', [100, 100, 1000, 650], Colormap=gray);

% Original surface elevation
subplot(2, 2, 1);
imagesc(snapshot);
title(sprintf('Surface Elevation at t = %d', t_index));
xlabel('X Coordinate');
ylabel('Y Coordinate');
colorbar;

% Wavelet coefficients
subplot(2, 2, 2);
imagesc(wavelet_coefficients);
title(sprintf('Wavelet Coefficients (Scale %d)', scales(selected_scale)));
xlabel('X Coordinate');
ylabel('Y Coordinate');
colorbar;

% Visualization of filtered coefficients
subplot(2, 2, 3);
imagesc(filtered_coefficients);
title('Filtered Wavelet Coefficients (W < W_{thr})');
xlabel('X Coordinate');
ylabel('Y Coordinate');
colorbar;

% Visualization of the filtered coefficients
subplot(2, 2, 4);
imagesc(filtered_by_eccentricity);
title('Filtered Coefficients (Eccentricity < 0.85)');
xlabel('X Coordinate');
ylabel('Y Coordinate');
colorbar;

% Calculate coverage after W-thresholding
nonzero_pixels = nnz(filtered_coefficients);  % Number of nonzero pixels
total_pixels = numel(filtered_coefficients); % Total number of pixels
coverage = nonzero_pixels / total_pixels;    % Fraction of nonzero pixels

% Display coverage
fprintf('Coverage after W-thresholding: %.2f%%\n', coverage * 100);

% Display eccentricity values of filtered regions
fprintf('Eccentricity of filtered regions:\n');
for i = 1:length(region_props)
    fprintf('Region %d: Eccentricity = %.4f\n', i, region_props(i).Eccentricity);
end

% %% Try to find out which region is which
% 
% % Get the labeled matrix of connected components
% labeled_regions = labelmatrix(connected_components);
% 
% % Display region information along with eccentricity
% fprintf('Region Eccentricity Mapping:\n');
% for i = 1:length(region_props)
%     fprintf('Region %d: Eccentricity = %.4f, Area = %d, Solidity = %d\n', i, region_props(i).Eccentricity, region_props(i).Area, region_props(i).Solidity);
% end
% 
% % Example: visualize one specific region by its index
% region_index = 48; % Change this to the desired region index
% specific_region_mask = (labeled_regions == region_index);
% 
% % Plot the specific region
% figure('Name', sprintf('Region %d Visualization', region_index), Colormap=gray);
% imagesc(specific_region_mask);
% title(sprintf('Region %d (Eccentricity = %.4f)', region_index, region_props(region_index).Eccentricity));
% xlabel('X Coordinate');
% ylabel('Y Coordinate');
% colorbar;


%% Test for time series

% Parameters
timesteps = 1:1000;  % Define the range of timesteps (100 timesteps)
scales = 1:10;  % Adjust scale range based on feature size
selected_scale = 7;  % Scale index to use
W_thr = 0.2;  % Threshold for wavelet coefficients
eccentricity_threshold = 0.85;  % Threshold for eccentricity
circularity_threshold = 0.8;
solidity_threshold = 0.6;

% Preallocate array for filtered snapshots
[x_dim, y_dim] = size(eta(:, :, 1));  % Dimensions of each snapshot
all_structures = zeros(x_dim, y_dim, length(timesteps));
filtered_snapshots = zeros(x_dim, y_dim, length(timesteps));  % 3D array to store results

% Loop through each timestep
for t_index = 1:length(timesteps)
    disp(t_index)
    t = timesteps(t_index);
    snapshot = eta(:, :, t);

    % Perform 2D continuous wavelet transform with the Mexican hat wavelet
    cwt_result = cwtft2(snapshot, 'Wavelet', 'mexh', 'Scales', scales);

    % Extract wavelet coefficients at the selected scale
    wavelet_coefficients = cwt_result.cfs(:, :, selected_scale);

    % Create a mask for regions where W < W_thr
    mask = wavelet_coefficients > W_thr;

    % Apply the mask to the wavelet coefficients
    filtered_coefficients = wavelet_coefficients .* mask;

    % Label connected regions in the binary mask
    connected_components = bwconncomp(mask);

    % Measure properties of connected regions
    region_props = regionprops(connected_components, 'Eccentricity', 'Solidity');

    % Create a new mask for regions with eccentricity < threshold
    eccentric_regions = ismember(labelmatrix(connected_components), ...
        find([region_props.Eccentricity] < eccentricity_threshold  & [region_props.Solidity] > solidity_threshold));

    % Apply the mask to the wavelet coefficients
    filtered_by_eccentricity = wavelet_coefficients .* eccentric_regions;
    %filtered_by_eccentricity = 1 - eccentric_regions; %binary version

    % Save the filtered snapshot
    all_structures(:, :, t_index) = filtered_coefficients;
    filtered_snapshots(:, :, t_index) = filtered_by_eccentricity;
end

% Save the filtered snapshots to a MAT file
%save('filtered_snapshots.mat', 'filtered_snapshots', '-v7.3');

%% Display the vortices over time

for t = 1:1000
    imagesc(filtered_snapshots(:, :, t));
    colormap 'gray';
    title('Timestep: ', t)
    pause(0.1);
end

%% Neste på programmet

% Forkaste alle strukturer som er mindre enn 5 tidssteg lange?
% Tracke virvlene, og kunne koordinatfeste dem?
% Hvordan kan jeg koble dette til hastighetsfeltet/virvlingsfeltet under

%% Track filtered regions using line-based distance approach
% Parameters for tracking
num_timesteps = size(filtered_snapshots, 3);
centroid_positions = cell(num_timesteps, 1);  % Store centroids for each timestep
structure_labels = cell(num_timesteps, 1);   % Store region labels for tracking
max_distance = 15;  % Maximum distance to associate centroids between frames

% Loop through each timestep to extract region centroids and labels
for t = 1:num_timesteps
    disp(t)
    % Get the binary mask for the current timestep
    binary_mask = filtered_snapshots(:, :, t) > 0;

    % Label connected components in the binary mask
    connected_components = bwconncomp(binary_mask);

    % Measure region properties (centroids)
    region_props = regionprops(connected_components, 'Centroid');

    % Store centroids and labels for the current timestep
    if ~isempty(region_props)
        centroids = cat(1, region_props.Centroid);
        centroid_positions{t} = centroids;

        % Assign labels to structures based on line-based distance to previous frame
        if t > 1 && ~isempty(centroid_positions{t-1})
            prev_centroids = centroid_positions{t-1};
            structure_labels{t} = zeros(size(centroids, 1), 1);

            % Check each centroid in the current frame
            for i = 1:size(centroids, 1)
                matched = false;
                for j = 1:size(prev_centroids, 1)
                    % Calculate line length between current and previous centroids
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
                    structure_labels{t}(i) = max(structure_labels{t-1}, [], 'omitnan') + 1;
                end
            end
        else
            structure_labels{t} = (1:size(centroids, 1))';  % Assign new labels for first frame
        end
    else
        centroid_positions{t} = [];
        structure_labels{t} = [];
    end
end

% %% Visualization of tracking over time
% figure;
% for t = 1:num_timesteps
%     % Display the filtered structure for the current timestep
%     subplot(1, 1, 1);
%     imagesc(filtered_snapshots(:, :, t));
%     colormap('gray');
%     hold on;
% 
%     % Plot centroids with unique colors for tracked structures
%     if ~isempty(centroid_positions{t})
%         for i = 1:size(centroid_positions{t}, 1)
%             label = structure_labels{t}(i);
%             color = lines(max(cellfun(@max, structure_labels)));  % Generate unique colors
%             scatter(centroid_positions{t}(i, 1), centroid_positions{t}(i, 2), ...
%                 20, color(label, :), 'filled');  % Smaller markers
%         end
%     end
% 
%     title(sprintf('Timestep: %d', t));
%     xlabel('X Coordinate');
%     ylabel('Y Coordinate');
%     hold off;
%     pause(0.41);  % Pause for visualization
% end

%% Visualize structures as small dots at each timestep
% Initialize figure
figure('Name', 'Structure Tracking with Dots', 'Position', [100, 100, 800, 600]);
hold on;
axis equal;
xlim([1, size(filtered_snapshots, 2)]);
ylim([1, size(filtered_snapshots, 1)]);
title('Tracked Structures Over Time');
xlabel('X Coordinate');
ylabel('Y Coordinate');
colormap('gray');

% Assign unique colors for structures
num_structures = max(cellfun(@max, structure_labels));
colors = lines(num_structures);

% Loop through timesteps to plot dots for structures
for t = 1:num_timesteps
    disp(['Processing timestep: ', num2str(t)]);  % Debug output

    % Get centroids and labels for the current frame
    if ~isempty(centroid_positions{t}) && ~isempty(structure_labels{t})
        for i = 1:size(centroid_positions{t}, 1)
            label = structure_labels{t}(i);
            if label > 0
                % Plot a small dot at the centroid
                scatter(centroid_positions{t}(i, 1), centroid_positions{t}(i, 2), ...
                    10, colors(label, :), 'filled');  % Small dots
            end
        end
    end
    pause(0.02);  % Pause for visualization
end

hold off;

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

%% Display the list of structures and their lifetimes
disp('Structure Lifetimes:');
for i = 1:num_structures
    fprintf('Structure %d: %d timesteps\n', i, structure_lifetimes(i));
end


%% Visualize a single structure and its active timesteps
structure_to_show = 536;  % Specify the structure label to visualize

% Initialize figure for visualization
figure('Name', sprintf('Structure %d Visualization', structure_to_show), 'Position', [100, 100, 800, 600]);

% Track the timesteps where the structure is present
active_timesteps = [];

for t = 1:num_timesteps
    if ~isempty(centroid_positions{t}) && ~isempty(structure_labels{t})
        % Check if the structure exists in this frame
        if ismember(structure_to_show, structure_labels{t})
            % Add to active timesteps
            active_timesteps = [active_timesteps, t];

            % Find the connected component for this structure
            binary_mask = filtered_snapshots(:, :, t) > 0;
            connected_components = bwconncomp(binary_mask);
            structure_idx = find(structure_labels{t} == structure_to_show, 1);

            if ~isempty(structure_idx)
                % Create a mask for the structure
                structure_mask = ismember(labelmatrix(connected_components), structure_idx);

                % Visualize the structure
                imagesc(structure_mask);
                colormap('gray');
                title(sprintf('Structure %d at Timestep %d', structure_to_show, t));
                xlabel('X Coordinate');
                ylabel('Y Coordinate');
                colorbar;
                pause(0.5);  % Pause for visualization
            end
        end
    end
end

% Display the active timesteps for the structure
disp(['Structure ', num2str(structure_to_show), ' is active during timesteps:']);
disp(active_timesteps);


%% Filter structures based on lifetime and visualize the full timeseries

lifetime_threshold = 6;  % Minimum lifetime for structures to be displayed

% Identify structures that meet the lifetime criterion
valid_structures = find(structure_lifetimes >= lifetime_threshold);

% Initialize figure for visualization
figure('Name', 'Filtered Structures Over Time', 'Position', [100, 100, 800, 600]);
hold on;
axis equal;
xlim([1, size(filtered_snapshots, 2)]);
ylim([1, size(filtered_snapshots, 1)]);
title('Filtered Structures Over Time');
xlabel('X Coordinate');
ylabel('Y Coordinate');
colormap('gray');

% Assign unique colors for valid structures
colors = lines(length(valid_structures));

% Loop through timesteps to plot valid structures
for t = 1:num_timesteps
    % Display the snapshot as background
    imagesc(filtered_snapshots(:, :, t));
    hold on;

    % Get centroids and labels for the current frame
    if ~isempty(centroid_positions{t}) && ~isempty(structure_labels{t})
        for i = 1:size(centroid_positions{t}, 1)
            label = structure_labels{t}(i);
            if ismember(label, valid_structures)
                % Plot a small dot for the valid structure
                structure_idx = find(valid_structures == label, 1);
                scatter(centroid_positions{t}(i, 1), centroid_positions{t}(i, 2), ...
                    20, colors(structure_idx, :), 'filled');
            end
        end
    end

    pause(0.1);  % Pause for visualization
    hold off;
end
