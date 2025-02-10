%% I. Read Video and Prepare Timestamps

data = load ('..\data\filtered_gray_5000t_indices.mat');

video = data.filteredFramesGray;
times = data.filteredTimeindeces;

% Short snippet to get data on correct form
[height, width] = size(video{1});

% Preallocate a 3D matrix for the frames
numFrames = length(video);
eta = zeros(height, width, numFrames, 'uint8'); % Use 'uint8' for grayscale images

% Populate the 3D matrix
for t = 1:numFrames
    eta(:, :, t) = video{t};
end

disp('Data read and converted to correct form.');

%% II. Wavelet Analysis and Filtering
% Parameters for wavelet filtering
scales = 1:15;
selected_scale = 15;
W_thr = 90;
eccentricity_threshold = 0.85;
solidity_threshold = 0.6;
[x_dim, y_dim, ~] = size(eta);
filtered_all_structures = zeros(x_dim, y_dim, numFrames);
filtered_dimples = zeros(x_dim, y_dim, numFrames);

%% III. Tracking Structures
% Initialize tracking structure (adding an "active" flag)
tracks = struct('id', {}, 'centroids', {}, 'frames', {}, 'active', {});
nextTrackId = 1;
baseYShift = 35;  % Base offset per time unit

for t_index = 1:200
    currentTime = times(t_index);
    disp(currentTime)
    snapshot = eta(:, :, t_index);
    
    % Wavelet transform and filtering (same as before)
    cwt_result = cwtft2(snapshot, 'Wavelet', 'mexh', 'Scales', scales);
    wavelet_coefficients = cwt_result.cfs(:, :, selected_scale);
    mask = wavelet_coefficients > W_thr;
    filtered_coefficients = wavelet_coefficients .* mask;
    connected_components = bwconncomp(mask);
    region_props = regionprops(connected_components, 'Eccentricity', 'Solidity', 'Centroid');
    validIdx = find([region_props.Eccentricity] < eccentricity_threshold & ...
                      [region_props.Solidity] > solidity_threshold);
    eccentric_regions = ismember(labelmatrix(connected_components), validIdx);
    filtered_by_eccentricity = wavelet_coefficients .* eccentric_regions;
    
    filtered_all_structures(:, :, t_index) = filtered_coefficients;
    filtered_dimples(:, :, t_index) = filtered_by_eccentricity;
    
    % Extract centroids of valid regions
    if isempty(validIdx)
        centroids = [];
    else
        centroids = cat(1, region_props(validIdx).Centroid);  % Each row: [x y]
    end
    numDetections = size(centroids, 1);
    
    % IV. Identify/Match Structures (Tracking)
    numTracks = length(tracks);
    costMatrix = Inf(numDetections, numTracks);
    for i = 1:numDetections
        for j = 1:numTracks
            if ~tracks(j).active, continue; end  % Skip dead tracks
            lastTime = tracks(j).frames(end);
            dt = currentTime - lastTime;
            % Predicted position: same x, y shifted by dt*baseYShift
            predicted = [tracks(j).centroids(end, 1), tracks(j).centroids(end, 2) + dt * baseYShift];
            allowedDistance = dt * baseYShift;
            d = norm(centroids(i, :) - predicted);
            if d <= allowedDistance
                costMatrix(i, j) = d;
            end
        end
    end
    
    % Greedy assignment of detections to tracks
    detectionTrackIDs = zeros(numDetections, 1);
    assignments = [];
    if ~isempty(costMatrix)
        while true
            [minVal, idx] = min(costMatrix(:));
            if isinf(minVal), break; end
            [detIdx, trackIdx] = ind2sub(size(costMatrix), idx);
            assignments = [assignments; detIdx, trackIdx, minVal];  %#ok<AGROW>
            detectionTrackIDs(detIdx) = tracks(trackIdx).id;
            costMatrix(detIdx, :) = Inf;
            costMatrix(:, trackIdx) = Inf;
        end
    end
    
    % Update assigned tracks with new detections
    if ~isempty(assignments)
        for k = 1:size(assignments, 1)
            detIdx = assignments(k, 1);
            trackIdx = assignments(k, 2);
            tracks(trackIdx).centroids(end+1, :) = centroids(detIdx, :);
            tracks(trackIdx).frames(end+1) = currentTime;
        end
    end
    
    % Start new tracks for unassigned detections
    for i = 1:numDetections
        if detectionTrackIDs(i) == 0
            tracks(nextTrackId).id = nextTrackId;
            tracks(nextTrackId).centroids = centroids(i, :);
            tracks(nextTrackId).frames = currentTime;
            tracks(nextTrackId).active = true;
            detectionTrackIDs(i) = nextTrackId;
            nextTrackId = nextTrackId + 1;
        end
    end
    
    % Declare lost tracks as dead if not updated in current frame (and if they had multiple updates)
    for j = 1:length(tracks)
        if tracks(j).active && tracks(j).frames(end) < currentTime && numel(tracks(j).frames) >= 2
            tracks(j).active = false;
        end
    end
end

%% V. Return Tracks with Lifetime Above a Threshold
lifetimeThreshold = 15;  % Adjust as needed (in number of frames)
numTracks = length(tracks);
trackInfo = struct('id', {}, 'lifetime', {}, 'coordinates', {});
for i = 1:numTracks
    if ~isempty(tracks(i).frames)
        % Lifetime defined as number of frames the structure was tracked.
        lifetime = numel(tracks(i).frames);
        trackInfo(end+1) = struct('id', tracks(i).id, 'lifetime', lifetime, 'coordinates', {tracks(i).centroids});
    end
end

disp('Track Information:');
for i = 1:length(trackInfo)
    if trackInfo(i).lifetime > lifetimeThreshold
    fprintf('Track %d: Lifetime = %d frames\n', trackInfo(i).id, trackInfo(i).lifetime);
    disp(trackInfo(i).coordinates);
    end
end

