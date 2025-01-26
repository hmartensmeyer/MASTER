% Import functions
addpath '../utils/'

% import dataset
load '..'\data\ceiling_vid_mobile.mat
data = grayscaleVideo_short;

%% Wavelet analysis

[original_flow, wavelet_coefficients_full, filtered_all_structures, filtered_dimples] = wavelet_func(data, 100, 15, 7, 0.2, 0.85, 0.6, 0.6);

%%

for t = 1:20
    imagesc(filtered_all_structures(:, :, t));
    colormap 'gray';
    title('Timestep: ', t)
    pause(0.1);
end

%% track dimples

[centroid_positions, structure_labels, structure_lifetimes] = dimpletracker(filtered_dimples, 15);

%% Visualization of tracking over time
figure;
for t = 1:100
    % Display the filtered structure for the current timestep
    subplot(1, 1, 1);
    imagesc(filtered_dimples(:, :, t));
    colormap('gray');
    hold on;

    % Plot centroids with unique colors for tracked structures
    if ~isempty(centroid_positions{t})
        for i = 1:size(centroid_positions{t}, 1)
            label = structure_labels{t}(i);
            color = lines(max(cellfun(@max, structure_labels)));  % Generate unique colors
            scatter(centroid_positions{t}(i, 1), centroid_positions{t}(i, 2), ...
                20, color(label, :), 'filled');  % Smaller markers
        end
    end

    title(sprintf('Timestep: %d', t));
    xlabel('X Coordinate');
    ylabel('Y Coordinate');
    hold off;
    pause(0.1);  % Pause for visualization
end

%% Display the list of structures and their lifetimes
disp('Structure Lifetimes:');
for i = 1:num_structures
    fprintf('Structure %d: %d timesteps\n', i, structure_lifetimes(i));
end
