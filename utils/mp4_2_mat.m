% MATLAB Script to Convert an MP4 Video to Grayscale, Crop, and Save as .mat
% File Input: MP4 video file
% File Output: MAT file containing cropped grayscale frames

% Specify the input and output files
inputFile = '..\data\ceiling_vid_mobile.mp4'; % Replace with the actual file name
outputFile = '..\data\ceiling_vid_mobile.mat'; % Output .mat file name

%% Create a VideoReader object
videoReader = VideoReader(inputFile);

% Get the total number of frames in the video
numFrames = floor(videoReader.Duration * videoReader.FrameRate);

% Preallocate a cell array to store cropped grayscale frames
grayscaleFrames = cell(1, numFrames);

% Loop through each frame of the video
maxFrames = 1000; % Limit the number of frames processed
frameIndex = 1;

while hasFrame(videoReader) && frameIndex <= maxFrames
    % Read the next frame
    frame = readFrame(videoReader);
    
    % Convert the frame to grayscale
    grayFrame = rgb2gray(frame);
    
    % Convert the grayscale frame to double
    grayFrame = im2double(grayFrame); % Normalizes values to the range [0, 1]
    
    % Crop the top 5% of the frame
    frameHeight = size(grayFrame, 1);
    cropHeight = floor(0.05 * frameHeight); % Calculate the top 5% height
    grayFrameCropped = grayFrame(cropHeight+1:end, :); % Remove the top 5% rows
    
    % Store the cropped grayscale frame
    grayscaleFrames{frameIndex} = grayFrameCropped;
    frameIndex = frameIndex + 1;
    disp(frameIndex)
end

% Convert the cell array to a 3D matrix (cropped height x width x numFrames)
frameHeightCropped = size(grayscaleFrames{1}, 1);
frameWidth = size(grayscaleFrames{1}, 2);
grayscaleVideo = zeros(frameHeightCropped, frameWidth, maxFrames);

for i = 1:numFrames
    grayscaleVideo(:, :, i) = grayscaleFrames{i};
end

%% Optionally, select only the first 1000 frames
grayscaleVideo_short = grayscaleVideo(:, :, 1:min(1000, numFrames));

%% Save the cropped grayscale video to a .mat file
save(outputFile, 'grayscaleVideo_short', '-v7.3');

disp(['Cropped grayscale video saved to ', outputFile]);
