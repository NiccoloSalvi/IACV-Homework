% Image Analysis and Computer Vision - Homework A.Y. 2024/25 - F1
% Author: Salvi Niccolò (10773726)

% if inputImage is a filename, otherwise skip imread if you already have the image in a variable.
Ioriginal = imread('images\Look-outCat.jpg');

% Convert to grayscale if it’s a color image
if size(Ioriginal, 3) == 3
    Igray = rgb2gray(Ioriginal);
else
    Igray = Ioriginal;
end

% Define the region of interest (ROI)
roi = [280, 200, 1200, 590];  % [x, y, width, heightt

Igray = imcrop(Igray, roi);
Ioriginal = imcrop(Ioriginal, roi);

%% 2) Edge Detection
% Use the Canny method for robust edge detection
BW = edge(Igray, 'canny', [0.05, 0.2]);

% save BW
imwrite(BW, 'BW.jpg');

%% 3) Hough Transform
% The output H is the accumulator array, theta is the angle vector,
% and rho is the distance vector in Hough space
[H, theta, rho] = hough(BW, 'RhoResolution', 1, 'ThetaResolution', 1);

%% 4) Find Peaks in the Hough Transform
% houghpeaks locates the peaks (local maxima) in the accumulator array.
% Adjust 'NumPeaks' to capture as many lines as you need, 
% or use 'Threshold' to set a minimum accumulator value.
numPeaks = 20;  % Example: detect up to 25 peaks
P = houghpeaks(H, numPeaks, 'Threshold', ceil(0.3 * max(H(:))));

%% 5) Extract Lines
% 'houghlines' uses the accumulator peaks to generate line segments.
% 'FillGap' = max gap between two line segments to link them
% 'MinLength' = minimum line segment length
lines = houghlines(BW, theta, rho, P, 'FillGap', 200, 'MinLength', 300);

%% 6) Visualize Detected Lines
figure('Name','Detected Lines');

% Show the original image
imshow(Ioriginal); 
hold on;

% Optionally overlay the edge map if you want (comment out if not needed):
% imshow(BW);

title('Detected Lines Overlaid on Original Image');

% Loop through all detected lines and plot them
for k = 1:length(lines)
    xy = [lines(k).point1; lines(k).point2];
    
    % Plot the line segment
    plot(xy(:,1), xy(:,2), 'LineWidth', 2, 'Color', 'red');
    
    % Optionally mark endpoints for clarity
    plot(xy(1,1), xy(1,2), 'x', 'LineWidth', 2, 'Color', 'yellow');
    plot(xy(2,1), xy(2,2), 'x', 'LineWidth', 2, 'Color', 'green');
end

% Save the original image with plotted lines
saveas(gcf, 'images\F1_detected_lines.jpg');