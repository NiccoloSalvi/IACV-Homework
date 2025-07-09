% Image Analysis and Computer Vision - Homework A.Y. 2024/25 - G5
% Author: Salvi Niccolò (10773726)

% use the metric rectified of the original image obtained from G2 (Hr.m)
metric_rectified_image = imread('images\G2_horizontal_metric_rectified.jpg');

% Read the curve S points from the file obtained in F3 (S.m)
try
    data = readmatrix('data\curve_S_points.csv');
catch
    error('Error reading curve_points.csv. Please check if the file exists and is in the correct format.');
end
xPoints = data(:, 1); % First column
yPoints = data(:, 2); % Second column

% known ground truth points in the original image (in pixels)
originalPoints = 1e03 * [
    0.4635, 0.4527;
    0.3835, 0.23075;
    1.3206, 0.3701;
    0.9432, 0.4868;
];

% known ground truth points in the rectified image (in pixels)
rectifiedPoints = 1e03 * [
    3.752, 9.4292;
    3.560, 8.338;
    6.522, 7.735;
    5.726, 9.048;
];

% compute the transformation between original and rectified points
tformRef = fitgeotrans(originalPoints, rectifiedPoints, 'projective');

% use this transformation to adjust rectified points
[xAdjusted, yAdjusted] = transformPointsForward(tformRef, xPoints, yPoints);

% Print a dozen random points from the rectified curve
numPoints = length(xAdjusted);
randomIndices = randperm(numPoints, 12);
randomPoints = [xAdjusted(randomIndices), yAdjusted(randomIndices)];

disp('Random points from the rectified curve:');
disp(randomPoints);

figure;
imshow(metric_rectified_image);
hold on;
plot(xAdjusted, yAdjusted, 'ro', 'MarkerSize', 2);

% save image to file
% saveas(gcf, 'images\G5_horizontal_metric_rectified_with_S.jpg');