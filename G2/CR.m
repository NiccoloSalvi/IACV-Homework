% Image Analysis and Computer Vision - Homework A.Y. 2024/25 - CROSS RATIO
% Author: Salvi Niccolò (10773726)

% use the metric rectified of the original image (look-outCat.jpg) obtained from G2 (Hr.m)
img = imread("images\look-outCat.jpg");

% visualize image
figure;
imshow(img);
hold on;
title("Metric Rectified Image");

[num_points, pts] = ginput(4);

A = pts(1,1);
B = pts(2,1);
C = pts(3,1);
D = pts(4,1);

% compute cross ratio(A, B, C, D)
cross_ratio_metric = ((C - A) / (C - B)) / ((D - A) / (D - B));
disp("CR: " + cross_ratio_metric);

close all;