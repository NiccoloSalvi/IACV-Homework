% Image Analysis and Computer Vision - Homework A.Y. 2024/25 - G2
% Author: Salvi Niccolò (10773726)

% use the metric rectified of the original image (look-outCat.jpg) obtained from G2 (Hr.m)
image = im2double(imread('images\G2_horizontal_metric_rectified.jpg'));
imshow(image);
hold on;

% draw two lines on the metric rectified image, to compute the depth of the parallelepiped
title("Draw l1 in the metric rectified image");
seg_l1 = drawline('Color', 'r');
l1 = segToLine(seg_l1.Position);

title("Draw hypotenuse in the metric rectified image");
seg_m1 = drawline('Color', 'b');
m1 = segToLine(seg_m1.Position);

% compute m
cos_alpha = (l1' * [1 0 0; 0 1 0; 0 0 0] * m1) / (sqrt((l1' * [1 0 0; 0 1 0; 0 0 0] * l1) * (m1' * [1 0 0; 0 1 0; 0 0 0] * m1)));
sin_alpha = sqrt(1 - cos_alpha^2);

length_l = 1;
value_m = length_l * sin_alpha;
disp("Estimated m value: " + value_m);