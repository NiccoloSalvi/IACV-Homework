% Image Analysis and Computer Vision - Homework A.Y. 2024/25 - Find h
% Author: Salvi Niccolò (10773726)

% image = im2double(imread('verticalReconstructionK.jpg')); % use the image reconstructed with K
image = im2double(imread('verticalReconstructionMetric.jpg')); % use the image reconstructed with metric stratified method
imshow(image);
hold on;

title("Draw l in the metric rectified image");
seg_l1 = drawline('Color', 'r');
l1 = segToLine(seg_l1.Position);

title("Draw hypotenuse in the metric rectified image");
seg_m1 = drawline('Color', 'b');
m1 = segToLine(seg_m1.Position);

% compute h
cos_alpha = (l1' * [1 0 0; 0 1 0; 0 0 0] * m1) / (sqrt((l1' * [1 0 0; 0 1 0; 0 0 0] * l1) * (m1' * [1 0 0; 0 1 0; 0 0 0] * m1)));
disp("Estimated cos_alpha value: " + cos_alpha);

sin_alpha = sqrt(1 - cos_alpha^2);
length_l = 1; % knwon value

value_h = length_l * sin_alpha;
disp("Estimated h value: " + value_h);

function lineH = segToLine(pts)
    a = [pts(1, :)'; 1];
    b = [pts(2, :)'; 1];
    lineH = cross(a,b);
    lineH = lineH./norm(lineH);
end