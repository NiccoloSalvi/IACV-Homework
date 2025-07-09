% Image Analysis and Computer Vision - Homework A.Y. 2024/25 - F2
% Author: Salvi Niccolò (10773726)

image = imread('images\look-outCat.jpg');

figure;
imshow(image);
hold on;

% instead of using getpts, we can use the following points
x = [417.9118; 437.6296; 456.3276; 585.1732; 721.4981; 724.8977]; 
y = [566.9941; 581.6125; 526.1986; 510.9003; 543.1967; 592.1513];

% conic parameters
A = [x.^2 x.*y y.^2 x y ones(size(x))];
[~, ~, V] = svd(A);
cc = V(:, end);

% conic matrix
[a, b, c, d, e, f] = deal(cc(1), cc(2), cc(3), cc(4), cc(5), cc(6));
conic = [a b/2 d/2; b/2 c e/2; d/2 e/2 f];
conic = conic ./ conic(3, 3);

% conic points
[height, width, ~] = size(image);
[X, Y] = meshgrid(1:width, 1:height);
Z = zeros(height, width);

for i = 1:height
    for j = 1:width
        Z(i,j) = [j i 1] * conic * [j i 1]';
    end
end

scatter(x, y, 50, 'filled', 'yellow');
contour(X, Y, Z, [0 0], 'r', 'LineWidth', 2);

% save a, b, c, d, e, f in a text file
fileID = fopen('data\conic_C_parameters.txt', 'w');
fprintf(fileID, 'a: %.15f\nb: %.15f\nc: %.15f\nd: %.15f\ne: %.15f\nf: %.15f\n', a, b, c, d, e, f);
fclose(fileID);

% save image with detected conic
% saveas(gcf, 'images\F2_detected_conic.jpg');