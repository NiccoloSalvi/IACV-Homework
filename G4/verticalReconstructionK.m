% Image Analysis and Computer Vision - Homework A.Y. 2024/25 - G4
% Author: Salvi Niccolò (10773726)

% load the provided image
image = im2double(imread('images\Look-outCat.jpg'));

% import values obtained in G3 (K.m)
K = [
    1.4896e03 0 0.0015e03;
    0 2.1845e03 -0.0002e03;
    0 0 0.0010e03
];

% import values obtained from G1 (line_at_infinity.m)
v_l = [3392.3; 666.4; 1]; % vanishing point of l lines
v_m = [581.4; 883.8; 1]; % vanishing point of m lines
v_h = [659.8; -1409.8; 1]; % vanishing point of h lines

inf_ml = cross(v_l, v_m);
inf_ml = inf_ml ./ inf_ml(3);

inf_lh = cross(v_l, v_h);
inf_lh = inf_lh ./ inf_lh(3); % use it

inf_mh = cross(v_m, v_h);
inf_mh = inf_mh ./ inf_mh(3);

omega = inv(K * K');

% Image of circular points
syms 'x'; syms 'y';
% define the equations
el = [omega(1,1), omega(1,2) * 2, omega(2,2), omega(1,3) * 2, omega(2,3) * 2, omega(3,3)]; % conic coefficients
eleq = el(1) .* x.^2 + el(2) .* x .* y + el(3) .* y.^2 + el(4) .* x + el(5) .* y + el(6);
lieq = inf_lh(1) .* x + inf_lh(2) .* y;
system = [eleq == 0, lieq == -1];

% find the two solutions, images of the circular points
S = solve(system, [x, y]);
s1 = double([S.x(1); S.y(1); 1]);
s2 = double([S.x(2); S.y(2); 1]);

Ccirc = s1 * s2' + s2 * s1';
Ccirc = Ccirc ./ norm(Ccirc);

[U, D, V] = svd(Ccirc);
D(3,3) = 1;
A = sqrt(D) \ V';

figure;
tformVert = projective2d(A');
Ivert = imwarp(image, tformVert);
imshow(Ivert);

% save image to file
% imwrite(Ivert, 'images\G4_vertical_reconstruction_K.jpg');