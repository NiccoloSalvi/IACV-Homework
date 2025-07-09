% Image Analysis and Computer Vision - Homework A.Y. 2024/25 - G1
% Author: Salvi Niccolò (10773726)

% read the provided image
original_image = imread('images\look-outCat.jpg');

% display the original image
figure(1);
imshow(original_image);
hold on;
title('Original Image');

num_l_lines = 3; % l_1, l_2, l_3 lines
l_lines = nan(num_l_lines, 3);
for count = 1:num_l_lines
    figure(1);
    title(['Draw l_', num2str(count), ' line']);

    % user draws a line on the image
    seg = drawline('Color', 'r');
    l = segToLine(seg.Position);

    % store the line coefficients
    l_lines(count, :) = l';
end

van_pts = [];
for i = 1:num_l_lines-1
    for j = i+1:num_l_lines
        % disp("Computing vanishing point between l_" + i + " and l_" + j); % debug
        vp = cross(l_lines(i, :)', l_lines(j, :)');
        vp = vp ./ vp(3);
        van_pts = [van_pts; vp(1), vp(2)]; %#ok<*AGROW>
    end
end
van_point_l = [mean(van_pts, 1), 1];
disp("Vanishing Point l lines:");
disp(van_point_l);

num_m_lines = 6; % m_1, m_2, m_3, m4, m5, m6 lines
m_lines = nan(num_m_lines, 3);
for count = 1:num_m_lines
    figure(1);
    title(['Draw m_', num2str(count), ' line']);

    % user draws a line on the image
    seg = drawline('Color', 'b');
    m = segToLine(seg.Position);

    % store the line coefficients
    m_lines(count, :) = m';
end

van_pts = [];
for i = 1:num_m_lines-1
    for j = i+1:num_m_lines
        % disp("Computing vanishing point between m_" + i + " and m_" + j); % debug
        vp = cross(m_lines(i, :)', m_lines(j, :)');
        vp = vp ./ vp(3);
        van_pts = [van_pts; vp(1), vp(2)];
    end
end
van_point_m = [mean(van_pts, 1), 1];
disp("Vanishing Point m lines:");
disp(van_point_m);

% Compute the line at infinity
l_infty = cross(van_point_l, van_point_m);
l_infty = l_infty ./ l_infty(3);

disp("Compute line at infinity between l lines and m lines");
disp(l_infty);

title('Line at Infinity');
line([van_point_l(1), van_point_m(1)], [van_point_l(2), van_point_m(2)], 'Color', 'cyan', 'Linewidth', 2);

function lineH = segToLine(pts)
    a = [pts(1, :)'; 1];
    b = [pts(2, :)'; 1];
    lineH = cross(a,b);
    lineH = lineH./norm(lineH);
end