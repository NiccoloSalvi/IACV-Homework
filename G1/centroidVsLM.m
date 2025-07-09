% Image Analysis and Computer Vision - Homework A.Y. 2024/25 - LM
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

van_point_l_lm = estimateVanishingPoint_LM(l_lines);
disp("Vanishing Point l lines (Levenberg-Marquardt):");
disp(van_point_l_lm);

% lines: Nx3 array of lines [a b c] in homogeneous form.
% find [vx, vy] that minimize sum of (a*vx + b*vy + c)^2.
function vx_vy = estimateVanishingPoint_LM(lines)
    vx_vy = [0; 0]; % initial guess for the vanishing point

    maxIterations = 100;
    lambda = 1e-3; % damping parameter
    tolerance = 1e-12;

    for iter = 1:maxIterations
        [residuals, J] = residualAndJacobian(vx_vy, lines);

        % compute the normal equations (J'J and J'r)
        A = J'*J;
        g = J'*residuals;

        % check for convergence
        if norm(g) < tolerance
            break;
        end

        % Levenberg-Marquardt step
        % solve (J'J + lambda * I) * delta = -J'r
        delta = -(A + lambda * eye(2)) \ g;

        % evaluate the new parameters
        new_vx_vy = vx_vy + delta;
        [new_residuals, ~] = residualAndJacobian(new_vx_vy, lines);

        % check if improvement
        if sum(new_residuals.^2) < sum(residuals.^2)
            % accept the step
            vx_vy = new_vx_vy;
            % decrease lambda to go more towards Gauss-Newton
            lambda = lambda / 10;
        else
            % reject the step and increase lambda
            lambda = lambda * 10;
        end

        % check for small step
        if norm(delta) < tolerance
            break;
        end
    end
end

function [r, J] = residualAndJacobian(vx_vy, lines)
    vx = vx_vy(1);
    vy = vx_vy(2);
    a = lines(:,1);
    b = lines(:,2);
    c = lines(:,3);

    % Residuals: r_i = a_i*vx + b_i*vy + c_i
    r = a*vx + b*vy + c;

    % Jacobian wrt [vx; vy] is: [a_i, b_i]
    J = [a, b];
end

van_point_l_lm = estimateVanishingPoint_LM(l_lines);

function lineH = segToLine(pts)
    a = [pts(1, :)'; 1];
    b = [pts(2, :)'; 1];
    lineH = cross(a,b);
    lineH = lineH./norm(lineH);
end