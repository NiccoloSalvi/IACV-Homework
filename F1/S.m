% Image Analysis and Computer Vision - Homework A.Y. 2024/25 - F3
% Author: Salvi Niccolò (10773726)

I = imread('images\look-outCat.jpg');

if size(I,3) == 3
    Igray = rgb2gray(I);
else
    Igray = I;
end

% Preprocess: Denoise and enhance contrast
IgraySmooth = imgaussfilt(Igray, 1); % Gaussian smoothing
Ienhanced = imadjust(IgraySmooth);   % Contrast enhancement

% Detect edges using Canny
BW = edge(Ienhanced, 'canny', [0.20, 0.40], 0.9); % Thresholds and sigma

% Postprocess: Clean and close gaps
BWclean = bwareaopen(BW, 50); % Remove small noisy edges
BWclosed = imclose(BWclean, strel('disk', 1)); % Close gaps in edges

%% 3) Extract (x, y) Edge Coordinates
% Find all nonzero pixels in BW
[yCoords, xCoords] = find(BWclosed);

% define a ROI for the S curve
x1 = 830;
x2 = 955;
y1 = 510;
y2 = 545;
mask = (xCoords >= x1) & (xCoords <= x2) & (yCoords >= y1) & (yCoords <= y2);

xCoords = xCoords(mask);
yCoords = yCoords(mask);

% Convert to double for computations
xCoords = double(xCoords);
yCoords = double(yCoords);

if numel(xCoords) < 10
    warning('Not enough edge points detected. Aborting.');
    return;
end

%% 4) RANSAC Fitting of Ellipse
% Main parameters
maxIter = 100;   % number of RANSAC iterations
distThresh = 10;    % inlier distance threshold (pixels, approx)

[bestParams, bestInliers] = ransacEllipse(xCoords, yCoords, maxIter, distThresh);

if isempty(bestParams)
    warning('No valid ellipse found via RANSAC.');
    return;
end

% If we want to refine by using all inliers in a final least-squares fit:
xIn = xCoords(bestInliers);
yIn = yCoords(bestInliers);
finalParams = fitEllipseDirect(xIn, yIn);  % might refine from bestParams

%% 5) Visualization
figure; imshow(I); hold on;
plot(xCoords, yCoords, 'b', 'MarkerSize', 1);   % raw edge points

% find the mean of x
xMean = mean(xCoords);

% Find the point with the maximum y value where xCoords > xMean
validPoints = xCoords > xMean;
if any(validPoints)
    [maxY, idx] = max(yCoords(validPoints));
    maxX = xCoords(validPoints);
    maxX = maxX(idx);
    lowRight = [maxX, maxY];
    % plot(maxX, maxY, 'go', 'MarkerSize', 10, 'LineWidth', 2);
end

validPoints = xCoords < xMean;
if any(validPoints)
    [maxY, idx] = max(yCoords(validPoints));
    maxX = xCoords(validPoints);
    maxX = maxX(idx);
    lowLeft = [maxX, maxY];
    % plot(maxX, maxY, 'go', 'MarkerSize', 10, 'LineWidth', 2);
end

% Plot the ellipse from finalParams:
plotEllipse(finalParams, 'r-', 2, lowLeft, lowRight);
title('RANSAC Ellipse Fit for Shelf Edge (S)');

%% 6) Extract Points on the Fitted Ellipse
% Extract fitted ellipse parameters
a = finalParams(1);  % Semi-major axis
b = finalParams(2);  % Semi-minor axis
cx = finalParams(3); % Center x-coordinate
cy = finalParams(4); % Center y-coordinate
theta = finalParams(5); % Rotation angle

% Choose 12 evenly spaced points on the curve
tValues = linspace(0, -pi, 100); % Adjust range if needed for partial curve

% Compute X-Y coordinates using parametric ellipse equations
xPoints = cx + a * cos(tValues) * cos(theta) - b * sin(tValues) * sin(theta);
yPoints = cy + a * cos(tValues) * sin(theta) + b * sin(tValues) * cos(theta);

% Display the results
% disp('X-Y Coordinates of Points on the Curve S:');
% for i = 1:numel(xPoints)
    % fprintf('Point %d: X = %.2f, Y = %.2f\n', i, xPoints(i), yPoints(i));
% end

% save points to file
points = [xPoints', yPoints'];
writematrix(points, 'data\curve_S_points.csv');

% plot(xPoints, yPoints, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5); % Overlay points
% title('Selected Points on Curve S');

% save iamge to file
% saveas(gcf, 'images\F3_detected_S.jpg');

% use spline to fit the curve from the xPoints and yPoints
x_s = linspace(min(xPoints), max(xPoints), 1000);
y_s = spline(xPoints, yPoints, x_s);
% plot the fitted curve
% plot(x_s, y_s, 'b', 'LineWidth', 2);

%==================================================================
function [bestParams, bestInliers] = ransacEllipse(xAll, yAll, maxIter, distThresh)
% RANSACELLIPSE  Robustly fits an ellipse using RANSAC.
%
%   Inputs:
%       xAll, yAll  - arrays of edge coordinates
%       maxIter     - number of RANSAC iterations
%       distThresh  - inlier distance threshold
%
%   Outputs:
%       bestParams  - ellipse parameters [a b cx cy theta] or any scheme
%       bestInliers - logical index of inliers in xAll,yAll

    N = numel(xAll);
    bestNumInliers = 0;
    bestParams = [];
    bestInliers = [];

    for i = 1:maxIter

        % Randomly sample a minimal set of points. 
        % For a unique ellipse, we ideally need 5 distinct points (non-collinear).
        subsetIdx = randperm(N, 5);  
        xSubset = xAll(subsetIdx);
        ySubset = yAll(subsetIdx);

        % Fit ellipse to these 5 points
        candParams = fitEllipseDirect(xSubset, ySubset);
        if isempty(candParams)
            continue;  % might have been a degenerate fit
        end

        % Compute distances of all points to the candidate ellipse
        distances = ellipseDistance(candParams, xAll, yAll);

        % Determine inliers
        inlierMask = (distances < distThresh);
        numInliers = sum(inlierMask);

        % Check if this is the best so far
        if numInliers > bestNumInliers
            bestNumInliers = numInliers;
            bestParams = candParams;
            bestInliers = inlierMask;
        end
    end

    % If bestNumInliers is too small, discard
    if bestNumInliers < 10
        bestParams = [];
        bestInliers = [];
    end
end

%==================================================================
function params = fitEllipseDirect(x, y)
% FITELLIPSEDIRECT  Fits an ellipse (in algebraic form) to point data
% using the "direct least squares" method (Fitzgibbon et al.).
%
%   Returns ellipse parameters in a convenient geometric form
%   [a b cx cy theta], where:
%       a = major axis radius
%       b = minor axis radius
%       (cx, cy) = center
%       theta = rotation angle of major axis
%
% If the data is degenerate or the fit fails, returns [].

    if numel(x) < 5
        params = [];
        return;
    end

    % Build the design matrix D for conic fitting: x^2, xy, y^2, x, y, 1
    D = [x.^2, x.*y, y.^2, x, y, ones(size(x))];
    
    % Solve the normal system D' * D * p = 0
    % subject to constraint 4ac - b^2 > 0 (for ellipse).
    % We'll do a simplified unconstrained approach here for demonstration.
    [~, ~, V] = svd(D, 0);
    p = V(:,end);  % last column => smallest singular value

    % p = [A B C D E F]' in "Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0"
    A = p(1); 
    B = p(2); 
    C = p(3); 
    D_ = p(4);
    E = p(5);
    F_ = p(6);

    % Check the discriminant for ellipse vs hyperbola vs parabola:
    % For ellipse, we need B^2 - 4AC < 0
    if (B^2 - 4*A*C) >= 0
        params = [];
        return;
    end

    % Convert conic parameters to geometric parameters.
    % This part is somewhat involved. A reference formula is:
    %   https://mathworld.wolfram.com/Ellipse.html (see "General Conic Section")
    %
    % Center (cx, cy)
    denom = (4*A*C - B^2);
    cx = (B*E - 2*C*D_) / denom;
    cy = (B*D_ - 2*A*E) / denom;

    % Translate to center for the axes computation
    % We'll shift coordinates so that (cx, cy) = (0,0).
    % Then the conic becomes A'x^2 + B'xy + C'y^2 = 1, after dividing by some factor.
    % There's a standard procedure for this. 
    % For brevity, let's do a simpler approach:
    % Evaluate translational invariants first:
    % (Detailed derivation omitted here—this is standard ellipse-fitting geometry.)
    %
    % The orientation angle:
    theta = 0.5 * atan2(B, (A - C));

    % Compute axis lengths by eigen decomposition or by direct formula:
    % We can rotate the conic to principal axes.
    % The length of semi-major axis (a) and semi-minor axis (b) can be found from:
    %   2 * sqrt(...) from certain expressions in A, B, C, D, E, F.
    %
    % Below is a condensed formula. For full derivation, see references:
    cosT = cos(theta);
    sinT = sin(theta);

    % Terms for rotation
    A_rot = A*cosT^2 + B*cosT*sinT + C*sinT^2;
    C_rot = A*sinT^2 - B*cosT*sinT + C*cosT^2;

    % We should compute the scale factor to normalize "=1".
    % at center => plugging (cx, cy) in:
    val = A*(cx^2) + B*(cx*cy) + C*(cy^2) + D_*cx + E*cy + F_;
    if val == 0
        params = [];
        return;
    end
    % Now a^2 = - (1/val) for the larger axis direction, etc.

    % The radius along principal axes
    % a^2 = -1/(2 * (some expression)) etc. For simplicity we do:
    %   x'^2/A_rot + y'^2/C_rot = 1 if A_rot, C_rot > 0
    scale = -1 / val;
    A_rot_scaled = A_rot * scale;
    C_rot_scaled = C_rot * scale;

    if (A_rot_scaled <= 0) || (C_rot_scaled <= 0)
        % Not an ellipse or something degenerate
        params = [];
        return;
    end

    a = sqrt(1 / A_rot_scaled);
    b = sqrt(1 / C_rot_scaled);

    % Make sure a >= b  (swap if needed)
    if b > a
        tmp = a; a = b; b = tmp;
        % Flip angle by 90 degrees
        theta = theta + pi/2;
    end

    % Normalize angle to range [0, 2*pi)
    theta = mod(theta, pi);

    % Return final geometric parameters
    params = [a, b, cx, cy, theta];
end

%==================================================================
function d = ellipseDistance(params, x, y)
% ELLIPSEDISTANCE  Compute approximate distance of points to an ellipse.
%
%   params = [a, b, cx, cy, theta] (major/minor radius, center, orientation)
%   x, y   = arrays of point coords
%
% The "distance" here can be an algebraic or approximate geometric distance.  
% For robust RANSAC gating, a simple approximate measure is enough.

    if isempty(params)
        d = inf(size(x));
        return;
    end

    a     = params(1);
    b     = params(2);
    cx    = params(3);
    cy    = params(4);
    theta = params(5);

    % Translate points to ellipse center
    xt = x - cx;
    yt = y - cy;

    % Rotate by -theta to align with ellipse principal axes
    cosT = cos(-theta);
    sinT = sin(-theta);
    xprime = xt*cosT - yt*sinT;
    yprime = xt*sinT + yt*cosT;

    % The ellipse in principal axes is x'^2/a^2 + y'^2/b^2 = 1.
    % A rough "distance" measure is the absolute difference from 1:
    % dist = | x'^2/a^2 + y'^2/b^2 - 1 | * a typical length scale
    eqVal = (xprime.^2)/(a^2) + (yprime.^2)/(b^2);
    d = abs(eqVal - 1) * min(a,b);  % scale so that "on ellipse" => 0

end

function plotEllipse(params, lineSpec, lineWidth, lowLeft, lowRight)
    if isempty(params)
        return;
    end

    % Extract ellipse parameters
    a = params(1);  % Semi-major axis
    b = params(2);  % Semi-minor axis
    cx = params(3); % Center x-coordinate
    cy = params(4); % Center y-coordinate
    theta = params(5); % Rotation angle

    % Parameterize the ellipse
    t = linspace(0, 2*pi, 500); % Full ellipse (adjust range below for partial)
    xUnit = a * cos(t); % X-coordinates before rotation
    yUnit = b * sin(t); % Y-coordinates before rotation

    % Rotate by theta and translate to center
    cosT = cos(theta);
    sinT = sin(theta);
    xRot = xUnit * cosT - yUnit * sinT + cx; % Rotated and translated x
    yRot = xUnit * sinT + yUnit * cosT + cy; % Rotated and translated y

    % Mask to select points within the desired range
    % Example: Include only points on the "lower half" of the ellipse
    validMask = (xRot >= lowLeft(1)-5) & (xRot <= lowRight(1)) & (yRot <= lowLeft(2));

    % Extract the valid points
    xPartial = xRot(validMask);
    yPartial = yRot(validMask);

    % Sort the points by x to avoid connecting backward
    [~, sortIdx] = sort(xPartial);
    xPartial = xPartial(sortIdx);
    yPartial = yPartial(sortIdx);

    % Plot the partial ellipse as an open figure
    plot(xPartial, yPartial, lineSpec, 'LineWidth', lineWidth);
end