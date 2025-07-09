% Image Analysis and Computer Vision - Homework A.Y. 2024/25 - G3
% Author: Salvi Niccolò (10773726)

% load rectification mapping obtained from G2 (Hr.m)
H_met = [
    1.4243 -6.9489e-1 0;
    -6.9489e-1 2.1380 0;
    -1.0000e-4 -1.1000e-3 1
];

% vanishing point h lines, computed in the same way as done in G1 (line_at_infinity.m) 
van_point_h = 1e03 * [0.6598; -1.4908; 0.0010];

% Compute Inverse Homography and Extract Directions
H_met_inv = inv(H_met);
h1 = H_met_inv(:, 1);
h2 = H_met_inv(:, 2);

% Normalize h1 and h2
h1 = h1 / norm(h1(1:2));
h2 = h2 / norm(h2(1:2));

% Define the Symbolic Variables and Constraints
syms x1 x2 x3 x4;
x = [
    x1 0 x2;
    0 1 x3;
    x2 x3 x4
];

% Constraints
eq1 = h1.' * x * h2 == 0; % Orthogonality
eq2 = h1.' * x * h1 - h2.' * x * h2 == 0; % Isotropy
eq3 = van_point_h.' * x * h1 == 0; % Vanishing Point Constraint 1
eq4 = van_point_h.' * x * h2 == 0; % Vanishing Point Constraint 2

% Solve the Equations
S = solve([eq1, eq2, eq3, eq4], [x1, x2, x3, x4], 'Real', true);

% Extract Solutions
s1 = double(S.x1);
s2 = double(S.x2);
s3 = double(S.x3);
s4 = double(S.x4);

% Construct Omega
omega = [
    s1 0 s2;
    0 1 s3;
    s2 s3 s4
];

fprintf('Constructed Omega:\n');
disp(omega);

% Compute Calibration Matrix K Correctly
try
    matrix_K = chol(inv(omega));
    fprintf('The calibration matrix K is:\n');
    matrix_K(1,2) = 0;
    matrix_K = matrix_K / matrix_K(3, 3); % Normalize K
    disp(matrix_K);
catch ME
    warning('Cholesky decomposition failed');
end
