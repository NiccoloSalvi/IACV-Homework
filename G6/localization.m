K = 1e+03 * [
    1.4896 0 0.0015;
    0 2.1845 -0.0002;
    0 0 0.0010
];

imagePoints = [
    300, 713;
    1405, 695;
    1191, 748;
    419, 777;
    397, 233;
    1354, 370;
];

worldPoints = [
    0, 0, 0;
    1, 0, 0;
    1, 1/3, 0;
    0, 1/3, 0;
    0, 0, 1/2;
    1, 0, 1/2;
];

[R, t] = localizeCamera_6Points(worldPoints, imagePoints, K);
centroid = mean(worldPoints, 1);  % Compute centroid of the parallelepiped
visualize_pose(R, t, centroid);

function visualize_pose(R, t, centroid)
    % Create a figure to visualize the camera pose and parallelepiped
    figure;
    hold on;
    grid on;
    axis equal;
    
    % Draw coordinate axes
    axisLength = 0.5;
    plot3([0 axisLength], [0 0], [0 0], 'r-', 'LineWidth', 2); % X-axis
    plot3([0 0], [0 axisLength], [0 0], 'g-', 'LineWidth', 2); % Y-axis
    plot3([0 0], [0 0], [0 axisLength], 'b-', 'LineWidth', 2); % Z-axis
    
    % Draw parallelepiped
    vertices = [
        0 0 0;
        1 0 0;
        1 1/3 0;
        0 1/3 0;
        0 0 1/2;
        1 0 1/2;
        1 1/3 1/2;
        0 1/3 1/2
    ];
    
    faces = [
        1 2 3 4;  % bottom
        5 6 7 8;  % top
        1 2 6 5;  % front
        2 3 7 6;  % right
        3 4 8 7;  % back
        4 1 5 8   % left
    ];
    
    patch('Vertices', vertices, 'Faces', faces, 'FaceColor', 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'k');
      
    % ---------------------------------------------------------------------
    % ENFORCE ORTHONORMALITY ON R:
    [U, ~, V] = svd(R);
    R_ortho = U * V';   % closest orthonormal matrix to R
    % If det < 0, fix the sign so we get a proper rotation (right-handed)
    if det(R_ortho) < 0
        R_ortho(:, 3) = -R_ortho(:, 3);
    end
    
    % 'Orientation' must be camera->world, so if your R is world->camera,
    % you pass R_ortho' to plotCamera. Usually for the standard equation:
    %    x_c = R * x_w + t,
    % the location in world coords is C = -R^T t 
    % and orientation is R^T if you want camera->world.
    % ---------------------------------------------------------------------
    cameraLocation = -R_ortho' * t;
    viewingDirection = (centroid - cameraLocation')';  % Vector from camera to centroid
    viewingDirection = viewingDirection / norm(viewingDirection);
    r3 = viewingDirection;
    r1 = cross([0; 0; 1], r3);  % Compute X-axis (arbitrary cross product to ensure orthogonality)
    r1 = r1 / norm(r1);  % Normalize
    r2 = cross(r3, r1);  % Compute Y-axis to complete the basis
    R_corrected = [r1, r2, r3];
    cameraOrientation = R_corrected';

    plotCamera('Location', cameraLocation, 'Orientation', cameraOrientation, 'Size', 0.15);
    
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('Camera Pose Relative to Parallelepiped');
    view(45, 30);
    hold off;
end


%   INPUTS:
%     worldPoints : Nx3  array with the 3D coordinates of each corner/feature 
%                    in the chosen world reference system (e.g., parallelepiped).
%     imagePoints : Nx2  array with the 2D coordinates [u_i, v_i] of each point 
%                    in the image (in pixels).
%     K           : 3x3  camera calibration matrix (zero-skew).
%
%   OUTPUTS:
%     R : 3x3 rotation matrix (world -> camera).
%     t : 3x1 translation vector (world -> camera).
%
%   STEPS:
%    1. We form a 3x4 camera matrix P via the DLT approach using the N >= 6 correspondences.
%    2. Factor out K:   P = K * [R | t].
%    3. From M = inv(K)*P, we get columns that correspond to R (up to scale) and t (up to scale).
%    4. Enforce that R is a proper rotation, and solve for scale.
function [R, t] = localizeCamera_6Points(worldPoints, imagePoints, K)
    %% Basic checks
    if size(worldPoints,1) < 6
        error('Need at least 6 non-coplanar points for a direct 3D DLT pose estimation.');
    end
    if size(worldPoints,1) ~= size(imagePoints,1)
        error('Number of 3D points must match number of 2D points.');
    end
    
    %% 1) Build the design matrix A for DLT
    %
    %   We want to solve  (u_i, v_i, 1)^T ~ P * (X_i, Y_i, Z_i, 1)^T,
    %   or in expanded form:
    %
    %      u_i = (p11 X_i + p12 Y_i + p13 Z_i + p14) / (p31 X_i + p32 Y_i + p33 Z_i + p34)
    %      v_i = (p21 X_i + p22 Y_i + p23 Z_i + p24) / (p31 X_i + p32 Y_i + p33 Z_i + p34)
    %
    %   This leads to the linear system in 12 unknowns (the entries of P, up to scale).
    %
    N = size(worldPoints,1);
    A = zeros(2*N, 12);
    for i = 1:N
        X = worldPoints(i,1);
        Y = worldPoints(i,2);
        Z = worldPoints(i,3);
        u = imagePoints(i,1);
        v = imagePoints(i,2);
        
        % Row for the 'u' equation
        A(2*i-1,:) = [ X, Y, Z, 1,  0, 0, 0, 0,  -u*X, -u*Y, -u*Z, -u ];
        % Row for the 'v' equation
        A(2*i,:)   = [ 0, 0, 0, 0,  X, Y, Z, 1,  -v*X, -v*Y, -v*Z, -v ];
    end
    
    %% 2) Solve A * p = 0 by SVD
    [~, ~, V] = svd(A, 'econ');
    p = V(:,end);         % last column is solution for homogeneous system
    P = reshape(p, 4, 3)';% put into 3x4 form. 
                          % Now P is only determined up to an overall scale factor.
    
    %% 3) Factor out the known intrinsics K:
    %      P = K * [R | t].
    %   Then:  M = inv(K) * P = [R | t]   (up to scale).
    M = K \ P;  % same as inv(K)*P
    
    % M is 3x4:  M = [m1 m2 m3 m4], where each m? is a 3x1 column.
    r1 = M(:,1);
    r2 = M(:,2);
    r3 = M(:,3);
    t_ = M(:,4);
    
    %% 4) Enforce that R = [r1 r2 r3] is orthonormal. 
    % We find a common scale factor from r1, r2, r3 norms. 
    % Typically we use the average of the first two column norms.
    s1 = norm(r1);
    s2 = norm(r2);
    s3 = norm(r3);  % you can average all three if you prefer
    s_avg = (s1 + s2 + s3)/3;  % or just (s1 + s2)/2 if r3 is too noisy
    
    r1 = r1 / s_avg;
    r2 = r2 / s_avg;
    % r3 = r3 / s_avg;
    t_ = t_ / s_avg;
    
    % Re-orthonormalize for numerical stability:
    % (1) ensure r2 is orthonormal to r1
    r2 = r2 - (dot(r1, r2))*r1; 
    r2 = r2 / norm(r2);
    % (2) set r3 = r1 x r2
    r3 = cross(r1, r2);
    r3 = r3 / norm(r3); % just to be safe
    
    R_ = [r1, r2, r3];
    
    % Check determinant. If negative => we must flip one axis
    if det(R_) < 0
        R_(:,3) = -R_(:,3);
    end
    
    %% 5) Return the final R, t
    R = R_;
    t = t_;

    disp('Camera pose estimated successfully.');
    disp('Rotation matrix:');
    disp(R);
    disp('Translation vector:');
    disp(t);
end

% Vertices of the parallelepiped (world coordinates)
vertices = [
    0, 0, 0;
    1, 0, 0;
    1, 1/3, 0;
    0, 1/3, 0;
    0, 0, 1/2;
    1, 0, 1/2;
    1, 1/3, 1/2;
    0, 1/3, 1/2
];

% Project each vertex to the image plane
for i = 1:size(vertices, 1)
    X_world = [vertices(i, :)'; 1];  % Homogeneous coordinates
    X_camera = R * X_world(1:3) + t;  % Camera coordinates
    X_image = K * X_camera;  % Image plane projection
    X_image = X_image / X_image(3);  % Normalize
    
    fprintf('Vertex %d projects to (u, v) = (%.2f, %.2f)\n', i, X_image(1), X_image(2));
    
    % Check if it falls within valid image bounds (e.g., [0, width] x [0, height])
    if X_image(1) < 0 || X_image(2) < 0 || isnan(X_image(1)) || isnan(X_image(2))
        fprintf('    -> Outside the image bounds or not visible!\n');
    else
        fprintf('    -> Visible in the image!\n');
    end
end
