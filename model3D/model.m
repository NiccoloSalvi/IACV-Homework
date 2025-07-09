% Dimensions of the overall outer box
width  = 1.0; % X dimension
height = 1.0/2; % Y dimension
depth  = 1.0/3; % Z dimension

% Thickness of each "board"
T = 0.02;

% Positions of the two vertical dividers (fraction of shelf width)
spacer1_x = width / 3;
spacer2_x = 2 * width / 3;

% Create a figure
figure('Name', 'Hollow Shelf with Thick Edges and Dividers');
hold on; grid on;
axis equal; view(3);
xlabel('X'); ylabel('Y'); zlabel('Z');

% Set specific viewing angle to match your image
% view(az, el) where:
% az = azimuth (rotation around z-axis)
% el = elevation (angle above xy-plane)
view(0, -80);  % You can adjust these values to match your desired orientation

% -----------
% (1) BACK BOARD
% -----------
% A thin plank at the rear, covering the entire back
% from Z = depth-T to Z = depth
[Vback, Fback] = rectPrism( ...
    0,         width, ...    % X range
    0,         height, ...   % Y range
    depth - T, depth );      % Z range
patch('Vertices', Vback, 'Faces', Fback, ...
        'FaceColor', [0.3 0.3 0.3], ...  % Dark gray for the back
        'EdgeColor', 'k');

% -----------
% (2) LEFT SIDE BOARD
% -----------
% A vertical plank on the left, from x=0 to x=T, full height & depth
[Vleft, Fleft] = rectPrism( ...
    0,    T, ...            % X range
    0,    height, ...       % Y range
    0,    depth );          % Z range
patch('Vertices', Vleft, 'Faces', Fleft, ...
        'FaceColor', [0.7 0.9 0.2], ...  % Light green
        'EdgeColor', 'k');
    
% -----------
% (3) RIGHT SIDE BOARD
% -----------
[Vright, Fright] = rectPrism( ...
    width - T,  width, ...  % X range
    0,          height, ...
    0,          depth );
patch('Vertices', Vright, 'Faces', Fright, ...
        'FaceColor', [0.7 0.9 0.2], ...
        'EdgeColor', 'k');
    
% -----------
% (4) BOTTOM BOARD
% -----------
% from y=0 to y=T
[Vbottom, Fbottom] = rectPrism( ...
    0,      width, ...
    0,      T, ...
    0,      depth );
patch('Vertices', Vbottom, 'Faces', Fbottom, ...
        'FaceColor', [0.7 0.9 0.2], ...
        'EdgeColor', 'k');
    
% -----------
% (5) TOP BOARD
% -----------
% from y=height - T to y=height
[Vtop, Ftop] = rectPrism( ...
    0,        width, ...
    height - T,  height, ...
    0,        depth );
patch('Vertices', Vtop, 'Faces', Ftop, ...
        'FaceColor', [0.7 0.9 0.2], ...
        'EdgeColor', 'k');
    
% -----------
% (6) FIRST VERTICAL DIVIDER BOARD
% -----------
% Center the divider at x = spacer1_x, with thickness T in x
% so x = [spacer1_x - T/2, spacer1_x + T/2]
[Vdiv1, Fdiv1] = rectPrism( ...
    spacer1_x - T/2, spacer1_x + T/2, ...
    0,               height, ...
    0,               depth );
patch('Vertices', Vdiv1, 'Faces', Fdiv1, ...
        'FaceColor', [0.7 0.9 0.2], ...
        'EdgeColor', 'k');
    
% -----------
% (7) SECOND VERTICAL DIVIDER BOARD
% -----------
[Vdiv2, Fdiv2] = rectPrism( ...
    spacer2_x - T/2, spacer2_x + T/2, ...
    0,               height, ...
    0,               depth );
patch('Vertices', Vdiv2, 'Faces', Fdiv2, ...
        'FaceColor', [0.7 0.9 0.2], ...
        'EdgeColor', 'k');
    
% Clean up and show
title('Hollow Shelf with Back Board and Thick Edges/Dividers');
hold off;

%% --- Helper Function: Rectangular Prism vertices & faces
function [V,F] = rectPrism(xmin, xmax, ymin, ymax, zmin, zmax)
% Returns a list of 8 vertices (V) and 6 faces (F) for a rectangular prism.
% You can then feed these to 'patch("Vertices", V, "Faces", F, ...)'.

    V = [ ...
        xmin, ymin, zmin;  % 1
        xmax, ymin, zmin;  % 2
        xmax, ymax, zmin;  % 3
        xmin, ymax, zmin;  % 4
        xmin, ymin, zmax;  % 5
        xmax, ymin, zmax;  % 6
        xmax, ymax, zmax;  % 7
        xmin, ymax, zmax;  % 8
    ];
    F = [ ...
        1 2 3 4;  % bottom face  (zmin)
        5 6 7 8;  % top face     (zmax)
        1 2 6 5;  % front face   (ymin side)
        2 3 7 6;  % right face   (xmax side)
        3 4 8 7;  % back face    (ymax side)
        4 1 5 8;  % left face    (xmin side)
    ];
end
